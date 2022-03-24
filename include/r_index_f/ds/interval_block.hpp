/* interval_block - Compresses intervals into a block using SDSL/template
    Copyright (C) 2021 Nathaniel Brown
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file interval_block.hpp
   \brief interval_block.hpp Compact representation of intervals as a block
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 09/07/2020
*/

#ifndef _I_BLOCK_HH
#define _I_BLOCK_HH

#include <common.hpp>
#include <ds/interval_pos.hpp>
#include <ds/heads_wt_w.hpp>
#include <ds/intervals_rank_w.hpp>
#include <ds/symbol_map.hpp>
#include <sdsl/dac_vector.hpp>

#include <ds/heads_bv_w.hpp>
#include <ds/ACGT_map.hpp>

using namespace sdsl;

template  < class heads_t = heads_bv_w<>,
            class intervals_t = intervals_rank_w<>,
            class lengths_t = dac_vector_dp<>,
            class offsets_t = dac_vector_dp<>,
            template<class> class char_map_t = ACGT_map >
class interval_block
{
private:
    typedef char_map_t<interval_pos> pos_map;

    // Heads supporting rank/select/access for character c
    heads_t heads;
    // Intervals which return the mapping given a character c and its rank in block
    intervals_t intervals;
    // Lengths supporting access
    lengths_t lengths;
    // Offsets supporting access
    offsets_t offsets;

    // Stores prior and next block LF mapping for a character (for block overruns)
    pos_map prior_block_LF;
    pos_map next_block_LF;

    // For a row k with offset d, character c of rank c_rank, compute it's LF
    interval_pos LF(ulint k, ulint d, uchar c, ulint c_rank)
    {
        ulint q = get_interval(c, c_rank);
        ulint d_prime = d + get_offset(k);

        return interval_pos(q, d_prime);
    }

public:
    interval_block() {}

    interval_block(std::vector<uchar> chars, std::vector<ulint> ints, std::vector<ulint> lens, std::vector<ulint> offs, std::unordered_map<uchar, interval_pos> prior_LF) {
        heads = heads_t(chars);
        intervals = intervals_t(chars, ints);
        lengths = lengths_t(lens);
        offsets = offsets_t(offs);

        prior_block_LF = pos_map(prior_LF);
        next_block_LF = pos_map();
    }

    // Return the character at row k
    const ulint get_char(const ulint k)
    {
        return heads[k];
    }

    // For a given character and it's rank, return the interval mapping (LF)
    const ulint get_interval(const uchar c, const ulint c_rank)
    {
        return intervals.get(c, c_rank);
    }

    // Get the length at row k
    const ulint get_length(const ulint k)
    {
        return lengths[k];
    }

    // Get the offset at row k
    const ulint get_offset(const ulint k)
    {
        return offsets[k];
    }

    bool has_prior_LF(const uchar c)
    {
        return prior_block_LF.contains(c);
    }

    bool has_next_LF(const uchar c)
    {
        return next_block_LF.contains(c);
    }

    void set_next_LF(const uchar c, interval_pos next_LF)
    {
        next_block_LF.insert(std::pair<uchar, interval_pos>(c, next_LF));
    }

    // For row k wih offset d, compute the LF mapping
    interval_pos LF(const ulint k, const ulint d)
    {
        const auto [c_rank, c] = heads.inverse_select(k);
        return LF(k, d, c, c_rank);
    }

    // Perform the LF mapping for character c prior or at position k with offset d
    interval_pos LF_prior(const ulint k, const ulint d, const uchar c)
    {
        // Look in row ahead so that the rank includes current position
        ulint c_rank = heads.rank(k + 1, c);
        // If there are no c prior to position in block, return LF of prior c in another block
        if (c_rank == 0) 
        {
            if (has_prior_LF(c))
            {
                return prior_block_LF[c];
            }
            else
            {
                return interval_pos();
            }
        }
        // We subtract 1 to maintain 0-based rank after ensuring it is not 0, since we use unsigned values
        else
        {
            c_rank -= 1;
        }

        ulint k_prime = heads.select(c_rank + 1, c);
        // If our k changed, set the offset to the last character in that prior run
        ulint d_prime = (k != k_prime) ? lengths[k_prime] - 1 : d;

        return LF(k_prime, d_prime, c, c_rank);
    }

    // Perform the LF mapping for character c succeding or at position k with offset d
    interval_pos LF_next(const ulint k, const ulint d, const uchar c)
    {
        // Count occ of c before position
        ulint c_rank = heads.rank(k, c);
        // If the c of rank at or succeding our position overruns the block, return LF of next c in another block
        if (c_rank + 1 > heads.rank(heads.size(), c))
        {
            if (has_next_LF(c))
            {
                return next_block_LF[c];
            }
            else
            {
                return interval_pos();
            }
        }
        ulint k_prime = heads.select(c_rank + 1, c);
        // If k changed, set it to the first character of the next run
        ulint d_prime = (k != k_prime) ? 0 : d;

        return LF(k_prime, d_prime, c, c_rank);
    }

    // Reduces position until offset shorter than length of interval, or returns if at end of block
    interval_pos reduce(interval_pos pos, ulint k)
    {
        ulint q = pos.run;
        ulint d = pos.offset;
        ulint next_len;
	    while (k < lengths.size() && d >= (next_len = get_length(k))) 
        {
            d -= next_len;
            ++k;
            ++q;
        }

	    return interval_pos(q, d);
    }

    /* serialize the interval block to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;


        written_bytes += heads.serialize(out,v,"Heads");
        written_bytes += intervals.serialize(out,v,"Intervals");
        written_bytes += lengths.serialize(out,v,"Lengths");
        written_bytes += offsets.serialize(out,v,"Offsets");

        written_bytes += prior_block_LF.serialize(out,v,"Prior_Block_LF");
        written_bytes += next_block_LF.serialize(out,v,"Next_Block_LF");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the interval block from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        heads.load(in);
        intervals.load(in);
        lengths.load(in);
        offsets.load(in);

        prior_block_LF.load(in);
        next_block_LF.load(in);
    }
};

#endif /* end of include guard: _I_BLOCK_HH */
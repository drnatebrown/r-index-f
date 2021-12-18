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

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/dac_vector.hpp>

using namespace sdsl;

template  < class wt_t = wt_huff<bit_vector>,
            class bit_vec = bit_vector,
            class dac_vec = dac_vector<> >
class interval_block
{
private:
    typedef typename bit_vec::select_1_type bv_select_1;

    // O(ALPHABET_SIZE)-space bits determining if a character was present in this block, used for serialization/error checks
    bit_vec char_present;

    // Keep Wavelet Tree for now
    // Interval Heads stored in Wavelet Tree
    wt_t heads;

    // Vectors indexed by byte (ASCII) to compute interval mapping
    std::vector<ulint> char_base_interval;
    std::vector<bv_select_1> char_diff_select;
    std::vector<bit_vec> char_diff_vec;

    // Compact Variable Length DAC for Length/Offset
    dac_vec lengths;
    dac_vec offsets;

    // Stores prior and next block LF mapping for a character (for block overruns)
    bit_vec char_prior;
    std::vector<interval_pos> prior_block_LF;
    bit_vec char_next;
    std::vector<interval_pos> next_block_LF;

    // For a row k with offset d, character c of rank c_rank, compute it's LF
    interval_pos LF(ulint k, ulint d, char c, ulint c_rank)
    {
        ulint q = get_interval(c, c_rank);
        ulint d_prime = d + offsets[k];

        return interval_pos(q, d_prime);
    }

public:
    interval_block() {}

    // Simple constructor for block, work to compute values done externally (i.e. no logic enforced here)
    interval_block(std::vector<char> chars, std::vector<ulint> base_map, std::vector<vector<bool>> diff_vec,
                    std::vector<ulint> lens, std::vector<ulint> offs, std::vector<bool> has_prior, std::vector<interval_pos> prior_LF) {
        assert(base_map.size() == ALPHABET_SIZE);
        assert(diff_vec.size() == ALPHABET_SIZE);
        assert(has_prior.size() == ALPHABET_SIZE);
        assert(prior_LF.size() == ALPHABET_SIZE);

        char_present = bit_vec(ALPHABET_SIZE, false);

        construct_im(heads, std::string(chars.begin(), chars.end()).c_str(), 1);

        char_base_interval = base_map;

        char_diff_vec = std::vector<bit_vec>(ALPHABET_SIZE);
        char_diff_select = std::vector<bv_select_1>(ALPHABET_SIZE);
        for(size_t i = 0; i < ALPHABET_SIZE; i++)
        {
            if(!diff_vec[i].empty())
            {
                char_diff_vec[i] = bool_to_bit_vec<bit_vec>(diff_vec[i]);
                char_diff_select[i] = bv_select_1(&char_diff_vec[i]);

                char_present[i] = true;
            }
        }

        lengths = dac_vec(lens);
        offsets = dac_vec(offs);

        char_prior = bool_to_bit_vec<bit_vec>(has_prior);
        prior_block_LF = prior_LF;
        char_next = bit_vec(ALPHABET_SIZE, false);
        next_block_LF = std::vector<interval_pos>(ALPHABET_SIZE, interval_pos());
    }

    // Return the character at row k
    const ulint get_char(const ulint k)
    {
        return heads[k];
    }

    // For a given character and it's rank, return the interval mapping (LF)
    const ulint get_interval(const char c, const ulint c_rank)
    {
        return char_base_interval[c] + char_diff_select[c](c_rank+1) - c_rank;
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

    bool has_prior_LF(const char c)
    {
        return char_prior[c];
    }

    bool has_next_LF(const char c)
    {
        return char_next[c];
    }

    void set_next_LF(const char c, interval_pos next_LF)
    {
        char_next[c] = true;
        next_block_LF[c] = next_LF;
    }

    // For row k wih offset d, compute the LF mapping
    interval_pos LF(const ulint k, const ulint d)
    {
        const auto [c_rank, c] = heads.inverse_select(k);
        return LF(k, d, c, c_rank);
    }

    // Perform the LF mapping for character c prior or at position k with offset d
    interval_pos LF_prior(const ulint k, const ulint d, const char c)
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
    interval_pos LF_next(const ulint k, const ulint d, const char c)
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
	    while (k < lengths.size() && d >= (next_len = lengths[k])) 
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

        // Serliaze character present bit vector
        char_present.serialize(out, v, "char_present");

        // Serialize wavelet tree (heads)
        written_bytes += heads.serialize(out,v,"Heads");

        for (size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(char_present[i])
            {
                // Serialize base map
                out.write((char*) &char_base_interval[i], sizeof(char_base_interval[i]));
                written_bytes += sizeof(char_base_interval[i]); 

                // Serialize diff vec
                written_bytes += char_diff_vec[i].serialize(out,v,"char_diff_vec_" + std::to_string(i)); 
            }
        }

        // Serialize DACs (Length/Offset)
        written_bytes += lengths.serialize(out,v,"Lengths");
        written_bytes += offsets.serialize(out,v,"Offsets");

        char_prior.serialize(out, v, "char_prior");
        // Serialize prior char LF (prior block)
        for (size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(char_prior[i])
            {
                written_bytes += prior_block_LF[i].serialize(out, v, "prior_LF_" + std::to_string(i));
            }
        }

        // Serialize next char LF (next block)
        char_next.serialize(out, v, "char_next");
        for (size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(char_next[i])
            {
                written_bytes += next_block_LF[i].serialize(out, v, "next_LF_" + std::to_string(i)); 
            }
        }

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the interval block from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        // Load char present bitvector
        char_present.load(in);

        // Load wavelet tree (heads)
        heads.load(in);

        char_base_interval = std::vector<ulint>(ALPHABET_SIZE);
        char_diff_vec = std::vector<bit_vec>(ALPHABET_SIZE);
        char_diff_select = std::vector<bv_select_1>(ALPHABET_SIZE);
        for(size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(char_present[i])
            {
                // Read base map
                in.read((char *)&char_base_interval[i], sizeof(char_base_interval[i]));

                // Load diff vec and select support
                char_diff_vec[i].load(in);
                char_diff_select[i] = bv_select_1(&char_diff_vec[i]);
            }
        }

        // Load DACs (Length/Offset)
        lengths.load(in);
        offsets.load(in);

        // Load prior char LF (prior block)
        char_prior.load(in);
        prior_block_LF = std::vector<interval_pos>(ALPHABET_SIZE);
        for(size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if (char_prior[i])
            {
                prior_block_LF[i].load(in);
            }
        }

        // Load next char LF (next block)
        char_next.load(in);
        next_block_LF = std::vector<interval_pos>(ALPHABET_SIZE);
        for(size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(char_next[i])
            {
                next_block_LF[i].load(in);
            }
        }
    }
};

#endif /* end of include guard: _I_BLOCK_HH */
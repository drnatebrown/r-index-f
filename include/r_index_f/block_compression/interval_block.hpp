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
    typedef bit_vector::select_1_type bv_select_1;

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

    // Index for start of block wrt. BWT
    ulint idx;

    // Stores prior and next block LF mapping for a character (for block overruns)
    std::vector<interval_pos> prior_block_LF;
    std::vector<interval_pos> next_block_LF;

    // For a row k with offset d, character c of rank c_rank, compute it's LF
    interval_pos LF(ulint k, ulint d, char c, ulint c_rank)
    {
        ulint q = get_interval(c, c_rank);
        ulint d_prime = d + offsets[k];

        return interval_pos(q, d_prime)
    }

    // Convert boolean vector to specified bit vector
    bit_vec bool_to_bit_vec(vector<bool> &b)
    {
		if(b.size()==0) return bit_vec();

		bit_vector bv(b.size());

		for(size_t i = 0; i < b.size(); ++i)
			bv[i] = b[i];

		return bit_vec(bv);
	}

public:
    interval_block() {}

    // Simple constructor for block, work to compute values done externally (i.e. no logic enforced here)
    interval_block(std::vector<ulint> chars, std::vector<ulint> base_map, std::vector<bit_vec> diff_vec, std::vector<bv_select_1> 
                    std::vector<ulint> lens, std::vector<ulint> offs, ulint i, std::vector<interval_pos> prior_LF, std::vector<interval_pos> next_LF) {
        
        construct_im(heads, std::string(chars.begin(), chars.end()).c_str(), 1);

        char_base_interval = base_map;
        char_diff_vec = diff_vec;
        char_diff_select = bv_select_1(&)

        lengths = dac_vec(lens);
        offsets = dac_vec(offsets);
                
        idx = i;

        prior_block_LF = prior_LF;
        next_block_LF = next_LF;
    }

    // Return the character at row k
    const ulint get_char(const ulint k)
    {
        return heads[k];
    }

    // For a given character and it's rank, return the interval mapping (LF)
    const ulint get_interval(const char c, const ulint c_rank)
    {
        return char_base_interval[c] + char_diff_vec[c](c_rank+1) - c_rank;
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
            return prior_block_LF[c];
        }
        k_prime = heads.select(c_rank);
        // If our k changed, set the offset to the last character in that prior run
        d_prime = (k != k_prime) ? lengths[k_prime] - 1 : d;

        return LF(k_prime, d_prime, c, c_rank);
    }

    // Perform the LF mapping for character c succeding or at position k with offset d
    interval_pos LF_next(const ulint k, const ulint d, const char c)
    {
        // Count occ of c before position
        ulint c_rank = heads.rank(k, c);
        // If the c of rank at or succeding our position overruns the block, return LF of next c in another block
        if (c_rank + 1 > heads.rank(k + 1, c))
        {
            return next_block_LF[c];
        }
        k_prime = heads.select(c_rank + 1);
        // If k changed, set it to the first character of the next run
        d_prime = (k != k_prime) ? 0 : d;

        return LF(k_prime, d_prime, c, c_rank);
    }

    /* serialize the interval block to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        // Serialize wavelet tree (heads)
        written_bytes += heads.serialize(out,v,"Heads");

        // Serialize base intervals
        size_t size = char_base_interval.size();
        out.write((char*) &size, sizeof(size_t));
        written_bytes += sizeof(size);
        for (const auto& val: char_base_interval)
        {
            out.write((char*) &val, sizeof(val));
            written_bytes += sizeof(val); 
        }

        // Serialize diff vector
        size_t size = char_diff_vec.size();
        out.write((char*) &size, sizeof(size_t));
        written_bytes += sizeof(size);
        char curr_c = 0;
        for (const auto& val: char_diff_vec)
        {
            written_bytes += val.serialize(out,v,"char_diff_vec_" + std::to_string(curr_c++)); 
        }

        // Serialize DACs (Length/Offset)
        written_bytes += lengths.serialize(out,v,"Lengths");
        written_bytes += offsets.serialize(out,v,"Offsets");

        // Serialize block index wrt. BWT
        out.write((char *)&idx, sizeof(idx));
        written_bytes += sizeof(idx);

        // Serialize prior char LF (prior block)
        size_t size = prior_block_LF.size();
        out.write((char*) &size, sizeof(size_t));
        written_bytes += sizeof(size);
        for (const auto& val: prior_block_LF)
        {
            out.write((char*) &val, sizeof(val));
            written_bytes += sizeof(val); 
        }

        // Serialize next char LF (next block)
        size_t size = next_block_LF.size();
        out.write((char*) &size, sizeof(size_t));
        written_bytes += sizeof(size);
        for (const auto& val: next_block_LF)
        {
            out.write((char*) &val, sizeof(val));
            written_bytes += sizeof(val); 
        }

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the interval block from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        // Load wavelet tree (heads)
        heads.load(in);

        // Load base interval mapping
        size_t size;
        in.read((char *)&size, sizeof(size));
        for(size_t i = 0; i < size; ++i)
        {
            ulint val;
            in.read((char *)&val, sizeof(val));
            char_base_interval[key] = val;
        }

        // Load char diff vectors
        in.read((char *)&size, sizeof(size));
        for(size_t i = 0; i < size; ++i)
        {
            bit_vec val;
            val.load(in);
            char_diff_vec[i] = val;
            char_diff_select[i] = bv_select_1(&char_diff_vec[i]);
        }

        // Load DACs (Length/Offset)
        lengths.load(in);
        offsets.load(in);

        // Load block index wrt. BWT
        in.read((char *)&idx, sizeof(idx));

        // Load prior char LF (prior block)
        in.read((char *)&size, sizeof(size));
        for(size_t i = 0; i < size; ++i)
        {
            iterval_pos val;
            in.read((char *)&val, sizeof(val));
            val.load(in);
            prior_block_LF[i] = val;
        }

        // Load next char LF (next block)
        in.read((char *)&size, sizeof(size));
        for(size_t i = 0; i < size; ++i)
        {
            interval_pos val;
            in.read((char *)&val, sizeof(val));
            val.load(in);
            next_block_LF[i] = val;
        }
    }
};

#endif /* end of include guard: _I_BLOCK_HH */
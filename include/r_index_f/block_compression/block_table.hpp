/* interval_block - Wrapper for vector holding interval blocks
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
   \file block_table.hpp
   \brief block_table.hpp Wrapper for vector holding interval blocks
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 19/11/2021
*/

#ifndef _BLOCK_TABLE_HH
#define _BLOCK_TABLE_HH

#include <common.hpp>
#include <interval_block.hpp>
#include <interval_pos.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/dac_vector.hpp>

using namespace sdsl;

template  < ulint block_size = 65536, // 2^16
            class wt_t = wt_huff<bit_vector>,
            class bit_vec = bit_vector,
            class dac_vec = dac_vector<> >
class block_table
{
private:
    vector<interval_block<block_size, wt_t, bit_vec, dac_vec>> blocks;
    ulint r;
    ulint n;

public:
    block_table() {}

    block_table(const LF_table LF_rows)
    {
        r = LF_rows.runs();
        n = LF.rows.size();

        // Round up if quotient not whole
        ulint B_len = (r / block_size) + ((r % block_size) != 0);
        blocks = vector<interval_block>(B_len);

        vector<char> block_chars = vector<char>();
        vector<ulint> block_intervals = vector<ulint>();
        vector<ulint> block_lens = vector<ulint>();
        vector<ulint> block_offsets = vector<ulint>();

        ulint block_idx = 0;
        ulint next_idx = 0;

        ulint b = 0;
        i = 0;
        while (i < r) 
        {
            LF_row curr = LF_rows.get(i);

            block_chars.push_back(curr.character);
            block_intervals.push_back(curr.interval);
            block_lens.push_back(curr.length);
            block_offsets.push_back(curr.offset);

            next_idx += l;

            if (!seen_c.count(curr.character)) {
                block_c_map.insert(std::pair<char, ulint>(curr.character, curr.interval));
                last_c_map.insert(std::pair<char, ulint>(curr.character, curr.interval));
                bit_diff[curr.character] = (std::pair<char, vector<bool>>(c, vector<bool>()));

                if (b > 0)
                {
                    interval_pos next_c;

                    ulint c_b = k;
                    ulint c_off = d;

                    while (c_off >= lens[c_b]) 
                    {
                        c_off -= lens[c_b];
                        ++c_b;
                    }
                    next_c = {c_b, c_off};

                    ulint b_curr = b;
                    while (b_curr > 0 && !next_valid[b_curr-1][c])
                    {
                        B_table[b_curr-1].else_next_LF.insert(std::pair<char, i_position>(c, next_c));
                        next_valid[b_curr-1][c] = true;
                        --b_curr;
                    }
                }
            }

            ulint diff = k - last_c_map[c];
            while (diff > 0) {
                bit_diff[c].push_back(false);
                --diff;
            }
            bit_diff[c].push_back(true);

            last_c_map[c] = k;

            ++i;

            // End of block of intervals, update block table
            if (i % block_size == 0 || i >= r)
            {
                i_block& curr = B_table[b];

                for (auto& kv: prior_c_map)
                {
                    if (b > 0)
                    {
                        i_position prior_c;
                        prior_valid[c] = true;

                        char c = kv.first;
                        ulint c_pos = kv.second;

                        ulint c_b = intervals[c_pos];
                        ulint c_off = (lens[c_pos] - 1) + offsets[c_pos];

                        while (c_off >= lens[c_b]) 
                        {
                            c_off -= lens[c_b];
                            ++c_b;
                        }

                        prior_c = i_position{c_b, c_off};
                    }
                }
                        
                

                curr.lengths = dac_vec(block_lens);
                curr.offsets = dac_vec(block_offsets);
                
                curr.idx = block_idx;
                
                curr.next_is_valid = bv(next_valid[b]);
                curr.prior_is_valid = bv(prior_valid);

                block_chars = vector<char>(block_size);
                block_lens = vector<ulint>(block_size);
                block_offsets = vector<ulint>(block_size);
                block_c_map = std::unordered_map<char, ulint>();
                for(auto& kv : last_c_map)
                {
                    prior_c_map[kv.first] = kv.second;
                }
                last_c_map = std::unordered_map<char, ulint>();
                bit_diff = std::unordered_map<char, vector<bool>>();
                block_idx = next_idx;
                prior_valid = vector<bool>(256, false);

                ++b;
            }
        }

        return B_table;
    }

    interval_block& get_block(interval_pos)
    {
        return blocks[pos.run / block_size]
    }

    char get_char(ulint run) 
    {
        return (char) blocks[run / block_size].get_char(run%block_size);
    }

    char get_char(i_position pos) 
    {
        return get_char(pos.run);
    }

    interval_pos begin()
    {
        return interval_pos(0, 0);
    }

    interval_pos end()
    {
        return interval_pos(r-1), blocks[(r-1)/block_size].get_length((r-1)%block_size]-1));
    }

    ulint i_position_to_idx(interval_pos pos)
    {
        ulint pos_b = pos.run/block_size;
        ulint pos_k = pos.run%block_size;

        i_block& block = blocks[pos_b];
        
        ulint idx;
        ulint k;
        idx = block.get_idx();
        k = 0;
        while (k < pos_k)
        {
            idx += block.get_length(k);
            ++k;
        }
        idx += pos.offset;
        return idx;
    }

    /* serialize the interval block to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
    {
        out.write((char *)&n, sizeof(n));
        written_bytes += sizeof(n);

        out.write((char *)&r, sizeof(r));
        written_bytes += sizeof(r);

        size_t size = B_table.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        for(size_t i = 0; i < size; ++i)
        {
            written_bytes += table[i].serialize(out,v,"block_table_" + std::to_string(i));
        }

        return written_bytes;
    }

    /* load the interval block from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        size_t size;

        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));

        in.read((char *)&size, sizeof(size));
        table = std::vector<interval_block>(size);
        for(size_t i = 0; i < size; ++i)
        {
            B_table[i].load(in);
        }
    }
};

#endif /* end of include guard: _BLOCK_TABLE_HH */
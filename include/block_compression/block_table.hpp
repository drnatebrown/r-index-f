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
#include <LF_table.hpp>
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
    typedef interval_block<block_size, wt_t, bit_vec, dac_vec> block;

    vector<block> blocks;
    ulint r;
    ulint n;

    ulint get_length(ulint run)
    {
        return get_block(run).get_length(row(run));
    }
    
    ulint get_length(interval_pos pos)
    {
        return get_length(pos.run);
    }

    ulint row(ulint run)
    {
        return run % block_size;
    }

    ulint row(interval_pos pos)
    {
        return row(pos.run);
    }

public:
    block_table() {}

    block_table(LF_table LF_rows)
    {
        r = LF_rows.runs();
        n = LF_rows.size();

        // Round up if quotient not whole
        ulint B_len = (r / block_size) + ((r % block_size) != 0);
        blocks = vector<block>(B_len);

        vector<char> block_chars = vector<char>();
        vector<ulint> block_intervals = vector<ulint>();
        vector<ulint> block_lens = vector<ulint>();
        vector<ulint> block_offsets = vector<ulint>();

        // Concerned with Interval sectioned by character (break into base pointer and difference from prior interval)
        std::unordered_map<char, ulint> last_c_pos = std::unordered_map<char, ulint>();
        std::vector<ulint> block_c_map = std::vector<ulint>(ALPHABET_SIZE, 0);
        std::vector<vector<bool>> bit_diff = std::vector<vector<bool>>(ALPHABET_SIZE, vector<bool>());
        std::vector<interval_pos> prior_LF = std::vector<interval_pos>(ALPHABET_SIZE, interval_pos());

        ulint block_idx = 0;
        ulint next_idx = 0;

        ulint b = 0;
        ulint i = 0;
        while (i < r) 
        {
            LF_table::LF_row curr = LF_rows.get(i);

            block_chars.push_back(curr.character);
            block_intervals.push_back(curr.interval);
            block_lens.push_back(curr.length);
            block_offsets.push_back(curr.offset);

            next_idx += curr.length;

            if (!last_c_pos.count(curr.character)) {
                last_c_pos.insert(std::pair<char, ulint>(curr.character, block_chars.size() - 1));

                block_c_map[curr.character] = curr.interval;
                bit_diff[curr.character] = vector<bool>();

                // For all blocks prior without set values to find the next c's mapping, loop back and set
                if (b > 0)
                {
                    auto[k, d] = LF_rows.LF(curr.interval, curr.offset);
                    interval_pos next_c = interval_pos(k, d);

                    ulint b_curr = b;
                    while (b_curr > 0 && !blocks[b_curr-1].get_next_LF(curr.character).is_set())
                    {
                        blocks[b_curr-1].set_next_LF(curr.character, next_c);
                        --b_curr;
                    }
                }
            }

            ulint diff = curr.interval - block_intervals[last_c_pos[curr.character]];
            while (diff > 0) {
                bit_diff[curr.character].push_back(false);
                --diff;
            }
            bit_diff[curr.character].push_back(true);

            last_c_pos[curr.character] = block_chars.size() - 1;

            ++i;

            // End of block of intervals, update block table
            if (i % block_size == 0 || i >= r)
            {        
                blocks[b] = block(block_chars, block_c_map, bit_diff, block_lens, block_offsets, block_idx, prior_LF);

                for(auto const& [c, pos] : last_c_pos)
                {
                    // Perform LF step from the last seen character in this run (which is at offset equal to last character, one minus length)
                    auto[k, d] = LF_rows.LF(block_intervals[pos], block_lens[pos] - 1);
                    prior_LF[c] = interval_pos(k, d);
                }

                block_chars = vector<char>();
                block_intervals = vector<ulint>();
                block_lens = vector<ulint>();
                block_offsets = vector<ulint>();

                last_c_pos = std::unordered_map<char, ulint>();
                block_c_map = std::vector<ulint>(ALPHABET_SIZE, 0);
                bit_diff = std::vector<vector<bool>>(ALPHABET_SIZE, std::vector<bool>());
                block_idx = next_idx;

                ++b;
            }
        }
    }

    block& get_block(ulint run)
    {
        assert(run < r);
        return blocks[run / block_size];
    }

    block& get_block(interval_pos pos)
    {
        return get_block(pos.run);
    }

    char get_char(ulint run) 
    {
        return (char) get_block(run).get_char(row(run));
    }

    char get_char(interval_pos pos) 
    {
        return get_char(pos.run);
    }

    ulint runs()
    {
        return r;
    }

    ulint size()
    {
        return n;
    }

    interval_pos LF(interval_pos pos)
    {
        return reduced_pos(get_block(pos).LF(row(pos), pos.offset));
    }

    interval_pos LF_prior(interval_pos pos, char c)
    {
        return reduced_pos(get_block(pos).LF_prior(row(pos), pos.offset, c));
    }

    interval_pos LF_next(interval_pos pos, char c)
    {
        return reduced_pos(get_block(pos).LF_next(row(pos), pos.offset, c));
    }

    interval_pos reduced_pos(interval_pos pos)
    {
        interval_pos curr = pos;
        while (curr.offset >= get_length(curr))
        {
            curr = get_block(curr).reduce(curr);
        }

        return curr;
    }

    interval_pos begin()
    {
        return interval_pos(0, 0);
    }

    interval_pos end()
    {
        return interval_pos(r-1, get_length(r-1)-1);
    }

    ulint i_position_to_idx(interval_pos pos)
    {
        ulint pos_k = row(pos);
        block& b = get_block(pos);
        
        ulint idx;
        ulint k;
        idx = b.get_idx();
        k = 0;
        while (k < pos_k)
        {
            idx += b.get_length(k);
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
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        out.write((char *)&n, sizeof(n));
        written_bytes += sizeof(n);

        out.write((char *)&r, sizeof(r));
        written_bytes += sizeof(r);

        size_t size = blocks.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        for(size_t i = 0; i < size; ++i)
        {
            written_bytes += blocks[i].serialize(out,v,"block_table_" + std::to_string(i));
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
        blocks = std::vector<block>(size);
        for(size_t i = 0; i < size; ++i)
        {
            blocks[i].load(in);
        }
    }
};

#endif /* end of include guard: _BLOCK_TABLE_HH */
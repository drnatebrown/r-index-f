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

#include <ds/interval_block.hpp>
#include <ds/interval_pos.hpp>
#include <ds/idx_bit_vector.hpp>

using namespace sdsl;

template  < ulint block_size = 1048576, // 2^20
            class idx_vec = idx_bit_vector<>,
            class block = interval_block<>>
class block_table
{
private:

    vector<block> blocks;
    idx_vec block_idx;

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

        std::vector<uchar> block_chars = std::vector<uchar>();
        std::vector<ulint> block_intervals = std::vector<ulint>();
        std::vector<ulint> block_lens = std::vector<ulint>();
        std::vector<ulint> block_offsets = std::vector<ulint>();

        // Where characters prior to block mapped to, in case we can't find that character when we LF
        std::unordered_map<uchar, interval_pos> prior_LF = std::unordered_map<uchar, interval_pos>();

        // True marks where blocks start in terms of absolute idx
        std::vector<bool> block_base_idx = std::vector<bool>();
        // idx of current block
        ulint curr_idx = 0;
        // absolute idx counter
        ulint absolute_idx = 0;

        // Where the last character's position was wrt. current block
        std::unordered_map<uchar, ulint> last_c_pos = std::unordered_map<uchar, ulint>();

        ulint b = 0;
        ulint i = 0;
        while (i < r) 
        {
            LF_table::LF_row curr = LF_rows.get(i);

            block_chars.push_back(curr.character);
            block_intervals.push_back(curr.interval);
            block_lens.push_back(curr.length);
            block_offsets.push_back(curr.offset);

            block_base_idx.push_back(i % block_size == 0);
            for (size_t j = 1; j < curr.length; j++)
            {
                block_base_idx.push_back(false);
            }
            absolute_idx += curr.length;

            if (!last_c_pos.count(curr.character)) {
                last_c_pos[curr.character] = block_chars.size() - 1;

                // For all blocks prior without set values to find the next c's mapping, loop back and set
                if (b > 0)
                {
                    auto[k, d] = LF_rows.LF(i, curr.offset);
                    interval_pos next_c = interval_pos(k, d);

                    ulint b_curr = b;
                    while (b_curr > 0 && !blocks[b_curr-1].has_next_LF(curr.character))
                    {
                        blocks[b_curr-1].set_next_LF(curr.character, next_c);
                        --b_curr;
                    }
                }
            }
            last_c_pos[curr.character] = block_chars.size() - 1;

            ++i;

            // End of block of intervals, update block table
            if (i % block_size == 0 || i >= r)
            {        
                blocks[b] = block(block_chars, block_intervals, block_lens, block_offsets, curr_idx, prior_LF);

                for(auto const& [c, pos] : last_c_pos)
                {
                    // Since pos is wrt. current block, add the runs seen prior before computing LF
                    ulint run = b*block_size + pos;
                    // Perform LF step from the last seen character in this run (which is at offset equal to last character, one minus length)
                    auto[k, d] = LF_rows.LF(run, block_lens[pos] - 1);

                    prior_LF[c] = interval_pos(k, d);
                }

                block_chars = std::vector<uchar>();
                block_intervals = std::vector<ulint>();
                block_lens = std::vector<ulint>();
                block_offsets = std::vector<ulint>();

                curr_idx = absolute_idx;

                last_c_pos = std::unordered_map<uchar, ulint>();

                ++b;
            }
        }

        block_idx = idx_vec(block_base_idx);
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

    uchar get_char(ulint run) 
    {
        return (uchar) get_block(run).get_char(row(run));
    }

    uchar get_char(interval_pos pos) 
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

    interval_pos LF_prior(interval_pos pos, uchar c)
    {
        return reduced_pos(get_block(pos).LF_prior(row(pos), pos.offset, c));
    }

    interval_pos LF_next(interval_pos pos, uchar c)
    {
        return reduced_pos(get_block(pos).LF_next(row(pos), pos.offset, c));
    }

    interval_pos reduced_pos(interval_pos pos)
    {
        if (!pos.is_set())
        {
            return pos;
        }

        interval_pos curr = pos;
        while (curr.offset >= get_length(curr))
        {
            curr = get_block(curr).reduce(curr, row(curr));
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

    // For a general interval position, return the idx wrt. the BWT
    ulint interval_to_idx(interval_pos pos)
    {
        // ensure pos is reduced to obtain correct value
        pos = reduced_pos(pos);
        return get_block(pos).get_idx(row(pos), pos.offset);
    }

    // For a general index on the BWT, return the corresponding interval position
    interval_pos idx_to_interval(ulint idx)
    {
        assert(idx < n);
        // Get first block with idx equal to or greater than its base idx
        ulint b = block_idx.predecessor(idx);
        // Get the interval from this block, passing how many which run its first position corresponds to
        return get_block(b).get_interval(idx, b*block_size);
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

        written_bytes += block_idx.serialize(out, v, "idx_samples");

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

        block_idx.load(in);
    }
};

#endif /* end of include guard: _BLOCK_TABLE_HH */
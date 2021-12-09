/* r-index-f - Computes the simple r-index-f block compressed table
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
   \file r_index_f.hpp
   \brief r_index_f.hpp Computes the r-Index-f block table from RLBWT
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 11/19/2021
*/

#ifndef _R_INDEX_F_HH
#define _R_INDEX_F_HH

#include <common.hpp>
#include <LF_table.hpp>
#include <interval_pos.hpp>
#include <interval_block.hpp>
#include <block_table.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/dac_vector.hpp>

using namespace sdsl;
using namespace std;

template  < ulint block_size = 65536,
            class wt_t = wt_huff<bit_vector>,
            class bit_vec = bit_vector,
            class dac_vec = dac_vector<> >
class r_index_f
{
public:
    typedef std::pair<i_position, i_position> range_t;

    r_index_f() {}

    r_index_f(std::string filename)
    {
        verbose("Building the R-Index-F using Block Table Compression");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        std::string bwt_heads_fname = bwt_fname + ".heads";
        std::ifstream ifs_heads(bwt_heads_fname);
        std::string bwt_len_fname = bwt_fname + ".len";
        std::ifstream ifs_len(bwt_len_fname);

        ifs_heads.seekg(0);
        ifs_len.seekg(0);
        LF_table temp(ifs_heads, ifs_len);
        B_table = block_table(temp);
        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Block-Table construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        mem_stats();
        bwt_stats();
    }

    /* 
     * LF step from the posision given
     * \param i_position of Interval-BWT
     * \return i_position of preceding character at pos
     */
    i_position LF(interval_pos pos)
    {
        i_block* next = &B_table[next_b];
        ulint next_len;
	    while (offset >= (next_len = next->lengths[next_k])) 
        {
            offset -= next_len;
            ++next_k;
            ++q;

            if (next_k >= block_size)
            {
                next = &B_table[++next_b];
                next_k = 0;
            }
        }

	    return i_position{q, offset};
    }

    /*
    * \param r inclusive range of a string w
    * \param c character
    * \return inclusive range of cw
    */
    range_t LF(range_t range, char c)
    {
        i_position start;
        assert(range.first.run < r);

        ulint b = range.first.run/block_size;
        ulint k = range.first.run%block_size;
        i_block* curr = &B_table[b];

        ulint offset = range.first.offset;

        auto [c_rank_f, bwt_c_f] = curr->heads.inverse_select(k);
        if (c != bwt_c_f)
        {
            c_rank_f = curr->heads.rank(k, c);
            if (curr->heads.rank(curr->heads.size(),c) < c_rank_f + 1)
            {
                //if (!curr->next_is_valid[c])
                //{
                //    return range_t(i_position{1,0}, i_position{0,0});
                //}
                //else
                //{
                    start = curr->get_next_LF(c);
                //}
            }
            else {
                k = curr->heads.select(c_rank_f + 1, c);
                offset = 0;

                start = LF_from_rank(curr, k, offset, c_rank_f, c);
            }
        }
        else
        {
            start = LF_from_rank(curr, k, offset, c_rank_f, c);
        }

        i_position end;
        assert(range.second.run < r);

        b = range.second.run/block_size;
        k = range.second.run%block_size;
        curr = &B_table[b];

        offset = range.second.offset;
      
        auto [c_rank_s, bwt_c_s] = curr->heads.inverse_select(k);
        if (c != bwt_c_s)
        {
            c_rank_s = curr->heads.rank(k, c);
            if (c_rank_s == 0)
            {
                //if (!curr->prior_is_valid[c])
                //{
                //    return range_t(i_position{1,0}, i_position{0,0});
                //}
                //else
                //{
                    end = curr->get_prior_LF(c);
                //}
            }
            else
            {
                k = curr->heads.select(c_rank_s, c);
                offset = curr->lengths[k] - 1;

                end = LF_from_rank(curr, k, offset, c_rank_s - 1, c);
            }
        }
        else
        {
            end = LF_from_rank(curr, k, offset, c_rank_s, c);
        }

        return range_t(start, end);
    }

    ulint runs()
    {
        return B_table.runs();
    }

    ulint size()
    {
        return B_table.size();
    }

    range_t full_range()
    {
        return range_t(B_table.begin(), B_table.end());
    }

    void mem_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("              Block table: ", serialize(ns));
    }

    void bwt_stats()
    {
        ulint n = size();
        ulint r = runs();
        verbose("Number of BWT equal-letter runs: r = ", r);
        verbose("Length of complete BWT: n = ", n);
        verbose("Rate n/r = ", double(n) / r);
        verbose("log2(r) = ", log2(double(r)));
        verbose("log2(n/r) = ", log2(double(n) / r));
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += B_table.serialize(out, v, "B_table");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    std::string get_file_extension() const
    {
        return ".rif";
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {
        B_table.load(in);
    }

private:
    block_table<block_size, wt_t, bit_vec, dac_vec> B_table;

    i_position LF_from_rank(i_block* curr, ulint k, ulint offset, ulint c_rank, char c)
    {
        ulint q = curr->get_interval(c, c_rank);

        offset += curr->offsets[k];

        ulint next_b = q/block_size;
        ulint next_k = q%block_size;
        i_block* next = &B_table[next_b];
        ulint next_len;
	    while (offset >= (next_len = next->lengths[next_k])) 
        {
            offset -= next_len;
            ++next_k;
            ++q;

            if (next_k >= block_size)
            {
                next = &B_table[++next_b];
                next_k = 0;
            }
        }

        return i_position{q, offset};
    }
};

#endif /* end of include guard: _R_INDEX_F_HH */
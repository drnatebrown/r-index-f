/* Lf_table - Uncompressed version of OptBWTR (LF table)
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
   \file alt_table.hpp
   \brief alt_table.hpp Columnar of OptBWTR (LF table) using run lengths and relative run-head LF destinations
   \author Nathaniel Brown
   \date 24/04/2024
*/

#ifndef _ALT_TABLE_HH
#define _ALT_TABLE_HH

#include <sdsl/dac_vector.hpp>
#include <common.hpp>

#include <ds/heads_wt_w.hpp>
#include <ds/intervals_rank_w.hpp>
#include <ds/interval_pos.hpp>

#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace std;

class alt_table
{
public:
    alt_table() {}

    alt_table(std::ifstream &bwt) {}

    alt_table(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        
        std::vector<uchar> chars = std::vector<uchar>();
        std::vector<ulint> lens = std::vector<ulint>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
        char c;
        ulint i = 0;
        n = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c <= TERMINATOR) c = TERMINATOR;

            chars.push_back(c);
            lens.push_back(length);
            L_block_indices[c].push_back(i++);

            n+=length;
        }
        r = i;

        std::vector<ulint> lf_dest = std::vector<ulint>(r);
        std::vector<ulint> offs = std::vector<ulint>(r);

        ulint curr_L_num = 0;
        ulint L_seen = 0;
        ulint F_seen = 0;
        for(size_t i = 0; i < L_block_indices.size(); ++i) 
        {
            for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
            {
                ulint pos = L_block_indices[i][j];

                lf_dest[pos] = curr_L_num;
                offs[pos] = F_seen - L_seen;

                F_seen += lens[pos];
            
                while (curr_L_num < r && F_seen >= L_seen + lens[curr_L_num]) 
                {
                    L_seen += lens[curr_L_num];
                    ++curr_L_num;
                }
            }
        }

        run_heads = run_heads_t(chars);
        dest_off = offsets_t(offs);
        dest_pred = dest_pred_t(chars, lf_dest);
        run_len = lengths_t(lens);

        mem_stats();
    }

    ulint size()
    {
        return n;
    }

    ulint runs()
    {
        return r;
    }

    uchar get_char(ulint i)
    {
        return run_heads[i];
    }

    uchar get_char(interval_pos pos) 
    {
        return get_char(pos.run);
    }

    interval_pos begin()
    {
        return interval_pos(0, 0);
    }

    interval_pos end()
    {
        return interval_pos(r-1, run_len[r-1] - 1);
    }

    interval_pos LF(interval_pos pos)
    {
        auto [c_rank, c] = run_heads.inverse_select(pos.run);
        return LF(pos.run, pos.offset, c_rank, c);
    }

    interval_pos LF_prior(interval_pos pos, uchar c)
    {
        ulint c_rank = run_heads.rank(pos.run + 1, c);
        if (c_rank == 0) {
            return interval_pos();
        }
        else {
            c_rank -= 1;
        }

        ulint prior_run = run_heads.select(c_rank + 1, c);
        ulint prior_off = (pos.run != prior_run) ? run_len[prior_run] - 1 : pos.offset;
        
        return LF(prior_run, prior_off, c, c_rank);
    }

    interval_pos LF_next(interval_pos pos, uchar c)
    {
        ulint c_rank = run_heads.rank(pos.run, c);
        if (c_rank + 1 > run_heads.rank(r, c)) return interval_pos();

        ulint next_run = run_heads.select(c_rank + 1, c);
        ulint next_off = (pos.run != next_run) ? 0 : pos.offset;

        return LF(next_run, next_off, c, c_rank);
    }

    interval_pos reduced_pos(interval_pos pos)
    {
        if (!pos.is_set())
        {
            return pos;
        }

        ulint run = pos.run;
        ulint off = pos.offset;

        while (off >= run_len[run]) 
        {
            off -= run_len[run++];
        }

        return interval_pos(run, off);
    }

    // For a general interval position, return the idx wrt. the BWT
    ulint interval_to_idx(interval_pos pos)
    {
        ulint curr = 0;
        ulint idx = 0;
        while(curr < pos.run) {
            idx += run_len[curr++];
        }
        return idx + pos.offset;
    }

    // For a general index on the BWT, return the corresponding interval position
    interval_pos idx_to_interval(ulint idx)
    {
        ulint curr = 0;
        ulint curr_idx = 0;
        while(idx - curr_idx > run_len[curr]) {
            curr_idx += run_len[curr++];
        }
        return interval_pos(curr, idx - curr_idx);
    }

    std::string get_file_extension() const
    {
        return ".alt_table";
    }

    void mem_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("  Columnar LF table: ", serialize(ns));
        verbose("                                 ");
        verbose("              Run_heads: ", run_heads.serialize(ns));
        verbose("                Run_len: ", run_len.serialize(ns));
        verbose("              Dest_pred: ", dest_pred.serialize(ns));
        verbose("               Dest_off: ", dest_off.serialize(ns));
    }

    void bwt_stats()
    {
        verbose("Number of BWT equal-letter runs: r = ", r);
        verbose("Length of complete BWT: n = ", n);
        verbose("Rate n/r = ", double(n) / r);
        verbose("log2(r) = ", log2(double(r)));
        verbose("log2(n/r) = ", log2(double(n) / r));
    }

    /* serialize to the ostream
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

        written_bytes += run_heads.serialize(out, v, "Run_heads");
        written_bytes += dest_off.serialize(out, v, "Offsets");
        written_bytes += dest_pred.serialize(out, v, "Dest_pred");
        written_bytes += run_len.serialize(out, v, "Lengths");

        return written_bytes;
    }

    /* load from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));

        run_heads.load(in);
        run_len.load(in);
        dest_pred.load(in);
        dest_off.load(in);
    }

private:
    ulint n; // Length of BWT
    ulint r; // Runs of BWT

    typedef heads_wt_w<> run_heads_t; // Huffman-Shaped WT
    typedef intervals_rank_w<base_bv<>, symbol_map> dest_pred_t; // Plain Bitvector per Character
    typedef sdsl::dac_vector_dp<> lengths_t; // Sparse Bitvector per Character
    typedef sdsl::dac_vector_dp<> offsets_t; // Sparse Bitvector

    run_heads_t run_heads;
    lengths_t run_len;
    offsets_t dest_off;
    dest_pred_t dest_pred;

    interval_pos LF(ulint k, ulint d, uchar c, ulint c_rank) {
        ulint next_interval = dest_pred.get(c, c_rank);
        ulint next_off = dest_off[k] + d;
        return reduced_pos(interval_pos(next_interval, next_off));
    }
};

#endif /* end of include guard: _LF_TABLE_HH */
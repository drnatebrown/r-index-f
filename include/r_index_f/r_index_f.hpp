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
#include <block_table.hpp>

#include <ds/interval_block.hpp>
#include <ds/interval_pos.hpp>
#include <ds/idx_bit_vector.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/dac_vector.hpp>

using namespace sdsl;
using namespace std;

template <class table = block_table<>>
class r_index_f
{
public:
    typedef std::pair<interval_pos, interval_pos> range_t;

    r_index_f() {}

    r_index_f(std::string filename, bool rle = true)
    {
        verbose("Building the R-Index-F using Block Table Compression");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        if (rle)
        {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            LF_table temp(ifs_heads, ifs_len);
            B_table = table(temp);
        }
        else
        {
            std::ifstream ifs_bwt(bwt_fname);

            ifs_bwt.seekg(0);
            LF_table temp(ifs_bwt);
            B_table = table(temp);
        }

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Block-Table construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        mem_stats();
        bwt_stats();
    }

    r_index_f(LF_table t) {
        verbose("Building the R-Index-F using Block Table Compression from LF Table Construction");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
        B_table = table(t);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Block-Table construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        mem_stats();
        bwt_stats();
    }

    ulint runs()
    {
        return B_table.runs();
    }

    ulint size()
    {
        return B_table.size();
    }

    void count(const std::string &pattern, vector<ulint> tnls, std::ofstream &outfile) {
        range_t range = full_range();
        ulint m = pattern.size();

        // TUNNEL STATS
        ulint tnl_hit = 0;
        ulint total_hit = 0;
        ulint pattern_len = 0;
        ulint tunnel_steps = 0;
        ulint tunnel = false;
        ulint width = 0;
        
        struct tunnel_stat {
            ulint run;
            ulint height;
            ulint width;
            ulint range_size;
            ulint steps;
            ulint pattern_len;
        };

        tunnel_stat curr;

        for (ulint i=0; i < m && range.second >= range.first; ++i){
            range = LF(range, pattern[m - i - 1]);
            
            pattern_len++;
            if (tunnel && (width == 0 || get_char(range.first) != pattern[m - i - 1])) {
                tunnel = false;
                width = 0;

                outfile << curr.run << "\t" << curr.height << "\t" << curr.width << "\t" << curr.range_size << "\t" << curr.steps << "\t" << curr.pattern_len + curr.steps << std::endl;
                range = full_range();
                pattern_len = 1;
            }
            if (range.first.run == range.second.run) {
                total_hit++;
                if (get_char(range.first) == pattern[m - i - 1] && tnls[range.first.run] > 1) {
                    tnl_hit++;
                    tunnel = true;
                    width = tnls[range.first.run];
                    curr = {range.first.run, B_table.get_length(range.first.run), width, (range.second.offset - range.first.offset) + 1, 0, pattern_len};
                }
                else {
                    // reset query if falls within a non-tunnel run
                    outfile << range.first.run << "\t" << B_table.get_length(range.first.run) << "\t" << 1 << "\t" << (range.second.offset - range.first.offset) << "\t" << 0 << "\t" << pattern_len << std::endl;
                    range = full_range();
                    pattern_len = 1;
                }
            }

            curr.steps += tunnel;
            if (width) --width;
        }
        if (tunnel) {
            outfile << curr.run << "\t" << curr.height << "\t" << curr.width << "\t" << curr.range_size << "\t" << curr.steps << "\t" << curr.pattern_len + curr.steps << std::endl;
        }
        else if (range.first.run == range.second.run) {
            total_hit++;
            outfile << range.first.run << "\t" << B_table.get_length(range.first.run) << "\t" << 1 << "\t" << (range.second.offset - range.first.offset) << "\t" << 0 << "\t" << curr.pattern_len << std::endl;
        }

        cout << tnl_hit << "/" << total_hit << " = " << std::setprecision(4) << ((1.0*tnl_hit)/total_hit)*100 << "%" << std::endl;
        // return interval_to_idx(range.second) - interval_to_idx(range.first) + 1;  
    }

    // void invert() {
    //     std::string outfile = args.filename + ".inverted";
    //     std::ofstream out(outfile);

    //     interval_pos i = {0,0};
    //     char c;
    //     while((c = bwt.get_char(i)) > TERMINATOR)
    //         out << c;
    //         i = bwt.LF(i), ++steps;
    //     out.close();
    // }

    interval_pos LF(interval_pos pos)
    {
        return B_table.LF(pos);
    }

    range_t LF(range_t range, uchar c)
    {
        return range_t(B_table.LF_next(range.first, c), B_table.LF_prior(range.second, c));
    }

    range_t full_range()
    {
        return range_t(B_table.begin(), B_table.end());
    }

    ulint interval_to_idx(interval_pos pos)
    {
        return B_table.interval_to_idx(pos);
    }

    uchar get_char(interval_pos pos)
    {
        return B_table.get_char(pos);
    }

    // Return underlying table (not recommended, add methods to access its capabilities)
    table get_table()
    {
        return B_table;
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
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

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
    table B_table;
};

#endif /* end of include guard: _R_INDEX_F_HH */
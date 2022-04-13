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
   \file LF_table.hpp
   \brief LF_table.hpp Uncompressed version of OptBWTR (LF table)
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 19/11/2021
*/

#ifndef _LF_TABLE_HH
#define _LF_TABLE_HH

#include <common.hpp>

#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace std;

class LF_table
{
public:
    // Row of the LF table
    typedef struct LF_row
    {
        char character;
        ulint length;
        ulint interval;
        ulint offset;

        size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
        {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_t written_bytes = 0;

            out.write((char *)&character, sizeof(character));
            written_bytes += sizeof(character);

            out.write((char *)&interval, sizeof(interval));
            written_bytes += sizeof(interval);

            out.write((char *)&length, sizeof(length));
            written_bytes += sizeof(length);

            out.write((char *)&offset, sizeof(offset));
            written_bytes += sizeof(offset);

            return written_bytes;
        }

        void load(std::istream &in)
        {
            in.read((char *)&character, sizeof(character));
            in.read((char *)&interval, sizeof(interval));
            in.read((char *)&length, sizeof(length));
            in.read((char *)&offset, sizeof(offset));
        }
    };

    LF_table() {}

    // TODO: Add builder for BWT (not heads/lengths)
    LF_table(std::ifstream &heads, std::ifstream &lengths, ulint max_run = 0)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        
        LF_runs = vector<LF_row>();
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
        char c;
        ulint i = 0;
        r = 0;
        n = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c <= TERMINATOR) c = TERMINATOR;

            if (max_run > 0 && length > max_run) {
                ulint max_splits = length/max_run;
                for (size_t split = 0; split < max_splits; ++split)
                {
                    LF_runs.push_back({c, max_run, 0, 0});
                    L_block_indices[c].push_back(i++);
                }

                if (length % max_run != 0)
                {
                    LF_runs.push_back({c, length % max_run, 0, 0});
                    L_block_indices[c].push_back(i++);
                }
            }
            else {
                LF_runs.push_back({c, length, 0, 0});
                L_block_indices[c].push_back(i++);
            }    
            n+=length;
        }
        r = LF_runs.size();

        ulint curr_L_num = 0;
        ulint L_seen = 0;
        ulint F_seen = 0;
        for(size_t i = 0; i < L_block_indices.size(); ++i) 
        {
            for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
            {
                ulint pos = L_block_indices[i][j];

                LF_runs[pos].interval = curr_L_num;
                LF_runs[pos].offset = F_seen - L_seen;

                F_seen += LF_runs[pos].length;
            
                while (curr_L_num < r && F_seen >= L_seen + LF_runs[curr_L_num].length) 
                {
                    L_seen += LF_runs[curr_L_num].length;
                    ++curr_L_num;
                }
            }
        }

        mem_stats();
    }

    const LF_row get(size_t i)
    {
        assert(i < LF_runs.size());
        return LF_runs[i];
    }

    ulint size()
    {
        return n;
    }

    ulint runs()
    {
        return r;
    }

    void invert(std::string outfile) 
    {
        std::ofstream out(outfile);

        ulint interval = 0;
        ulint offset = 0;

        char c;
        while((c = get_char(interval)) > TERMINATOR) 
        {
            out << c;
            std::pair<ulint, ulint> pos = LF(interval, offset);
            interval = pos.first;
            offset = pos.second;
        }
    }

    /*
     * \param Run position (RLE intervals)
     * \param Current character offset in block
     * \return block position and offset of preceding character
     */
    std::pair<ulint, ulint> LF(ulint run, ulint offset)
    {
        ulint next_interval = LF_runs[run].interval;
        ulint next_offset = LF_runs[run].offset + offset;

        while (next_offset >= LF_runs[next_interval].length) 
        {
            next_offset -= LF_runs[next_interval++].length;
        }

        return std::make_pair(next_interval, next_offset);
    }

    uchar get_char(ulint i)
    {
        return get(i).character;
    }

    std::string get_file_extension() const
    {
        return ".LF_table";
    }

    void mem_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("              LF table: ", serialize(ns));
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

        size_t size = LF_runs.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        for(size_t i = 0; i < size; ++i)
        {
            written_bytes += LF_runs[i].serialize(out, v, "LF_run_" + std::to_string(i));
        }

        return written_bytes;
    }

    /* load from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        size_t size;

        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));

        in.read((char *)&size, sizeof(size));
        LF_runs = std::vector<LF_row>(size);
        for(size_t i = 0; i < size; ++i)
        {
            LF_runs[i].load(in);
        }
    }

private:
    ulint n; // Length of BWT
    ulint r; // Runs of BWT

    vector<LF_row> LF_runs;
};

#endif /* end of include guard: _LF_TABLE_HH */
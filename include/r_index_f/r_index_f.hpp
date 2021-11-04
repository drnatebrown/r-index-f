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
   \date 09/07/2020
*/

#ifndef _R_INDEX_F_HH
#define _R_INDEX_F_HH

#include <common.hpp>

#include <utility>
#include <iostream> 
#include <algorithm>
#include <random>
#include <vector>
#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/dac_vector.hpp>

using namespace sdsl;
using namespace std;

static const uint8_t TERMINATOR = 1;
typedef unsigned long int ulint;

template    <ulint block_size = 1048576,
            class wt_t = wt_huff<bit_vector>,
            class bit_vec = bit_vector,
            class dac_vec = dac_vector<>>
class r_index_f
{
public:
    typedef size_t size_type;
    typedef bit_vector::select_1_type bv_select_1;

    struct i_block
    {
        wt_t heads;

        ulint A_map;
        bv_select_1 A_diff;
        bit_vec A_bv;

        ulint C_map;
        bv_select_1 C_diff;
        bit_vec C_bv;

        ulint G_map;
        bv_select_1 G_diff;
        bit_vec G_bv;

        ulint T_map;
        bv_select_1 T_diff;
        bit_vec T_bv;

        std::unordered_map<char, ulint> else_map;
        std::unordered_map<char, bv_select_1> else_diff;
        std::unordered_map<char, bit_vec> else_bv;

        dac_vec lengths;
        dac_vec offsets;

        const ulint get_interval(const char c, const ulint d)
        {
            switch(c)
            {
                case 'A':
                    return A_map + A_diff(d+1) - d;

                case 'C':
                    return C_map + C_diff(d+1) - d;

                case 'G':
                    return G_map + G_diff(d+1) - d;

                case 'T':
                    return T_map + T_diff(d+1) - d;

                default:
                    return else_map[c] + else_diff[c](d+1) - d;
            }
        }

        /* serialize the structure to the ostream
        * \param out     the ostream
        */
        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
        {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            written_bytes += heads.serialize(out,v,"Heads");
            
            out.write((char *)&A_map, sizeof(A_map));
            written_bytes += sizeof(A_map);
            written_bytes += A_bv.serialize(out,v,"A_bv");
            
            out.write((char *)&C_map, sizeof(C_map));
            written_bytes += sizeof(C_map);
            written_bytes += C_bv.serialize(out,v,"C_bv");

            out.write((char *)&G_map, sizeof(G_map));
            written_bytes += sizeof(G_map);
            written_bytes += G_bv.serialize(out,v,"G_bv");

            out.write((char *)&T_map, sizeof(T_map));
            written_bytes += sizeof(T_map);
            written_bytes += T_bv.serialize(out,v,"T_bv");

            size_t size = else_map.size();
            out.write((char *)&size, sizeof(size_t));
            written_bytes += sizeof(size);
            for(auto const& [key, val] : else_map)
            {
                out.write((char *)&key, sizeof(key));
                written_bytes += sizeof(key);                
                out.write((char *)&val, sizeof(val));
                written_bytes += sizeof(val);                
            }

            size = else_bv.size();
            out.write((char *)&size, sizeof(size_t));
            written_bytes += sizeof(size);
            for(auto const& [key, val] : else_bv)
            {
                out.write((char *)&key, sizeof(key));
                written_bytes += sizeof(key);                
                written_bytes += val.serialize(out,v,"else_bv_" + std::to_string(key));                
            }

            written_bytes += lengths.serialize(out,v,"Lengths");
            written_bytes += offsets.serialize(out,v,"Offsets");

            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        /* load the structure from the istream
        * \param in the istream
        */
        void load(std::istream &in)
        {
            heads.load(in);
            
            in.read((char *)&A_map, sizeof(A_map));
            A_bv.load(in);
            A_diff = bv_select_1(&A_bv);

            in.read((char *)&C_map, sizeof(C_map));
            C_bv.load(in);
            C_diff = bv_select_1(&C_bv);

            in.read((char *)&G_map, sizeof(G_map));
            G_bv.load(in);
            G_diff = bv_select_1(&G_bv);

            in.read((char *)&T_map, sizeof(T_map));
            T_bv.load(in);
            T_diff = bv_select_1(&T_bv);

            size_t size;
            in.read((char *)&size, sizeof(size));
            for(size_t i = 0; i < size; ++i)
            {
                char key;
                ulint val;
                in.read((char *)&key, sizeof(key));
                in.read((char *)&val, sizeof(val));
                else_map[key] = val;
            }

            in.read((char *)&size, sizeof(size));
            for(size_t i = 0; i < size; ++i)
            {
                char key;
                bit_vec val;
                in.read((char *)&key, sizeof(key));
                val.load(in);
                else_bv[key] = val;
                else_diff[key] = bv_select_1(&else_bv[key]);
            }

            lengths.load(in);
            offsets.load(in);
        }
    };

    ulint n;
    ulint r;
    vector<i_block> B_table; 

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
        build_B_table(ifs_heads, ifs_len);
        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Block-Table construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        mem_stats();
        bwt_stats();
    }

    vector<i_block> build_B_table(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(256);
        vector<char> chars = vector<char>(); 
        vector<ulint> lens = vector<ulint>();
        
        char c;
        ulint i = 0;
        n = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
            {
                chars.push_back(c);
                lens.push_back(length);
                L_block_indices[c].push_back(i);
            }
            else
            {
                chars.push_back(TERMINATOR);
                lens.push_back(length);
                L_block_indices[TERMINATOR].push_back(i);
            }
            ++i;
            n+=length;
        }
        
        r = chars.size();

        vector<ulint> intervals = vector<ulint>(r);
        vector<ulint> offsets = vector<ulint>(r);

        ulint curr_L_num = 0;
        ulint L_seen = 0;
        ulint F_seen = 0;
        for(size_t i = 0; i < L_block_indices.size(); ++i) 
        {
            for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
            {
                ulint pos = L_block_indices[i][j];

                intervals[pos] = curr_L_num;
                offsets[pos] = F_seen - L_seen;

                F_seen += lens[pos];
            
                while (curr_L_num < r && F_seen >= L_seen + lens[curr_L_num]) 
                {
                    L_seen += lens[curr_L_num];
                    ++curr_L_num;
                }
            }
        }

        ulint B_len = (r/block_size) + ((r % block_size) != 0);
        B_table = vector<i_block>(B_len);

        vector<char> block_chars = vector<char>(block_size);
        vector<ulint> block_lens = vector<ulint>(block_size);
        vector<ulint> block_offsets = vector<ulint>(block_size);
        std::unordered_map<char, ulint> block_c_map = std::unordered_map<char, ulint>();
        std::unordered_map<char, ulint> last_c_map = std::unordered_map<char, ulint>();
        std::unordered_map<char, vector<bool>> bit_diff = std::unordered_map<char, vector<bool>>();

        ulint b = 0;
        ulint b_i = 0;
        i = 0;
        while (i < r) 
        {
            char c = chars[i];
            ulint l = lens[i];
            ulint k = intervals[i];
            ulint d = offsets[i];

            block_chars[b_i] = c;
            block_lens[b_i] = l;
            block_offsets[b_i] = d;

            if (!block_c_map.count(c)) {
                block_c_map.insert(std::pair<char, ulint>(c, k));
                last_c_map.insert(std::pair<char, ulint>(c, k));
                bit_diff.insert(std::pair<char, vector<bool>>(c, vector<bool>()));
            }

            ulint diff = k - last_c_map[c];
            while (diff > 0) {
                bit_diff[c].push_back(false);
                --diff;
            }
            bit_diff[c].push_back(true);

            last_c_map[c] = k;

            ++i;
            ++b_i;

            // End of block of intervals, update block table
            if (b_i >= block_size || i >= r)
            {
                i_block& curr = B_table[b];

                construct_im(curr.heads, std::string(block_chars.begin(), (i >= r) ? block_chars.begin()+b_i : block_chars.end()).c_str(), 1);

                for (auto& kv: bit_diff) 
                {
                    switch(kv.first)
                    {
                        case 'A':
                            curr.A_map = block_c_map['A'];
                            curr.A_bv = bv(kv.second);
                            curr.A_diff = bv_select_1(&curr.A_bv);
                            break;

                        case 'C':
                            curr.C_map = block_c_map['C'];
                            curr.C_bv = bv(kv.second);
                            curr.C_diff = bv_select_1(&curr.C_bv);
                            break;

                        case 'G':
                            curr.G_map = block_c_map['G'];
                            curr.G_bv = bv(kv.second);
                            curr.G_diff = bv_select_1(&curr.G_bv);
                            break;

                        case 'T':
                            curr.T_map = block_c_map['T'];
                            curr.T_bv = bv(kv.second);
                            curr.T_diff = bv_select_1(&curr.T_bv);
                            break;

                        default:
                            curr.else_map.insert(std::pair<char, ulint>(kv.first, block_c_map[kv.first]));
                            curr.else_bv.insert(std::pair<char, bit_vec>(kv.first, bv(kv.second)));
                            curr.else_diff.insert(std::pair<char, bv_select_1>(kv.first, bv_select_1(&curr.else_bv[kv.first])));
                            break;
                    }
                }

                curr.lengths = dac_vec(block_lens);
                curr.offsets = dac_vec(block_offsets);
                
                block_chars = vector<char>(block_size);
                block_lens = vector<ulint>(block_size);
                block_offsets = vector<ulint>(block_size);
                block_c_map = std::unordered_map<char, ulint>();
                last_c_map = std::unordered_map<char, ulint>();
                bit_diff = std::unordered_map<char, vector<bool>>();

                ++b;
                b_i = 0;
            }
        }

        return B_table;
    }

    bit_vec bv(vector<bool> &b){

		if(b.size()==0) return bit_vector();

		bit_vector bv(b.size());

		for(uint64_t i=0;i<b.size();++i)
			bv[i] = b[i];

		return bit_vector(bv);
	}

    /*
     * \param Run position (RLBWT)
     * \param Current offset within interval
     * \param Character to step from
     * \return run position and offset of preceding character
     */
    std::pair<ulint, ulint> LF(ulint run, ulint offset, char c)
    {
        ulint b = run/block_size;
        ulint k = run%block_size;
        i_block* curr = &B_table[b];

        auto [c_rank, bwt_c] = curr->heads.inverse_select(k);
        if (c != bwt_c)
        {
            c_rank = curr->heads.rank(k, c);
            while (c_rank == 0)
            {
                if (b == 0)
                {
                    error("No preceding character for position, cannot LF step");
                    throw std::logic_error("Character" + util::to_string(c) + "does not occur anywhere preceding run" + util::to_string(run) +".");
                }
                curr =  &B_table[--b];
                c_rank = curr->heads.rank(block_size, c);
            }

            k = curr->heads.select(k+1, c);
            offset = curr->lengths[k] - 1;
        }

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

	    return std::make_pair(q, offset);
    }

    char get_char(ulint run) 
    {
        return (char) B_table[run/block_size].heads[run%block_size];
    }

    void mem_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("              Block table: ", serialize(ns));
    }

    void bwt_stats()
    {
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

        size_t size = B_table.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        for(size_t i = 0; i < size; ++i)
            written_bytes += B_table[i].serialize(out,v,"B_table_" + std::to_string(i));

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
        size_t size;
        in.read((char *)&size, sizeof(size));
        B_table = std::vector<i_block>(size);
        for(size_t i = 0; i < size; ++i)
        {
            B_table[i].load(in);
        }
    }
};

#endif /* end of include guard: _R_INDEX_F_HH */
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

template  < ulint block_size = 1048576,
            class wt_t = wt_huff<bit_vector>,
            class bit_vec = bit_vector,
            class dac_vec = dac_vector<> >
class r_index_f
{
public:
    typedef size_t size_type;

    struct i_position
    {
        ulint run;
        ulint offset;

        i_position& operator++() 
        {
            ++offset;
            if (offset >= B_table[run/block_size].lengths[run%block_size])
            {
                ++run;
                offset = 0;
            }
            return *this;
        }

        i_position operator++(int) 
        {
            i_position old = *this;
            operator++();
            return old;
        }

        inline bool operator< (const i_position& pos){ return (run == pos.run) ? (offset < pos.offset) : (run < pos.run); }
        inline bool operator> (const i_position& pos){ return pos < *this; }
        inline bool operator<=(const i_position& pos){ return !(*this > pos); }
        inline bool operator>=(const i_position& pos){ return !(*this < pos); }
        inline bool operator==(const i_position& pos){ return run == pos.run && offset == pos.offset; }
        inline bool operator!=(const i_position& pos){ return !(*this == pos); }

        /* serialize the structure to the ostream
        * \param out     the ostream
        */
        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
        {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            out.write((char *)&run, sizeof(run));
            written_bytes += sizeof(run);

            out.write((char *)&offset, sizeof(offset));
            written_bytes += sizeof(offset);

            return written_bytes;
        }

        /* load the structure from the istream
        * \param in the istream
        */
        void load(std::istream &in)
        {
            
        }
    };

    typedef bit_vector::select_1_type bv_select_1;
    typedef std::pair<i_position, i_position> range_t;

    struct i_block
    {
        wt_t heads;

        ulint A_map;
        bv_select_1 A_diff;
        bit_vec A_bv;
        i_position next_A_LF;
        i_position prior_A_LF;

        ulint C_map;
        bv_select_1 C_diff;
        bit_vec C_bv;
        i_position next_C_LF;
        i_position prior_C_LF;

        ulint G_map;
        bv_select_1 G_diff;
        bit_vec G_bv;
        i_position next_G_LF;
        i_position prior_G_LF;

        ulint T_map;
        bv_select_1 T_diff;
        bit_vec T_bv;
        i_position next_T_LF;
        i_position prior_T_LF;

        std::unordered_map<char, ulint> else_map;
        std::unordered_map<char, bv_select_1> else_diff;
        std::unordered_map<char, bit_vec> else_bv;
        std::unordered_map<char, i_position> else_next_LF;
        std::unordered_map<char, i_position> else_prior_LF;

        dac_vec lengths;
        dac_vec offsets;

        ulint idx;

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

        i_position get_prior_LF(const char c)
        {
            switch(c)
            {
                case 'A':
                    return prior_A_LF;

                case 'C':
                    return prior_C_LF;

                case 'G':
                    return prior_G_LF;

                case 'T':
                    return prior_T_LF;

                default:
                    return else_prior_LF[c];
            }
        }

        i_position get_next_LF(const char c)
        {
            switch(c)
            {
                case 'A':
                    return next_A_LF;

                case 'C':
                    return next_C_LF;

                case 'G':
                    return next_G_LF;

                case 'T':
                    return next_T_LF;

                default:
                    return else_next_LF[c];
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
            written_bytes += next_A_LF.serialize(out,v,"Next_A_LF");
            written_bytes += prior_A_LF.serialize(out,v,"Prior_A_LF");
            
            out.write((char *)&C_map, sizeof(C_map));
            written_bytes += sizeof(C_map);
            written_bytes += C_bv.serialize(out,v,"C_bv");
            written_bytes += next_C_LF.serialize(out,v,"Next_C_LF");
            written_bytes += prior_C_LF.serialize(out,v,"Prior_C_LF");

            out.write((char *)&G_map, sizeof(G_map));
            written_bytes += sizeof(G_map);
            written_bytes += G_bv.serialize(out,v,"G_bv");
            written_bytes += next_G_LF.serialize(out,v,"Next_G_LF");
            written_bytes += prior_G_LF.serialize(out,v,"Prior_G_LF");

            out.write((char *)&T_map, sizeof(T_map));
            written_bytes += sizeof(T_map);
            written_bytes += T_bv.serialize(out,v,"T_bv");
            written_bytes += next_T_LF.serialize(out,v,"Next_T_LF");
            written_bytes += prior_T_LF.serialize(out,v,"Prior_T_LF");

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

            size = else_next_LF.size();
            out.write((char *)&size, sizeof(size_t));
            written_bytes += sizeof(size);
            for(auto const& [key, val] : else_next_LF)
            {
                out.write((char *)&key, sizeof(key));
                written_bytes += sizeof(key);                
                out.write((char *)&val, sizeof(val));
                written_bytes += sizeof(val);                
            }

            size = else_prior_LF.size();
            out.write((char *)&size, sizeof(size_t));
            written_bytes += sizeof(size);
            for(auto const& [key, val] : else_prior_LF)
            {
                out.write((char *)&key, sizeof(key));
                written_bytes += sizeof(key);                
                out.write((char *)&val, sizeof(val));
                written_bytes += sizeof(val);                
            }

            written_bytes += lengths.serialize(out,v,"Lengths");
            written_bytes += offsets.serialize(out,v,"Offsets");

            out.write((char *)&idx, sizeof(idx));
            written_bytes += sizeof(idx);

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
            next_A_LF.load(in);
            prior_A_LF.load(in);

            in.read((char *)&C_map, sizeof(C_map));
            C_bv.load(in);
            C_diff = bv_select_1(&C_bv);
            next_C_LF.load(in);
            prior_C_LF.load(in);

            in.read((char *)&G_map, sizeof(G_map));
            G_bv.load(in);
            G_diff = bv_select_1(&G_bv);
            next_G_LF.load(in);
            prior_G_LF.load(in);

            in.read((char *)&T_map, sizeof(T_map));
            T_bv.load(in);
            T_diff = bv_select_1(&T_bv);
            next_T_LF.load(in);
            prior_T_LF.load(in);

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

            in.read((char *)&size, sizeof(size));
            for(size_t i = 0; i < size; ++i)
            {
                char key;
                i_position val;
                in.read((char *)&key, sizeof(key));
                val.load(in);
                else_next_LF[key] = val;
            }

            in.read((char *)&size, sizeof(size));
            for(size_t i = 0; i < size; ++i)
            {
                char key;
                i_position val;
                in.read((char *)&key, sizeof(key));
                val.load(in);
                else_prior_LF[key] = val;
            }

            lengths.load(in);
            offsets.load(in);

            in.read((char *)&idx, sizeof(idx));
        }
    };

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
        std::unordered_map<char, ulint> prior_last_c_map = std::unordered_map<char, ulint>();
        std::unordered_map<char, vector<bool>> bit_diff = std::unordered_map<char, vector<bool>>();
        ulint block_idx = 0;
        ulint next_idx = 0;

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

            next_idx += l;

            if (!block_c_map.count(c)) {
                block_c_map.insert(std::pair<char, ulint>(c, k));
                last_c_map.insert(std::pair<char, ulint>(c, k));
                bit_diff.insert(std::pair<char, vector<bool>>(c, vector<bool>()));

                if (b > 0)
                {
                    ulint c_b = k;
	                ulint c_off = d;

	                while (c_off >= lens[c_b]) 
                    {
                        c_off -= lens[c_b];
                        ++c_b;
                    }

                    switch(c)
                    {
                        case 'A':
                            B_table[b-1].next_A_LF = i_position{c_b, c_off};
                            break;

                        case 'C':
                            B_table[b-1].next_C_LF = i_position{c_b, c_off};
                            break;

                        case 'G':
                            B_table[b-1].next_G_LF = i_position{c_b, c_off};
                            break;

                        case 'T':
                            B_table[b-1].next_T_LF = i_position{c_b, c_off};
                            break;

                        default:
                            B_table[b-1].else_next_LF.insert(std::pair<char, i_position>(c, i_position{c_b, c_off}));
                            break;
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
            ++b_i;

            // End of block of intervals, update block table
            if (b_i >= block_size || i >= r)
            {
                i_block& curr = B_table[b];

                construct_im(curr.heads, std::string(block_chars.begin(), (i >= r) ? block_chars.begin()+b_i : block_chars.end()).c_str(), 1);

                for (auto& kv: bit_diff) 
                {
                    i_position prior_c;
                    if (b > 0)
                    {
                        ulint c_pos = prior_last_c_map[kv.first];
                        ulint c_b = intervals[c_pos];
                        ulint c_off = (lens[c_pos] - 1) + offsets[c_pos];

                        while (c_off >= lens[c_b]) 
                        {
                            c_off -= lens[c_b];
                            ++c_b;
                        }

                        prior_c = i_position{c_b, c_off};
                    }

                    switch(kv.first)
                    {
                        case 'A':
                            curr.A_map = block_c_map['A'];
                            curr.A_bv = bv(kv.second);
                            curr.A_diff = bv_select_1(&curr.A_bv);
                            curr.prior_A_LF = prior_c;
                            break;

                        case 'C':
                            curr.C_map = block_c_map['C'];
                            curr.C_bv = bv(kv.second);
                            curr.C_diff = bv_select_1(&curr.C_bv);
                            curr.prior_C_LF = prior_c;
                            break;

                        case 'G':
                            curr.G_map = block_c_map['G'];
                            curr.G_bv = bv(kv.second);
                            curr.G_diff = bv_select_1(&curr.G_bv);
                            curr.prior_G_LF = prior_c;
                            break;

                        case 'T':
                            curr.T_map = block_c_map['T'];
                            curr.T_bv = bv(kv.second);
                            curr.T_diff = bv_select_1(&curr.T_bv);
                            curr.prior_T_LF = prior_c;
                            break;

                        default:
                            curr.else_map.insert(std::pair<char, ulint>(kv.first, block_c_map[kv.first]));
                            curr.else_bv.insert(std::pair<char, bit_vec>(kv.first, bv(kv.second)));
                            curr.else_diff.insert(std::pair<char, bv_select_1>(kv.first, bv_select_1(&curr.else_bv[kv.first])));
                            curr.else_prior_LF.insert(std::pair<char, i_position>(c, prior_c));
                            break;
                    }
                }

                curr.lengths = dac_vec(block_lens);
                curr.offsets = dac_vec(block_offsets);
                
                curr.idx = block_idx;
                block_idx = next_idx;
                
                block_chars = vector<char>(block_size);
                block_lens = vector<ulint>(block_size);
                block_offsets = vector<ulint>(block_size);
                block_c_map = std::unordered_map<char, ulint>();
                prior_last_c_map = last_c_map;
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
     * LF step from the posision given
     * \param i_position of Interval-BWT
     * \return i_position of preceding character at pos
     */
    i_position LF(i_position pos)
    {
        assert(pos.run < r);

        i_block* curr = &B_table[pos.run/block_size];
        const ulint k = pos.run%block_size;

        const auto [d, c] = curr->heads.inverse_select(k);
        ulint q = curr->get_interval(c, d);
        ulint offset = pos.offset + curr->offsets[k];

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
                if (b == B_table.size()-1)
                {
                    return range_t(i_position{1,0}, i_position{0,0});
                }
                else
                {
                    start = curr.get_next_LF(c);
                }
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
                if (b == 0)
                {
                    return range_t(i_position{1,0}, i_position{0,0});
                }
                else
                {
                    end = curr.get_prior_LF(c);
                }
            }
            else
            {
                k = curr->heads.select(c_rank_s, c);
                offset = curr->lengths[k] - 1;

                end = LF_from_rank(curr, k, offset, c_rank_s, c);
            }
        }
        else
        {
            end = end = LF_from_rank(curr, k, offset, c_rank_s, c);
        }

        return range_t(start, end);
    }

    char get_char(ulint run) 
    {
        return (char) B_table[run/block_size].heads[run%block_size];
    }

    char get_char(i_position pos) 
    {
        return get_char(pos.run);
    }

    range_t full_range()
    {
        i_position first = {0, 0};
        i_position second = {(r-1), B_table[(r-1)/block_size].offsets[(r-1)%block_size]};
        return range_t(first, second);
    }

    ulint i_position_to_idx(i_position pos)
    {
        ulint pos_b = pos.run/block_size;
        ulint pos_k = pos.run%block_size;

        i_block* block = &B_table[pos_b];
        
        ulint idx;
        ulint k;
        // Take the stored idx of the start of the next block
        // Walk back to solution
        if (pos_k > block_size/2 && pos_b != B_table.size() - 1)
        {
            idx = B_table[pos_b+1].idx;
            k = block_size - 1;
            while (k >= pos_k)
            {
                idx -= block->lengths[k];
                --k;
            }
        }
        // Take the stored idx of the current block
        // Walk forward to solution
        else
        {
            idx = block->idx;
            k = 0;
            while (k < pos_k)
            {
                idx += block->lengths[k];
                ++k;
            }
        }
        idx += pos.offset;
        return idx;
    }

    // ulint idx_to_i_position(i_position pos)
    // {
    //
    // }

    ulint number_of_runs()
    {
        return r;
    }

    ulint size()
    {
        return n;
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

        out.write((char *)&n, sizeof(n));
        written_bytes += sizeof(n);

        out.write((char *)&r, sizeof(r));
        written_bytes += sizeof(r);

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
        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));
        in.read((char *)&size, sizeof(size));
        B_table = std::vector<i_block>(size);
        for(size_t i = 0; i < size; ++i)
        {
            B_table[i].load(in);
        }
    }

private:
    ulint n;
    ulint r;
    vector<i_block> B_table;

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
    }
};

#endif /* end of include guard: _R_INDEX_F_HH */
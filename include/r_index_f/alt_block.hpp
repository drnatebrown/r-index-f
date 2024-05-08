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

#ifndef _ALT_BLOCK_HH
#define _ALT_BLOCK_HH

#include <sdsl/dac_vector.hpp>
#include <common.hpp>

#include <ds/heads_wt_w.hpp>
#include <ds/heads_bv_w.hpp>
#include <ds/intervals_rank_w.hpp>
#include <ds/interval_pos.hpp>

#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

constexpr uchar A = 0b00;
constexpr uchar C = 0b01;
constexpr uchar G = 0b10;
constexpr uchar T = 0b11;

using namespace std;

template  < ulint block_size = 256 >
class alt_block
{
public:
    alt_block() {}

    alt_block(std::ifstream &bwt) {}

    alt_block(std::ifstream &heads, std::ifstream &lengths)
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
            lens.push_back((ulint) length);
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

        blocks = std::vector<block>(num_blocks()); // use ::move
        for (size_t i = 0; i < (r / block_size) * block_size; i += block_size) {
            block& b = get_block(i);
            b.run_heads = run_heads_t(std::vector<uchar>(chars.begin() + i, chars.begin() + i + block_size));
            b.dest_off = offsets_t(std::vector<ulint>(offs.begin() + i, offs.begin() + i + block_size));
            b.dest_pred = dest_pred_t(std::vector<uchar>(chars.begin() + i, chars.begin() + i + block_size), std::vector<ulint>(lf_dest.begin() + i, lf_dest.begin() + i + block_size));
            b.run_len = lengths_t(std::vector<ulint>(lens.begin() + i, lens.begin() + i + block_size));
        }
        if (r % block_size != 0) {
            size_t p = (r / block_size) * block_size;
            block& b = get_block(p);
            b.run_heads = run_heads_t(std::vector<uchar>(chars.begin() + p, chars.end()));
            b.dest_off = offsets_t(std::vector<ulint>(offs.begin() + p, offs.end()));
            b.dest_pred = dest_pred_t(std::vector<uchar>(chars.begin() + p, chars.end()), std::vector<ulint>(lf_dest.begin() + p, lf_dest.end()));
            b.run_len = lengths_t(std::vector<ulint>(lens.begin() + p, lens.end()));
        }
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
        return run_heads(i);
    }

    uchar get_char(interval_pos pos) 
    {
        return run_heads(pos.run);
    }

    interval_pos begin()
    {
        return interval_pos(0, 0);
    }

    interval_pos end()
    {
        return interval_pos(r-1, run_len(r-1) - 1);
    }

    interval_pos LF(interval_pos pos)
    {
        auto [c_rank, c] = get_block(pos.run).run_heads.inverse_select(row(pos.run));
        return LF(pos.run, pos.offset, c_rank, c);
    }

    interval_pos LF_prior(interval_pos pos, uchar c)
    {
        ulint curr_run = pos.run;
        block& b = get_block(curr_run);
        ulint c_rank = b.run_heads.rank(row(curr_run) + 1, c);
        while(c_rank == 0 && (curr_run / block_size) > 0) {
            curr_run = ((curr_run/block_size) * block_size) - 1;
            b = get_block(curr_run);
            c_rank = b.run_heads.rank(row(curr_run) + 1, c);
        }
        if (c_rank == 0) {
            return interval_pos();
        }
        c_rank -= 1;

        ulint prior_run = first_block_run(curr_run) + b.run_heads.select(c_rank + 1, c);
        ulint prior_off = (pos.run != prior_run) ? run_len(prior_run) - 1 : pos.offset;
        
        return LF(prior_run, prior_off, c, c_rank);
    }

    interval_pos LF_next(interval_pos pos, uchar c)
    {
        ulint curr_run = pos.run;
        ulint c_rank = get_block(curr_run).run_heads.rank(row(curr_run), c);
        while(c_rank + 1 > get_block(curr_run).run_heads.rank(block_size, c) && (curr_run / block_size) < blocks.size() - 1) {
            curr_run = ((curr_run/block_size) + 1) * block_size;
            c_rank = get_block(curr_run).run_heads.rank(row(curr_run) + 1, c);
        }
        if (c_rank + 1 > get_block(curr_run).run_heads.rank(block_size, c)) return interval_pos();

        ulint next_run = first_block_run(curr_run) + get_block(curr_run).run_heads.select(c_rank + 1, c);
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

        while (off >= run_len(run)) 
        {
            off -= run_len(run++);
        }

        return interval_pos(run, off);
    }

    // For a general interval position, return the idx wrt. the BWT
    ulint interval_to_idx(interval_pos pos)
    {
        ulint curr = 0;
        ulint idx = 0;
        while(curr < pos.run) {
            idx += run_len(curr++);
        }
        return idx + pos.offset;
    }

    // For a general index on the BWT, return the corresponding interval position
    interval_pos idx_to_interval(ulint idx)
    {
        ulint curr = 0;
        ulint curr_idx = 0;
        while(idx - curr_idx > run_len(curr)) {
            curr_idx += run_len(curr++);
        }
        return interval_pos(curr, idx - curr_idx);
    }

    // ulint distance(interval_pos p1, interval_pos p2) {

    // }

    std::string get_file_extension() const
    {
        return ".alt_block";
    }

    void mem_stats()
    {
        sdsl::nullstream ns;

        size_t rh_bytes = 0;
        size_t rl_bytes = 0;
        size_t dp_bytes = 0;
        size_t do_bytes = 0;
        for (size_t i = 0; i < blocks.size(); ++i) {
            rh_bytes += blocks[i].run_heads.serialize(ns);
            rl_bytes += blocks[i].run_len.serialize(ns);
            dp_bytes += blocks[i].dest_pred.serialize(ns);
            do_bytes += blocks[i].dest_off.serialize(ns);
        }
        verbose("Blocks: ", num_blocks());
        verbose("Memory consumption (bytes).");
        verbose("  Columnar LF table: ", serialize(ns));
        verbose("         Mean block: ", (rh_bytes + rl_bytes + dp_bytes + do_bytes)/blocks.size());
        verbose("                                 ");
        verbose("              Run_heads: ", rh_bytes);
        verbose("                Run_len: ", rl_bytes);
        verbose("              Dest_pred: ", dp_bytes);
        verbose("               Dest_off: ", do_bytes);
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

        for(size_t i = 0; i < blocks.size(); ++i)
        {
           written_bytes += blocks[i].serialize(out,v,"block_table_" + std::to_string(i));
        }

        return written_bytes;
    }

    /* load from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));

        blocks = std::vector<block>(num_blocks());
        for(size_t i = 0; i < blocks.size(); ++i)
        {
            blocks[i].load(in);
        }
    }

private:
    ulint n; // Length of BWT
    ulint r; // Runs of BWT

    // struct vec_heads {
    //     std::vector<uchar> data;
    //     size_t cursor = 0;
    //     size_t last_rank = 0;
    //     uchar last_c = A; // Default to A

    //     vec_heads() : data() {}

    //     vec_heads(const std::vector<uchar>& input) {
    //         data.reserve((input.size() / 4) + (input.size() % 4 != 0)); // Reserve space for packed 2-bit representations
    //         uchar packed = 0;
    //         int count = 0;
    //         for (char c : input) {
    //             packed |= map_to_2bit(c) << (count * 2);
    //             if (++count == 4) {
    //                 data.push_back(packed);
    //                 packed = 0;
    //                 count = 0;
    //             }
    //         }
    //         if (count > 0) {
    //             // Pad with A if the last byte is incomplete
    //             packed |= A << (count * 2);
    //             data.push_back(packed);
    //         }
    //     }

    //     uchar map_to_2bit(char c) {
    //         switch (c) {
    //             case 'A': return A;
    //             case 'C': return C;
    //             case 'G': return G;
    //             case 'T': return T;
    //             default: return A; // Default to A for other characters
    //         }
    //     }

    //     size_t rank(char c, size_t pos) {
    //         size_t count = 0;
    //         size_t i;
    //         for (i = 0; i < pos; ++i) {
    //             if (get_char(i) == map_to_2bit(c)) {
    //                 count++;
    //             }
    //         }
    //         cursor = i;
    //         last_rank = count;
    //         last_c = map_to_2bit(c);
    //         return count;
    //     }

    //     size_t select(char c, size_t i) {
    //         size_t count = 0;
    //         for (size_t pos = cursor; pos < data.size() * 4; ++pos) {
    //             if (get_char(pos) == map_to_2bit(c)) {
    //                 count++;
    //                 if (count == i) {
    //                     return pos;
    //                 }
    //             }
    //         }
    //         return data.size() * 4; // Return size of data if not found
    //     }

    //     std::pair<size_t, char> inverse_select(size_t i) {
    //         uchar packed = data[i / 4];
    //         size_t pos_in_byte = (i % 4) * 2;
    //         uchar c = (packed >> pos_in_byte) & 0b11;
    //         size_t rank = rank_from_start(c, i);
    //         return {rank, char_from_2bit(c)};
    //     }

    //     char operator[](size_t i) {
    //         return char_from_2bit(get_char(i));
    //     }

    //     size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") {
    //         sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    //         size_t written_bytes = 0;

    //         // Serialize the data size
    //         size_t dataSize = data.size();

    //         // Serialize the data
    //         out.write(reinterpret_cast<char*>(data.data()), dataSize * sizeof(uchar));
    //         written_bytes += dataSize * sizeof(uchar);

    //         return written_bytes;
    //     }

    //     void load(std::istream &in, size_t r, size_t curr_i, size_t bs) {
    //         size_t dataSize = (r % bs != 0 && curr_i == (r / bs)) ? ((r / bs) + 1)*bs - r : bs;
    //         dataSize = (dataSize / 4) + (dataSize % 4 != 0);
    //         data.resize(dataSize);
    //         in.read(reinterpret_cast<char*>(data.data()), dataSize * sizeof(uchar));
    //     }

    // private:
    //     size_t rank_from_start(uchar c, size_t i) {
    //         size_t count = 0;
    //         for (size_t pos = 0; pos <= i; ++pos) {
    //             if (get_char(pos) == c) {
    //                 count++;
    //             }
    //         }
    //         return count;
    //     }

    //     char char_from_2bit(uchar c) {
    //         switch (c) {
    //             case A: return 'A';
    //             case C: return 'C';
    //             case G: return 'G';
    //             case T: return 'T';
    //             default: return 'A'; // Default to A for unknown characters
    //         }
    //     }

    //     uchar get_char(size_t pos) {
    //         size_t byte_index = pos / 4;
    //         size_t bit_offset = (pos % 4) * 2;
    //         return (data[byte_index] >> bit_offset) & 0b11;
    //     }
    // };

    // typedef struct vec_heads
    // {
    //     std::vector<uchar> data;
    //     size_t cursor = 0;
    //     size_t last_rank = 0;
    //     uchar last_c = 0;
    //     vec_heads() : data() {}

    //     vec_heads(const std::vector<uchar>& input) : data(input) {}

    //     size_t rank(uchar c, size_t pos) {
    //         size_t count = 0;
    //         size_t i;
    //         for (i = 0; i < pos; ++i) {
    //             if (data[i] == c) {
    //                 count++;
    //             }
    //         }
    //         cursor = i;
    //         last_rank = count;
    //         last_c = c;
    //         return count;
    //     }

    //     // Select function: finds the position of the ith character c
    //     size_t select(uchar c, size_t i) {
    //         size_t count = 0;
    //         if (c == last_c && i > last_rank) {
    //             count = last_rank;
    //         }
    //         else {
    //             cursor = 0;
    //         }
    //         for (size_t pos = cursor; pos < data.size(); ++pos) {
    //             if (data[pos] == c) {
    //                 count++;
    //                 if (count == i) {
    //                     return pos;
    //                 }
    //             }
    //         }
    //         return data.size(); // Return size of data if not found
    //     }

    //     // Inverse Select function: returns the rank at position i and its character
    //     std::pair<size_t, uchar> inverse_select(size_t i) {
    //         uchar c = data[i];
    //         size_t rnk = rank(c, i);
    //         return {rnk, c};
    //     }

    //     // Access function: returns the character at position i
    //     uchar operator[](size_t i) {
    //         return data[i];
    //     }

    //     size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") {
    //         sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    //         size_t written_bytes = 0;

    //         // Serialize the data size
    //         size_t dataSize = data.size();

    //         // Serialize the data
    //         out.write(reinterpret_cast<char*>(data.data()), dataSize * sizeof(uchar));
    //         written_bytes += dataSize * sizeof(uchar);

    //         return written_bytes;
    //     }

    //     void load(std::istream &in, size_t r, size_t curr_i, size_t bs) {
    //         size_t dataSize = (r % bs != 0 && curr_i == (r / bs)) ? ((r / bs) + 1)*bs - r : bs;
    //         data.resize(dataSize);
    //         in.read(reinterpret_cast<char*>(data.data()), dataSize * sizeof(uchar));
    //     }
    // };

    // typedef heads_wt_w<> run_heads_t; // Huffman-Shaped WT
    // typedef vec_heads run_heads_t;
    typedef heads_bv_w<bit_vector, symbol_map> run_heads_t;
    typedef intervals_rank_w<base_bv<>, symbol_map> dest_pred_t; // Plain Bitvector per Character
    typedef sdsl::dac_vector_dp<> lengths_t; // Sparse Bitvector per Character
    typedef sdsl::dac_vector_dp<> offsets_t; // Sparse Bitvector

    typedef struct block
    {
        run_heads_t run_heads;
        lengths_t run_len;
        offsets_t dest_off;
        dest_pred_t dest_pred;

        size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name ="")
        {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_t written_bytes = 0;

            written_bytes += run_heads.serialize(out, v, "Run_heads");
            written_bytes += run_len.serialize(out, v, "Lengths");
            written_bytes += dest_pred.serialize(out, v, "Dest_pred");
            written_bytes += dest_off.serialize(out, v, "Offsets");

            return written_bytes;
        }

        // void load(std::istream &in, size_t r, size_t i, size_t bs)
        // {
        //     run_heads.load(in, r, i, bs);
        //     run_len.load(in);
        //     dest_pred.load(in);
        //     dest_off.load(in);
        // }

        void load(std::istream &in)
        {
            run_heads.load(in);
            run_len.load(in);
            dest_pred.load(in);
            dest_off.load(in);
        }
    };

    std::vector<block> blocks;

    ulint row(ulint run)
    {
        return run % block_size;
    }

    block& get_block(ulint run)
    {
        assert(run < r);
        return blocks[run / block_size];
    }

    ulint first_block_run(ulint run)
    {
        return (run / block_size) * block_size;
    }

    uchar run_heads (ulint i) {
        return get_block(i).run_heads[row(i)];
    }

    ulint run_len (ulint i) {
        return get_block(i).run_len[row(i)];
    }
    
    ulint dest_pred (ulint i, uchar c, ulint c_rank) {
        return get_block(i).dest_pred.get(c, c_rank);
    }

    ulint dest_off (ulint i) {
        return get_block(i).dest_off[row(i)];
    }

    ulint num_blocks() {
        return r / block_size + ((r % block_size) != 0);
    }

    interval_pos LF(ulint k, ulint d, uchar c, ulint c_rank) {
        ulint next_interval = dest_pred(k, c, c_rank);
        ulint next_off = dest_off(k) + d;
        return reduced_pos(interval_pos(next_interval, next_off));
    }
};

#endif /* end of include guard: _LF_TABLE_HH */
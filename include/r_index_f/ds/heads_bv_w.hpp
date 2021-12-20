/* heads_bv_w - Wrapper to store the heads in full bit vectors (access is bad, trying all bit vectors)
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
   \file heads_bv_w.hpp
   \brief heads_bv_w Wrapper to store the heads in a wavelet tree
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _HEADS_BV_W_HH
#define _HEADS_BV_W_HH

#include <common.hpp>

#include <ds/ACGT_map.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace sdsl;

template< class bv_t = bit_vector,
          template<class> class char_map_t = ACGT_map >
class heads_bv_w
{
private:
    typedef typename bv_t::select_1_type bv_select_1;
    typedef typename bv_t::rank_1_type bv_rank_1;

    struct rank_select_bv {
        bv_t bv;
        bv_select_1 select;
        bv_rank_1 rank;

        rank_select_bv() {}

        rank_select_bv(ulint size)
        {
            bv = bv_t(size, false);
        }

        size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
        {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_t written_bytes = 0;

            written_bytes += bv.serialize(out, v, "bv");

            return written_bytes;
        }

        void load(std::istream &in)
        {
            bv.load(in);
            select = bv_select_1(&bv);
            rank = bv_rank_1(&bv);
        }
    };

    typedef char_map_t<rank_select_bv> bv_map;

    bv_map bit_vecs;
    ulint bv_size;

    uchar scan(ulint idx)
    {
        for (size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(bit_vecs.contains(i))
            {
                if (bit_vecs[i].bv[idx])
                {
                    return i;
                }
            }
        }

        return 0;
    }

public:
    heads_bv_w() {}

    heads_bv_w(std::vector<uchar> chars) {
        bit_vecs = bv_map();
        bv_size = chars.size();

        for (size_t i = 0; i < chars.size(); ++i)
        {
            uchar c = chars[i];
            if (!bit_vecs.contains(c))
            {
                bit_vecs.insert(std::pair<uchar, rank_select_bv>(c, rank_select_bv(bv_size)));
            }

            if (bit_vecs.contains(c))
            {
                bit_vecs[c].bv[i] = true;
            }
        }

        for (size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(bit_vecs.contains(i))
            {
                bit_vecs[i].select = bv_select_1(&bit_vecs[i].bv);
                bit_vecs[i].rank = bv_rank_1(&bit_vecs[i].bv);
            }
        }
    }

    ulint rank(ulint idx, uchar c)
    {
        return bit_vecs[c].rank(idx);
    }

    ulint select(ulint idx, uchar c)
    {
        return bit_vecs[c].select(idx);
    }

    std::pair<ulint, uchar> inverse_select(ulint idx)
    {
        uchar c = scan(idx);
        ulint r = 0;
        if (bit_vecs.contains(c))
        {
            r = rank(idx, c);
        }
        return std::pair<ulint, uchar>(r, c);
    }

    uchar operator[](size_t idx) {
        scan(idx);
    }

    ulint size()
    {
        return bv_size;;
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += bit_vecs.serialize(out, v, "symbols");
        out.write((char *)&bv_size, sizeof(bv_size));
        written_bytes += sizeof(bv_size);

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        bit_vecs.load(in);
        in.read((char *)&bv_size, sizeof(bv_size));
    }
};

#endif /* end of include guard: _HEADS_WT_W_HH */
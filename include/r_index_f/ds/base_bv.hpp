/* base_bv - Holds interval as base and diff computed using bit vector
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
   \file base__bv.hpp
   \brief base_bv.hpp Returns intervals by calculating a difference from a stored base (intervals non-decreasing sequence wrt. character)
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _BASE_BV_HH
#define _BASE_BV_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace sdsl;

template<class bv_t = bit_vector >
class base_bv
{
private:
    typedef typename bv_t::select_1_type bv_select_1;

    ulint base;
    bv_t diff_bv;
    bv_select_1 diff_select;

public:

    base_bv() {}

    base_bv(ulint b, std::vector<ulint> diffs) 
    {
        std::vector<bool> bit_diff = std::vector<bool>();
        for(size_t i = 0; i < diffs.size(); ++i)
        {
            ulint diff = diffs[i];
            while (diff > 0) {
                bit_diff.push_back(false);
                --diff;
            }
            bit_diff.push_back(true);
        }

        base = b;
        diff_bv = bool_to_bit_vec<bv_t>(bit_diff);
        diff_select = bv_select_1(&diff_bv);
    }

    ulint get(ulint rank) const
    {
        return base + diff_select(rank+1) - rank;
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        out.write((char *)&base, sizeof(base));
        written_bytes += sizeof(base);

        written_bytes += diff_bv.serialize(out, v, "diff_bv");

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        in.read((char *)&base, sizeof(base));
        diff_bv.load(in);
        diff_select= bv_select_1(&diff_bv);
    }
};

#endif /* end of include guard: _BASE_BV_HH */
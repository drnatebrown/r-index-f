/* heads_wt_w - Wrapper to store the heads in a wavelet tree
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
   \file heads_wt_w.hpp
   \brief heads_wt_w Wrapper to store the heads in a wavelet tree
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _HEADS_WT_W_HH
#define _HEADS_WT_W_HH

#include <common.hpp>

#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace sdsl;

template< class wt_t = wt_huff<bit_vector> >
class heads_wt_w
{
private:
wt_t symbols;

public:
    heads_wt_w() {}

    heads_wt_w(std::vector<char> chars) {
        construct_im(symbols, std::string(chars.begin(), chars.end()).c_str(), 1);
    }

    ulint rank(ulint idx, char c)
    {
        return symbols.rank(idx, c);
    }

    ulint select(ulint idx, char c)
    {
        return symbols.select(idx, c);
    }

    std::pair<ulint, ulint> inverse_select(ulint idx)
    {
        return symbols.inverse_select(idx);
    }

    const char& operator[](const ulint idx) const {
        return symbols[idx];
    }

    ulint size()
    {
        return symbols.size();
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        symbols.serialize(out, v, "symbols");

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        symbols.load(in);
    }
};

#endif /* end of include guard: _HEADS_WT_W_HH */
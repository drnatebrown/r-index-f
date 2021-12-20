/* idx_bit_vector.hpp - Sampling idx using bit vector approach
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
   \file idx_bit_vector.hpp
   \brief idx_bit_vector.hpp template class wrapper used to access idx sampling using bit vector approach
   \author Nathaniel Brown
   \date 16/12/2021
*/

#ifndef _IDX_BV_HH
#define _IDX_BV_HH

#include <common.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

template < class bit_vec = sd_vector<> >
class idx_bit_vector
{
private:
    typedef typename bit_vec::rank_1_type idx_rank;
    typedef typename bit_vec::select_1_type idx_select;

    bit_vec samples;
    idx_rank pred;
    idx_select run_sample;

public:

    idx_bit_vector() {}

    idx_bit_vector(vector<bool> vec) {
        samples = bool_to_bit_vec<bit_vec>(vec);
        pred = idx_rank(&samples);
        run_sample = idx_select(&samples);
    }

    ulint sample(ulint rank)
    {
        return run_sample(rank + 1);
    }

    ulint predecessor(ulint idx) {
        return pred(idx + 1);
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += samples.serialize(out, v, "idx_bit_vec");

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        samples.load(in);
        pred = idx_rank(&samples);
        run_sample = idx_select(&samples);
    }
};

#endif /* end of include guard: _IDX_BV_HH */
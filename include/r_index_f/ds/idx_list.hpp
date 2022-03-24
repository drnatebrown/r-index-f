/* idx_list.hpp - Sampling idx using bexplicit list (vector)
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
   \file idx_list.hpp
   \brief idx_list.hpp template class wrapper used to access idx sampling using list approach
   \author Nathaniel Brown
   \date 16/12/2021
*/

#ifndef _IDX_LIST_HH
#define _IDX_LIST_HH

#include <common.hpp>

#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

class idx_list
{
private:
    std::vector<ulint> samples;

public:

    idx_list() {}

    idx_list(std::vector<bool> vec) {
        samples = std::vector<ulint>();

        for(size_t i = 0; i < vec.size(); ++i)
            if (vec[i]) samples.push_back(i);
    }

    ulint sample(ulint rank)
    {
        assert(rank < samples.size());
        return samples[rank];
    }

    ulint predecessor(ulint idx) {
        // Get first element equal to or greater than idx (runs are sorted, so O(lg n) using binary search)
        auto pred = std::lower_bound(samples.begin(), samples.end(), idx);
        if(*pred != idx)
        {
            // Index in sampling array of predecessor (minus 1, since it is first element greater)
            pred -= 1;
        }

        std::distance(samples.begin(), pred);
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        size_t size = samples.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        for(size_t i = 0; i < size; ++i)
        {
            out.write((char *)&samples[i], sizeof(samples[i]));
            written_bytes += sizeof(samples[i]);
        }

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        size_t size;
        in.read((char *)&size, sizeof(size));
        samples = std::vector<ulint>(size);
        for(size_t i = 0; i < size; ++i)
        {
            in.read((char *)&samples[i], sizeof(samples[i]));
        }
    }
};

#endif /* end of include guard: _IDX_LIST_HH */
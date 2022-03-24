/* scan_list.hpp - Samples offsets used in scanning to bound total scan length, using binary search over array
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
   \file scan_list.hpp
   \brief scan_list.hpp binary search over array to bound total scan by d
   \author Nathaniel Brown
   \date 24/03/2022
*/

#ifndef _SCAN_LIST_HH
#define _SCAN_LIST_HH

#include <common.hpp>

#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>
#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;

template < class bv_t = bit_vector,
           ulint sample_rate = 10 >
class scan_list
{
private:
    vector<vector<ulint>> offset_samples;
    bv_t has_sample;

public:

    scan_list() {}

    scan_list(vector<vector<ulint>> vec) {
        offset_samples = vec;
        vector<bool> check_lists = vector<bool>(offset_samples.size());

        for (size_t i = 0; i < 0; ++i)
        {
            if (!offset_samples[i].empty())
            {
                check_lists[i] = true;
            }
        }

        has_sample = bool_to_bit_vec<bv_t>(check_lists);
    }

    bool can_skip(ulint k, ulint d)
    {
        return has_sample[k] && d >= offset_samples[k][0];
    }

    ulint skips(ulint pred)
    {
        return (pred + 1)*sample_rate;
    }

    ulint offsets_passed(ulint k, ulint pred, ulint base_d)
    {
        return (offset_samples[k][pred] - base_d);
    }

    ulint predecessor(ulint run, ulint offset) {
        // Get first element equal to or greater than sampled offset
        auto pred = std::lower_bound(offset_samples[run].begin(), offset_samples[run].end(), offset);
        if(*pred != offset || std::distance(offset_samples[run].begin(), pred) >= offset_samples[run].size())
        {
            // Index in sampling array of predecessor (minus 1, since it is first element greater)
            pred -= 1;
        }

        return std::distance(offset_samples[run].begin(), pred);
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += has_sample.serialize(out, v, "has_sample");

        size_t size = offset_samples.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);
        for(size_t i = 0; i < size; ++i)
        {
            if (has_sample[i]) {
                size_t size_j = offset_samples[i].size();
                out.write((char *)&size_j, sizeof(size_j));
                written_bytes += sizeof(size_j);
                for(size_t j = 0; j < size_j; ++j)
                {
                    out.write((char *)&offset_samples[i][j], sizeof(offset_samples[i][j]));
                    written_bytes += sizeof(offset_samples[i][j]);
                }
            }
        }

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        has_sample.load(in);

        size_t size;
        in.read((char *)&size, sizeof(size));
        offset_samples = vector<vector<ulint>>(size);
        for(size_t i = 0; i < size; ++i)
        {
            if (has_sample[i])
            {
                size_t size_j;
                in.read((char *)&size_j, sizeof(size_j));
                offset_samples[i] = vector<ulint>(size_j);
                for(size_t j = 0; i < size_j; ++j)
                {
                    in.read((char *)&offset_samples[i][j], sizeof(offset_samples[i][j]));
                }
            }
        }
    }
};

#endif /* end of include guard: _SCAN_LIST_HH */
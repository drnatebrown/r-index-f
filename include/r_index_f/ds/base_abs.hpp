/* base_sample - Holds interval as base and diff computed using dac and sampled diffs
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
   \file base_sample.hpp
   \brief base_sample.hpp Returns interval by computing from a base and diff, sampling absolute diff positions
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _BASE_SAMPLE_HH
#define _BASE_SAMPLE_HH

#include <common.hpp>

#include <sdsl/dac_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace sdsl;

template < ulint sample_rate = 10,
           class vec_t = dac_vector<> >
class base_sample
{
private:
    ulint base;
    std::vector<ulint> sampled_diffs;
    vec_t partial_diffs;

    ulint get_diff(ulint rank) const
    {
        ulint sample_rank = rank / sample_rate;
        ulint sample_next = sample_rank*sample_rate + 1;
        ulint diff = sampled_diffs[sample_rank];
        while (sample_next <= rank)
        {
            diff += partial_diffs[sample_next++];
        }

        return diff;
    }

public:

    base_abs() {}

    base_abs(ulint b, std::vector<ulint> diffs) 
    {
        base = b;

        ulint absolute_diff = 0;
        
        for(size_t i = 0; i < diffs.size(); ++i)
        {
            absolute_diff += diffs[i];
            if (i % sample_rate == 0)
            {
                sampled_diffs.push_back(absolute_diff);
            }
        }

        partial_diffs = vec_t(diffs);
    }

    ulint get(ulint rank) const
    {
        return base + get_diff(rank);
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

        size_t size = sampled_diffs.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);
        for(size_t i = 0; i < size; ++i)
        {
            out.write((char *)&sampled_diffs[i], sizeof(sampled_diffs[i]));
            written_bytes += sizeof(sampled_diffs[i]);
        }

        written_bytes += partial_diffs.serialize(out, v, "partial_diffs");

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        in.read((char *)&base, sizeof(base));
       
        size_t size;
        in.read((char *)&size, sizeof(size));
        sampled_diffs = std::vector<ulint>(size);
        for(size_t i = 0; i < size; ++i)
        {
            in.read((char *)&sampled_diffs[i], sizeof(sampled_diffs[i]));
        }

        partial_diffs.load(in);
    }
};

#endif /* end of include guard: _BASE_SAMPLE_HH */
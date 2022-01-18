/* base_interpolate - Holds interval as base and diff computed using dac/interpolative coding
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
   \file base_interpolate.hpp
   \brief base_interpolate.hpp Returns interval by computing from a base and diff retrieved from interpolation, sampling absolute diff positions
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _BASE_INTERPOLATE_HH
#define _BASE_INTERPOLATE_HH

#include <common.hpp>

#include <sdsl/dac_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace sdsl;

template < ulint sample_rate = 16,
           class bv_t = bit_vector,
           class vec_t = dac_vector_dp<> >
class base_interpolate
{
private:
    ulint base;
    std::vector<ulint> sampled_diffs;
    bv_t interp_neg;
    vec_t interp_diffs;

    ulint get_diff(ulint rank) const
    {
        // If value is sampled (at sample rate or last value) then return it
        if (rank % sample_rate == 0) return sampled_diffs[rank/sample_rate];
        if (rank == interp_diffs.size() - 1) return sampled_diffs[sampled_diffs.size() - 1];

        // Find the last and next samples for the given position
        ulint x = sampled_diffs[rank/sample_rate];
        ulint z = sampled_diffs[rank/sample_rate+1];

        int sign = (interp_neg[rank]) ? -1 : 1;
        // Add weighted sum of samples to stored encoding, which returns y
        return sign*interp_diffs[rank] + (x + ((z - x)*(rank - sample_rate*(rank/sample_rate))/sample_rate));
    }

public:

    base_interpolate() {}

    base_interpolate(ulint b, std::vector<ulint> diffs) 
    {
        base = b;

        ulint absolute_diff = 0;
        vector<long int> encoding = vector<long int>(diffs.size(), 0);
        interp_neg = bv_t(diffs.size(), 0);
        vector<ulint> full_diffs = vector<ulint>(diffs.size());
        for(size_t i = 0; i < diffs.size(); ++i)
        {
            absolute_diff += diffs[i];
            full_diffs[i] = absolute_diff;
            if (i % sample_rate == 0 || i == diffs.size() - 1)
            {
                sampled_diffs.push_back(absolute_diff);
                if (i != 0)
                {
                    // Last sample
                    ulint x = sampled_diffs[sampled_diffs.size() - 2];
                    // Next sample
                    ulint z = absolute_diff;

                    // If not at last entry, the last sampled is always a distance sample_rate away, otherwise the distance past the last sample
                    ulint last_dist = i % sample_rate;
                    if (last_dist == 0) last_dist = sample_rate;

                    // Iterate between last and next sample
                    for (ulint j = i - last_dist + 1; j < i; ++j)
                    {
                        // current diff
                        ulint y = full_diffs[j];
                        // Get the difference from the weighted average of the next/last samples and the current diff
                        long int encode = y - (x + ((z - x)*(j - sample_rate*(j/sample_rate))/sample_rate));
                        if (encode < 0) {
                            interp_neg[j] = true;
                            encode*=-1;
                        }

                        encoding[j] = encode;
                    }
                }
            }
        }

        interp_diffs = vec_t(encoding);
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

        written_bytes += interp_neg.serialize(out, v, "interp_neg");
        written_bytes += interp_diffs.serialize(out, v, "interp_diffs");

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

        interp_neg.load(in);
        interp_diffs.load(in);
    }
};

#endif /* end of include guard: _BASE_INTERPOLATIVE_HH */
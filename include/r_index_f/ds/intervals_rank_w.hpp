/* intervals_rank_w - Wrapper which accesses intervals by seperating their mapping wrt. character, i.e. compute given character rank
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
   \file intervals_rank_w.hpp
   \brief intervals_rank_w Returns intervals for respective character and rank of that character
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _INTERVALS_RANK_W_HH
#define _INTERVALS_RANK_W_HH

#include <common.hpp>
#include <ds/base_bv.hpp>
#include <ds/base_sample.hpp>
#include <ds/symbol_map.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

#include <ds/ACGT_map.hpp>

using namespace sdsl;

template< class interval_t = base_bv<>,
          template<class> class char_map_t = ACGT_map >
class intervals_rank_w
{
private:
    typedef char_map_t<interval_t> interval_map;

    interval_map char_map;

public:
    intervals_rank_w() {}

    intervals_rank_w(std::vector<uchar> characters, std::vector<ulint> intervals) {
        assert(characters.size() == intervals.size());

        // Concerned with Interval sectioned by character (break into base pointer and difference from prior interval)
        std::unordered_map<uchar, ulint> last_c_map = std::unordered_map<uchar, ulint>();
        std::unordered_map<uchar, ulint> block_c_map = std::unordered_map<uchar, ulint>();
        std::unordered_map<uchar, std::vector<ulint>> diff = std::unordered_map<uchar, std::vector<ulint>>();

        for(size_t i = 0; i < characters.size(); ++i) 
        {
            uchar character = characters[i];
            ulint interval = intervals[i];

            if (!block_c_map.count(character)) {
                block_c_map[character] = interval;
                last_c_map[character] = interval;
                diff[character] = std::vector<ulint>();
            }

            diff[character].push_back(interval - last_c_map[character]);
            last_c_map[character] = interval;
        }

        for(size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(diff.count(i))
            {
                char_map.insert(std::pair<uchar, interval_t>(i, interval_t(block_c_map[i], diff[i])));
            }
        }
    }

    ulint get(uchar c, ulint c_rank)
    {
        return char_map[c].get(c_rank);
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += char_map.serialize(out, v, "char_map");

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        char_map.load(in);
    }
};

#endif /* end of include guard: _BASE_BV_HH */
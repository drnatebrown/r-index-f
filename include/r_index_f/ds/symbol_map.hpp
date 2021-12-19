/* symbol_map - Implements a simple map taking character (byte) positions
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
   \file symbol_map.hpp
   \brief symbol_map.hpp Bitvector for contains and vector for access
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _SYMBOL_MAP_HH
#define _SYMBOL_MAP_HH

#include <common.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

template<class T>
class symbol_map
{
private:
    bit_vector map_contains;
    std::vector<T> access;

public:

    symbol_map() {
        map_contains = bit_vector(ALPHABET_SIZE, false);
        access = std::vector<T>(ALPHABET_SIZE);
    }

    symbol_map(std::unordered_map<uchar, T> map) {
        map_contains = bit_vector(ALPHABET_SIZE, false);
        access = std::vector<T>(ALPHABET_SIZE);

        for(auto const& [c, val] : map)
        {
            map_contains[c] = true;
            access[c] = val;
        }
    }

    bool contains(uchar c) const
    {
        return map_contains[c];
    }

    bool insert(std::pair<uchar, T> kv)
    {
        if (!contains(kv.first)) {
            map_contains[kv.first] = true;
            access[kv.first] = kv.second;

            return true;
        }
        else {
            return false;
        }
    }

    T at(const uchar key)
    {
        if (!contains(key))
        {
            throw std::out_of_range("Symbol" + std::to_string(key) + "not in map");
        }

        return access[key];
    }

    const T& operator[](const uchar key ) const {
        return access[key];
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        map_contains.serialize(out, v, "char_next");
        for (size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if(contains(i))
            {
                written_bytes += access[i].serialize(out, v, "access_" + std::to_string(i)); 
            }
        }

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        map_contains.load(in);
        access = std::vector<T>(ALPHABET_SIZE);
        for(size_t i = 0; i < ALPHABET_SIZE; ++i)
        {
            if (contains(i))
            {
                access[i].load(in);
            }
        }
    }
};

#endif /* end of include guard: _SYMBOL_MAP_HH */
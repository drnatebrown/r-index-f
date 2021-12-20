/* ACGT_map Map which only accepts values for characters ACGT
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
   \file ACGT_map.hpp
   \brief ACGT_map.hpp maps only for characters ACGT
   \author Nathaniel Brown
   \date 18/12/2021
*/

#ifndef _ACGT_MAP_HH
#define _ACGT_MAP_HH

#include <common.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

#define ACGT_SIZE 4

template<class T>
class ACGT_map
{
private:
    bit_vector map_contains;
    T a_type;
    T c_type;
    T g_type;
    T t_type;

    T& find(uchar c)
    {
        switch(c)
        {
            case 'A': return a_type;
            case 'C': return c_type;
            case 'G': return g_type;
            case 'T': return t_type;
            default: throw std::out_of_range("Symbol" + std::to_string(c) + "not in map");
        }
    }

    bool allowed(uchar c)
    {
        switch(c)
        {
            case 'A': return true;
            case 'C': return true;
            case 'G': return true;
            case 'T': return true;
            default:  return false;
        }
    }

    void set_contains(uchar c, bool val)
    {
        switch(c)
        {
            case 'A': 
                map_contains[0] = val;
                break;
            case 'C':
                map_contains[1] = val;
                break;
            case 'G': 
                map_contains[2] = val;
                break;
            case 'T': 
                map_contains[3] = val;
                break;
            default:
                return;
        }
    }

public:

    ACGT_map() {
        map_contains = bit_vector(ACGT_SIZE, false);
    }

    ACGT_map(std::unordered_map<uchar, T> map) 
    {
        map_contains = bit_vector(ACGT_SIZE, false);

        for(auto const& [c, val] : map)
        {
            if(allowed(c))
            {
                set_contains(c, true);
                find(c) = val;
            }
        }
    }

    bool contains(uchar c)
    {
        switch(c)
        {
            case 'A': return map_contains[0];
            case 'C': return map_contains[1];
            case 'G': return map_contains[2];
            case 'T': return map_contains[3];
            default:  return false;
        }
    }

    bool insert(std::pair<uchar, T> kv)
    {
        if (allowed(kv.first) && !contains(kv.first)) {
            set_contains(kv.first, true);
            find(kv.first) = kv.second;

            return true;
        }
        else {
            return false;
        }
    }

    T& at(const uchar c)
    {
        if (!contains(c))
        {
            throw std::out_of_range("Symbol" + std::to_string(c) + "not in map");
        }

        return find(c);
    }

    T& operator[](const uchar c) 
    {
        return find(c);
    }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += map_contains.serialize(out, v, "map_contains");
        written_bytes += a_type.serialize(out, v, "a_type");
        written_bytes += c_type.serialize(out, v, "c_type");
        written_bytes += g_type.serialize(out, v, "g_type");
        written_bytes += t_type.serialize(out, v, "t_type");

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        map_contains.load(in);
        a_type.load(in);
        c_type.load(in);
        g_type.load(in);
        t_type.load(in);
    }
};

#endif /* end of include guard: _ACGT_MAP_HH */
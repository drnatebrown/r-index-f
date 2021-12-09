/* interval_pos - Pair describing run/offset access of r-index-f
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
   \file interval_pos.hpp
   \brief interval_pos.hpp Pair describing run/offset access of r-index-f
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 11/19/2021
*/

#ifndef _INTERVAL_POS_F_HH
#define _INTERVAL_POS_F_HH

#include <common.hpp>

#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

class interval_pos
{
private:
    ulint run;
    ulint offset;
    bool set;

public:
    interval_pos() {
        run = 0;
        offset = 0;
        set = false;
    }

    interval_pos(ulint r, ulint o) {
        run = r;
        offset = o;
    }

    bool is_set()
    {
        return set;
    }

    interval_pos& operator++() 
    {
        ++offset;
    }

    interval_pos operator++(int) 
    {
        interval_pos old = *this;
        operator++();
        return old;
    }

    inline bool operator< (const interval_pos& pos){ return (run == pos.run) ? (offset < pos.offset) : (run < pos.run); }
    inline bool operator> (const interval_pos& pos){ return (run == pos.run) ? (offset > pos.offset) : (run > pos.run); }
    inline bool operator<=(const interval_pos& pos){ return !(*this > pos); }
    inline bool operator>=(const interval_pos& pos){ return !(*this < pos); }
    inline bool operator==(const interval_pos& pos){ return run == pos.run && offset == pos.offset; }
    inline bool operator!=(const interval_pos& pos){ return !(*this == pos); }

    /* serialize the structure to the ostream
    * \param out     the ostream
    */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        out.write((char *)&run, sizeof(run));
        written_bytes += sizeof(run);

        out.write((char *)&offset, sizeof(offset));
        written_bytes += sizeof(offset);

        out.write((char *)&set, sizeof(set));
        written_bytes += sizeof(set);

        return written_bytes;
    }

    /* load the structure from the istream
    * \param in the istream
    */
    void load(std::istream &in)
    {
        in.read((char *)&run, sizeof(run));
        in.read((char *)&offset, sizeof(offset));
        in.read((char *)&set, sizeof(set));
    }
};

#endif /* end of include guard: _INTERVAL_POS_HH */
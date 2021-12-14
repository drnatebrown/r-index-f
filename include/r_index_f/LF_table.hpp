/* Lf_table - Uncompressed version of OptBWTR (LF table)
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
   \file LF_table.hpp
   \brief LF_table.hpp Uncompressed version of OptBWTR (LF table)
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 19/11/2021
*/

#ifndef _LF_TABLE_HH
#define _LF_TABLE_HH

#include <common.hpp>

using namespace std;

class LF_table
{
private:
    ulint n; // Length of BWT
    ulint r; // Runs of BWT

    vector<char> chars;
    vector<ulint> intervals;
    vector<ulint> lens;
    vector<ulint> offsets;

public:
    // Row of the LF table
    struct LF_row
    {
        char character;
        ulint interval;
        ulint length;
        ulint offset;
    };

    // TODO: Add builder for BWT (not heads/lengths)
    LF_table(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        chars = vector<char>(); 
        lens = vector<ulint>();
        
        char c;
        ulint i = 0;
        n = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
            {
                chars.push_back(c);
                lens.push_back(length);
                L_block_indices[c].push_back(i);
            }
            else
            {
                chars.push_back(TERMINATOR);
                lens.push_back(length);
                L_block_indices[TERMINATOR].push_back(i);
            }
            ++i;
            n+=length;
        }
        
        r = chars.size();

        intervals = vector<ulint>(r);
        offsets = vector<ulint>(r);

        ulint curr_L_num = 0;
        ulint L_seen = 0;
        ulint F_seen = 0;
        for(size_t i = 0; i < L_block_indices.size(); ++i) 
        {
            for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
            {
                ulint pos = L_block_indices[i][j];

                intervals[pos] = curr_L_num;
                offsets[pos] = F_seen - L_seen;

                F_seen += lens[pos];
            
                while (curr_L_num < r && F_seen >= L_seen + lens[curr_L_num]) 
                {
                    L_seen += lens[curr_L_num];
                    ++curr_L_num;
                }
            }
        }
    }

    const LF_row get(size_t i)
    {
        assert(i < r);
        return LF_row{chars[i], intervals[i], lens[i], offsets[i]};
    }

    ulint size()
    {
        return n;
    }

    ulint runs()
    {
        return r;
    }

    /*
     * \param Run position (RLE intervals)
     * \param Current character offset in block
     * \return block position and offset of preceding character
     */
    std::pair<ulint, ulint> LF(ulint run, ulint offset)
    {
        ulint next_interval = intervals[run];
	    ulint next_offset = offsets[run] + offset;

	    while (next_offset >= lens[next_interval]) 
        {
            next_offset -= lens[next_interval++];
        }

	    return std::make_pair(next_interval, next_offset);
    }
};

#endif /* end of include guard: _LF_TABLE_HH */
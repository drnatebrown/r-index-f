/* build_rif - Build the simple R-Index-F tablke
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
   \file build_LF_table.cpp
   \brief build_LF_table.cpp Build the simple R-Index-F table.
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 02/11/2021
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>
#include <sdsl/io.hpp>
#include <r_index_f.hpp>
#include <malloc_count.h>


int main(int argc, char *const argv[])
{
    Args args;
    parseArgs(argc, argv, args);

    if (args.pattern == "") {
        error("-p flag is required for count query");
        return 1;
    }

    verbose("Loading the R-Index-F from B-Table");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    r_index_f<> rif;
    std::string filename_rif = args.filename + rif.get_file_extension();

    ifstream fs_rif(filename_rif);
    rif.load(fs_rif);
    fs_rif.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("R-Index-F load complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    rif.mem_stats();
    rif.bwt_stats();

    verbose("Computing count query");
    t_insert_start = std::chrono::high_resolution_clock::now();

    if (args.is_fasta) {
        cout << "\tCOUNT: " << rif.count(args.pattern) << endl;
    }
    else {
        std::string pattern_file = args.pattern;
        ifstream fs_pattern(pattern_file);
        std::string pattern;
        size_t count = 0;
        while (std::getline(fs_pattern, pattern)) {
            cout << "P_LINE: " << count << "\tCOUNT: " << rif.count(pattern) << endl;
            ++count;
        }
    }
    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Count query complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    return 0;
}
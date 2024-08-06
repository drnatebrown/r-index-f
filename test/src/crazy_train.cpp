/* rif_tests Performs benchmarks on the constructed R-Index-F --- must serialize both r-index-f and LF_table
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
   \file rif_tests.cpp
   \brief rif_tests Benchmark tests on the R-Index-F
   \author Nathaniel Brown
   \author Massimiliano Rossi
   \date 02/11/2021
*/


#define VERBOSE

#include "LF_table.hpp"
#include <iostream>
#include <fstream> 
#include <common.hpp>
#include <r_index_f.hpp>
#include <malloc_count.h>

int main(int argc, char *const argv[])
{
    Args args;
    parseArgs(argc, argv, args);

    verbose("Loading the R-Index-F from B-Table");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    LF_table rif;
    std::string filename_rif = args.filename + rif.get_file_extension();
    if(args.d) {
        filename_rif += "." + std::to_string(args.d) + "_col";
    }

    ifstream fs_rif(filename_rif);
    rif.load(fs_rif);
    fs_rif.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("R-Index-F load complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    t_insert_start = std::chrono::high_resolution_clock::now();

    rif.get_order();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("MEOV reordering query complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    return 0;
}
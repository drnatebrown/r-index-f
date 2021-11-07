/* rif_tests Performs benchmarks on the constructed R-Index-F
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

#include <iostream>
#include <fstream> 

#define VERBOSE

#include <common.hpp>
#include <r_index_f.hpp>
#include <malloc_count.h>

void invert_bwt(r_index_f<> rif) 
{
    verbose("Inverting BWT using R-Index-F (B table)");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
    ulint steps = 0;
    ulint run = 0;
    ulint offset = 0;
    char c;
    while((c = rif.get_char(run)) > TERMINATOR) 
    {
        std::pair<ulint, ulint> block_pair = rif.LF(run, offset, c);
        run = block_pair.first;
        offset = block_pair.second;

        ++steps;
    }
    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
    verbose("BWT Inverted using B Table");
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    verbose("Average step (ns): ", std::chrono::duration<double, std::ratio<1, 1000000000>>((t_insert_end - t_insert_start)/steps).count());
}

int main(int argc, char *const argv[])
{
    Args args;
    parseArgs(argc, argv, args);

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
    invert_bwt(rif);

    return 0;
}
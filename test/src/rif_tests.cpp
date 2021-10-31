/* rif_tests Performs benchmarks on the constructed R-Index-F
    Copyright (C) 2020 Massimiliano Rossi
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
   \author Massimiliano Rossi
   \author Nathaniel Brown
   \date 03/03/2021
*/

#include <iostream>
#include <fstream> 

#define VERBOSE
#define SAMPLES 100000000
#define SEED 23

#include <common.hpp>

#include <r_index_f.hpp>

#include <malloc_count.h>


int main(int argc, char *const argv[])
{
    Args args;
    parseArgs(argc, argv, args);

    verbose("Loading the R-Index-F from LF-Table");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    r_index_f rif;

    std::string filename_rif = args.filename + rif.get_file_extension();

    ifstream fs_rif(filename_rif);
    rif.load(fs_rif);
    fs_rif.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("R-Index-F load complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    rif.mem_stats();

    rif.invert_bwt(args.filename);
    //rif.sample_LF(SAMPLES, SEED);
    //rif.print_table();

    return 0;
}

/* Can move test only methods here, as extensions of base class
class r_index_f_test : public ToBeTested, public testing::Test
{
   // Empty - bridge to protected members for unit-testing
}
*/
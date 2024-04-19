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

#include <iostream>
#include <fstream> 

#define VERBOSE

#include <common.hpp>
#include <r_index_f.hpp>
#include <malloc_count.h>

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

    t_insert_start = std::chrono::high_resolution_clock::now();

    FILE *fd;
    std::string tunnel_file = args.filename + ".tnl";
    if ((fd = fopen(tunnel_file.c_str(), "r")) == nullptr)
        error("open() file " + tunnel_file + " failed");

    vector<ulint> tunnels = vector<ulint>(rif.runs());
    for (size_t i = 0; i < rif.runs(); ++i) {
        if ((fread(&tunnels[i], 5, 1, fd)) != 1)
            error("fread() file " + tunnel_file + " failed");
    }

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Tunnels load complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    std::ifstream t("/home/nbrown99/vast/data/chr19.1.trim.raw");
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::string pattern = buffer.str();

    t_insert_start = std::chrono::high_resolution_clock::now();

    std::ofstream outfile(args.filename + ".tnl_stats");
    rif.count(pattern, tunnels, outfile);
    outfile.close();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Tunnel pattern query complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    return 0;
}
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
#include <columnar_table.hpp>
#include <alt_block.hpp>
#include <alt_table.hpp>
#include <hybrid_table.hpp>
#include <LF_table.hpp>
#include <malloc_count.h>


int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  verbose("Building the Alt Table");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  r_index_f<alt_table> alt(args.filename);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  // lf_table.mem_stats();
  // lf_table.bwt_stats();

  verbose("Construction Complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Serializing Table");

  std::string outfile = args.filename + alt.get_file_extension();
  std::ofstream out(outfile);
  alt.serialize(out);
  out.close();
  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // Building the r-index-f table
  // verbose("Building the R-Index-F");
  // t_insert_start = std::chrono::high_resolution_clock::now();

  // r_index_f<> rif(args.filename);

  // t_insert_end = std::chrono::high_resolution_clock::now();

  // verbose("Construction Complete");
  // verbose("Memory peak: ", malloc_count_peak());
  // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  //  verbose("Serializing Table");

  // outfile = args.filename + rif.get_file_extension();
  // std::ofstream out_1(outfile);
  // rif.serialize(out_1);
  // out_1.close();
  // t_insert_end = std::chrono::high_resolution_clock::now();

  // verbose("Memory peak: ", malloc_count_peak());
  // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // verbose("Building the Columnar Table");
  // t_insert_start = std::chrono::high_resolution_clock::now();

  // r_index_f<columnar_table> cc(args.filename);

  // t_insert_end = std::chrono::high_resolution_clock::now();

  // cc.mem_stats();
  // cc.bwt_stats();

  // verbose("Construction Complete");
  // verbose("Memory peak: ", malloc_count_peak());
  // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // verbose("Serializing Table");

  // outfile = args.filename + cc.get_file_extension();
  // std::ofstream out_2(outfile);
  // cc.serialize(out_2);
  // out_2.close();
  // t_insert_end = std::chrono::high_resolution_clock::now();

  // verbose("Memory peak: ", malloc_count_peak());
  // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Building the Alt-Block Table");
  t_insert_start = std::chrono::high_resolution_clock::now();

  r_index_f<alt_block<>> alt_b(args.filename);

  t_insert_end = std::chrono::high_resolution_clock::now();

  // cc.mem_stats();
  // cc.bwt_stats();

  verbose("Construction Complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Serializing Table");

  outfile = args.filename + alt_b.get_file_extension();
  std::ofstream out_3(outfile);
  alt_b.serialize(out_3);
  out_3.close();
  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  return 0;
}
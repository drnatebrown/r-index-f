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
#include <alt_table.hpp>
#include <LF_table.hpp>
#include <malloc_count.h>


int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  std::string bwt_heads_fname = args.filename + ".bwt.heads";
  std::ifstream ifs_heads(bwt_heads_fname);
  std::string bwt_len_fname = args.filename + ".bwt.len";
  std::ifstream ifs_len(bwt_len_fname);

  ifs_heads.seekg(0);
  ifs_len.seekg(0);

  verbose("Building the Alt-LF Table");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  alt_table lf_table(ifs_heads, ifs_len);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  // lf_table.mem_stats();
  // lf_table.bwt_stats();

  verbose("Construction Complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // Building the r-index-f table

  // verbose("Building the R-Index-F");
  // t_insert_start = std::chrono::high_resolution_clock::now();

  // r_index_f<> rif(lf_table);

  // t_insert_end = std::chrono::high_resolution_clock::now();

  // verbose("Construction Complete");
  // verbose("Memory peak: ", malloc_count_peak());
  // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // ifs_heads.seekg(0);
  // ifs_len.seekg(0);

  verbose("Building the Columnar Table");
  t_insert_start = std::chrono::high_resolution_clock::now();

  columnar_table cc(ifs_heads, ifs_len);

  t_insert_end = std::chrono::high_resolution_clock::now();

  // cc.mem_stats();
  // cc.bwt_stats();

  verbose("Construction Complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // verbose("Serializing Table");

  // std::string outfile = args.filename + rif.get_file_extension();
  // std::ofstream out(outfile);
  // rif.serialize(out);
  // out.close();
  // t_insert_end = std::chrono::high_resolution_clock::now();

  // verbose("Memory peak: ", malloc_count_peak());
  // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  return 0;
}
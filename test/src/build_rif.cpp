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

#define VERBOSE

#include "LF_table.hpp"
#include <common.hpp>
#include <sdsl/io.hpp>
#include <r_index_f.hpp>
#include <malloc_count.h>


int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index-f table

  verbose("Building the R-Index-F");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  std::string bwt_fname = args.filename + ".bwt";
  std::string bwt_heads_fname = bwt_fname + ".heads";
  std::ifstream ifs_heads(bwt_heads_fname);
  std::string bwt_len_fname = bwt_fname + ".len";
  std::ifstream ifs_len(bwt_len_fname);

  LF_table rif;

  if (args.d) {
    std::string splitting_filename = args.filename + "." + std::to_string(args.d) + "_col";
    cout << splitting_filename << endl;
    std::ifstream ifs_split(splitting_filename);
    bit_vector run_splits;
    run_splits.load(ifs_split);

    rif = LF_table(ifs_heads, ifs_len, run_splits);
  }
  else {
    rif = LF_table(ifs_heads, ifs_len);
  }

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  rif.mem_stats();
  rif.bwt_stats();

  verbose("Construction Complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Serializing Table");

  std::string outfile = args.filename + rif.get_file_extension();
  if (args.d) {
    outfile += "." + std::to_string(args.d) + "_col";
  }
  std::ofstream out(outfile);
  rif.serialize(out);
  out.close();
  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  return 0;
}
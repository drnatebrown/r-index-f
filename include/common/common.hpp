/*  common
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
   \file common.hpp
   \brief common.hpp contains common features.
   \author Massimiliano Rossi
   \date 12/03/2020
*/

#ifndef _COMMON_HH
#define _COMMON_HH

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <sys/time.h>

#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <sstream>      // std::stringstream
#include <vector>      // std::vector
#include <chrono>       // high_resolution_clock
#include <algorithm> 
#include <utility> 
#include <assert.h>

#include <sdsl/int_vector.hpp>

#define ALPHABET_SIZE 256

static const uint8_t TERMINATOR = 1;
typedef unsigned long int ulint;
typedef unsigned char uchar;

std::string NowTime();
void _internal_messageInfo(const std::string message);
void _internal_messageWarning( const std::string file, const unsigned int line, const std::string message);
void _internal_messageError( const std::string file, const unsigned int line,const std::string message);


std::string NowTime()
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    char buffer[100];
    tm r;
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&tv.tv_sec, &r));
    char result[100];
    snprintf(result, 100, "%s"/*.%06ld"*/, buffer/*, (long)tv.tv_usec*/);
    return result;
}


template<typename T>
inline void _internal_message_helper(std::stringstream &ss, T const &first) { ss << first; }
template<typename T, typename... Args>
inline void _internal_message_helper(std::stringstream &ss, T const &first, const Args&... args) { ss << first << " "; _internal_message_helper(ss,args...); }
template<typename T, typename... Args>
inline std::string _internal_message(T const &first, const Args&... args) { std::stringstream ss; _internal_message_helper(ss,first,args...); return ss.str(); }


void _internal_messageInfo(const std::string message)
{
  std::cout << "[INFO] " << NowTime() << " - " << "Message: " << message << std::endl;
}

void _internal_messageWarning( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cout << "[WARNING] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
}

void _internal_messageError( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cerr << "[ERROR] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
  assert( false );
  exit( 1 );
}



#define info( args... ) \
    _internal_messageInfo( _internal_message(args) )

#ifdef VERBOSE
  #define verbose( args... ) \
      _internal_messageInfo( _internal_message(args) )
#else
  #define verbose( args... )
#endif

#define warning( args... ) \
    _internal_messageWarning( __FILE__, __LINE__, _internal_message(args) )

#define error( args... ) \
    _internal_messageError( __FILE__, __LINE__, _internal_message(args) )

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  bool store = false; // store the data structure in the file
  bool memo  = false; // print the memory usage
  bool rle   = false; // outpt RLBWT
  size_t th = 1; // number of threads
  bool is_fasta = false; // read a fasta file
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string sarg;
  while ((c = getopt(argc, argv, "w:smcfl:rhp:t:")) != -1)
  {
    switch (c)
    {
    case 's':
      arg.store = true;
      break;
    case 'm':
      arg.memo = true;
      break;
    case 'r':
      arg.rle = true;
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'f':
      arg.is_fasta = true;
      break;
    case '?':
      error("Unknown option.\n");
      break;
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 1)
  {
    arg.filename.assign(argv[optind]);
  }
  else
  {
    error("Invalid number of arguments\n");
  }
}

//********** end argument options ********************

// Convert boolean vector to specified bit vector
template<class B>
B bool_to_bit_vec(std::vector<bool> &b)
{
  if(b.size()==0) return B();

  sdsl::bit_vector bv(b.size());

  for(size_t i = 0; i < b.size(); ++i)
    bv[i] = b[i];

  return B(bv);
}

#endif /* end of include guard: _COMMON_HH */
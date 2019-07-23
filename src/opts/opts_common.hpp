/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <iostream>
// command-line option handling
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

// some globals for option handling
extern int ac; 
extern char** av; // TODO: write it down to a file as in icicle ... write the default (i.e. not specified) values as well!
extern po::options_description opts_main; 

void handle_opts(
  po::options_description &opts_micro,
  po::variables_map &vm 
);


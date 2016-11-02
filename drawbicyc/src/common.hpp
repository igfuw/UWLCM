#pragma once

#include <iostream>
#include <iomanip>
#include <set>
#include <string>
#include <sstream>
#include <vector>
#include "gnuplot.hpp"
#include <boost/units/systems/si.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

// error reporting
#define error_macro(msg) \
{ \
  cerr << "error: " << msg << endl; \
  throw exception(); \
}

#define notice_macro(msg) \
{ \
  cerr << " info: " << msg << endl; \
}

using std::set;
using std::string;
using std::ostringstream;
using std::cerr;
using std::endl;
using std::exception;
using std::vector;
namespace si = boost::units::si;
using boost::units::quantity;
namespace po = boost::program_options;

// some globals for option handling
int ac;
char** av;
po::options_description opts_main("General options");

// command-line option handling
void handle_opts(
  po::options_description &opts_new,
  po::variables_map &vm
)
{
  opts_main.add(opts_new);
//  po::store(po::parse_command_line(ac, av, opts_main), vm); // could be exchanged with a config file     parser
  po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); //         ignores unknown


  // handling the "help" option
  if (vm.count("help"))
  {
    std::cout << opts_main;
    exit(EXIT_SUCCESS);
  }
  po::notify(vm); // includes checks for required options
}

string zeropad(int n, int w=3)
{
  std::ostringstream tmp;
  tmp << std::setw(w) << std::setfill('0') << n;
  return tmp.str();
}

using namespace blitz;

double iscloudy_rc(double x)
{
  return x > 1.e-5 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(iscloudy_rc)

double is_th_prtrb(double x)
{
  return x > 300.1 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(is_th_prtrb)

double iscloudy(double x)
{
  return x > 20. ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(iscloudy)

double isdowndraught(double x)
{
  return  x < -0.2 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(isdowndraught)

double isupdraught(double x)
{
  return  x > 0.2 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(isupdraught)



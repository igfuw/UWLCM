/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// command-line option handling
#include "opts_common.hpp"

// define globals for option handling
int ac; 
char** av; // TODO: write it down to a file as in icicle ... write the default (i.e. not specified) values as well!
po::options_description opts_main("General options"); 

void handle_opts(
  po::options_description &opts_micro,
  po::variables_map &vm,
  const bool check_help
)
{
  opts_main.add(opts_micro);
    po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); // ignores unknown, could be exchanged with a config file parser

  // hendling the "help" option
  if (check_help && vm.count("help"))
  {
    std::cout << opts_main;
    exit(EXIT_SUCCESS);
  }
  po::notify(vm); // includes checks for required options
}


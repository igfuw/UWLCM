/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <set>

// signal handling (kill, Ctrl+c)
#if defined(__linux__)
#  include <signal.h>
#endif

namespace
{
  bool *panic;
  
  void panic_handler(int)
  {
    *panic = true;
  }
  
  void set_sigaction()
  {
  #if defined(__linux__)
    //const struct sigaction sa({.sa_handler = panic_handler}); // gcc fails to compile it (TODO: report it)
    struct sigaction sa = {0};
    sa.sa_handler = panic_handler;
    for (auto &s : std::set<int>({SIGTERM, SIGINT})) sigaction(s, &sa, NULL);
  #endif
  }
};

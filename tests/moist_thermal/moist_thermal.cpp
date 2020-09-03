// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <cstdlib> // system()
#include <sstream> // std::ostringstream

#include <drawbicyc/notice_macros.hpp>
#include "common.hpp"

using std::ostringstream;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  for (auto &opts_m : opts_micro)
  {
    // run the simulation
    ostringstream cmd;
    cmd << av[1] <<  "/../../build/uwlcm " << opts_common << " " << opts_m.second;
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}

#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string outdir;
  string opts_common = 
    "--outfreq=1000 --nt=2 --spinup=1 --nx=4 --nz=4 --dt=1";
  set<string> opts_micro({
    "--adv_serial=true --async=false --micro=lgrngn --outdir=out_lgrngn --backend=OpenMP --sd_conc=8 z_rlx_sclr=100"  
  });

  for (auto &opts_m : opts_micro)
  {
    ostringstream cmd;
    cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m;  
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}

#include <cstdlib> // system()
#include <unordered_set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"

using std::ostringstream;
using std::unordered_set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string outdir;
  string opts_common = 
    "--outfreq=1000 --nt=2 --spinup=1 --dt=1 --adv_serial=true";
  unordered_set<string> opts_dim({
    "--nx=4 --nz=4",
    "--nx=4 --ny=4 --nz=4"
  });
  unordered_set<string> opts_micro({
    "--async=false --micro=lgrngn --outdir=out_lgrngn --backend=serial --sd_conc=8 --z_rlx_sclr=100 --unit_test=true",
    "--micro=blk_1m --outdir=out_blk_1m"  
  });

  for (auto &opts_d : opts_dim)
    for (auto &opts_m : opts_micro)
    {
      ostringstream cmd;
      cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m << " " << opts_d;
      notice_macro("about to call: " << cmd.str())

      if (EXIT_SUCCESS != system(cmd.str().c_str()))
        error_macro("model run failed: " << cmd.str())
    }
}

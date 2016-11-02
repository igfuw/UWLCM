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
    "--outfreq=60 --nt=600 --spinup=1000 --dt=1 --nx=181 --nz=121 --th_src=false --uv_src=false --rv_src=false --w_src=true --adv_serial=true";
  unordered_set<string> opts_micro({
//    "--async=false --micro=lgrngn --outdir=out_lgrngn --backend=serial --sd_conc=8 --z_rlx_sclr=100 --unit_test=true",
    "--micro=blk_1m --outdir=tmp --cond=true --cevp=true --revp=false --conv=false --accr=false --sedi=false"  
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

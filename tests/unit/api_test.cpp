#include <cstdlib> // system()
#include <vector>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"

using std::ostringstream;
using std::vector;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string outdir;
  string opts_common = 
    "--outfreq=1000 --nt=2 --spinup=1 --dt=0.1 --serial=true"; // dt=1 caused blk1m dycoms to freeze on pressure solver
  vector<string> opts_dim({
    "--nx=4 --ny=4 --nz=4",
    "--nx=4 --nz=4"
  });
  vector<string> opts_micro({
    "--async=false --micro=lgrngn --outdir=out_lgrngn --backend=serial --sd_conc=8 --z_rlx_sclr=100 --unit_test=0",
    "--micro=blk_1m --outdir=out_blk_1m"  
  });
  vector<string> opts_case({
    "--case=dry_thermal --cond=0 --coal=0",
    "--case=moist_thermal",
    "--case=dycoms"
  });

  for (auto &opts_d : opts_dim)
    for (auto &opts_m : opts_micro)
      for (auto &opts_c : opts_case)
      {
        if((opts_c ==  opts_case[0] || opts_c == opts_case[1]) && opts_d == opts_dim[0])
        {
          std::cout << "skipping 3d thermal tests" << std::endl;
          continue; 
        }
        ostringstream cmd;
        cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m << " " << opts_d << " " << opts_c;
        notice_macro("about to call: " << cmd.str())

        if (EXIT_SUCCESS != system(cmd.str().c_str()))
          error_macro("model run failed: " << cmd.str())
      }
}

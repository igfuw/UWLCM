#include <cstdlib> // system()
#include <vector>
#include <string>
#include <sstream> // std::ostringstream
#include <fstream>
#include <unordered_map>

#include "../common.hpp"

using std::ostringstream;
using std::vector;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2 && ac != 3) error_macro("expecting one or two arguments: 1. CMAKE_BINARY_DIR 2. additional command line options (optional)");
  string opts_additional = ac == 3 ? av[2] : "";

  string opts_common = 
    "--outfreq=1000 --nt=2 --spinup=1 --dt=1 --serial=true --prs_tol=1e-3 --rng_seed=44"; 
  vector<string> opts_dim({
    "--nx=4 --nz=4",
    "--nx=4 --ny=4 --nz=4"
  });
  string opts_dim_lgrngn_3d =  "--nx=40 --ny=40 --nz=40"; // = 4 caused multiplicity overflows in the Lagrangian 3D
  vector<string> opts_micro({
    "--micro=blk_1m"  ,
    "--async=false --micro=lgrngn --backend=serial --sd_conc=8 --z_rlx_sclr=100"
  });
  vector<string> opts_case({
    "--case=moist_thermal",
    "--case=dry_thermal --cond=0 --coal=0",
//    "--case=dycoms_rf01",
    "--case=dycoms_rf02",
    "--case=dycoms_rf02 --gccn=1 --out_dry_spec=1 --out_wet_spec=1",
//    "--case=lasher_trapp"
  });
  vector<string> opts_piggy({
    "--piggy=0",
    "--piggy=0 --save_vel=1",
    "--piggy=1 --vel_in=velocity_out.dat"  // take vel file from blk, cause it's ran first
  });
  string opts_pig_from_lgrngn = "--piggy=1 --vel_in=velocity_out.dat";

  system("mkdir output");

  // file with a dict translating hashed outdir into options
  std::ofstream ofdict("hash_dict"+opts_additional+".txt");

  for (auto &opts_d : opts_dim)
    for (auto &opts_m : opts_micro)
      for (auto &opts_c : opts_case)
        for (auto &opts_p : opts_piggy) // piggy has to be last to prevent overwriting of vel_out
        {
          if((opts_c == opts_case[1]) && opts_d == opts_dim[1])
          {
            std::cout << "skipping 3d dry thermal tests" << std::endl;
            continue; 
          }

          // more cells in 3d lgrngn to avoid n overflow
          if(opts_d == opts_dim[1] && opts_m == opts_micro[1])
          {
            opts_d = opts_dim_lgrngn_3d;
            if(opts_p == opts_piggy[2])
              opts_p = opts_pig_from_lgrngn;
          }
          ostringstream cmd, opts;
          opts << opts_common << " " << opts_m << " " << opts_d << " " << opts_c << " " << opts_p << " " << opts_additional;
          // we want outdir=opts.str(), but that gives a long outdir that h5diff has trouble with reading and gives error when comparing const.h5 with refdata
          // hence we hash opts to get outdir
          auto outdir = std::hash<std::string>{}(opts.str());
          ofdict << outdir << " : " << opts.str() << std::endl;

          cmd << av[1] << "/src/bicycles " << opts.str() << " --outdir=\"output/" << outdir << "\"";
          notice_macro("about to call: " << cmd.str())
  
          if (EXIT_SUCCESS != system(cmd.str().c_str()))
            error_macro("model run failed: " << cmd.str())

          // copy the stored velocity for the next run
          if(opts_p == opts_piggy[1])
          {
            ostringstream cpcmd;
            cpcmd << "cp \"output/" << outdir << "/velocity_out.dat\" .";
            notice_macro("about to call: " << cpcmd.str())
            system(cpcmd.str().c_str());
          }
        }
}

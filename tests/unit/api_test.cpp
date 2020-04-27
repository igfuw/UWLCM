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
  vector<string> opts_micro({
    "--micro=blk_1m"  ,
    "--async=false --micro=lgrngn --backend=serial --sd_conc=8" 
  });
  vector<string> opts_case({
    "--case=moist_thermal_api_test",
    "--case=dry_thermal_api_test --cond=0 --coal=0",
    "--case=dycoms_rf02_api_test",
    "--case=dycoms_rf02_api_test --gccn=1 --out_dry_spec=1 --out_wet_spec=1",
    "--case=rico11_api_test",
    "--case=dycoms_rf01_api_test",
    "--case=lasher_trapp_api_test"
  });
  vector<string> opts_piggy({
    "--piggy=0",
    "--piggy=0 --save_vel=1",
    "--piggy=1 --vel_in=velocity_out.dat"  // take vel file from blk, cause it's ran first
  });

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
          if((opts_c == opts_case[1]) && opts_m == opts_micro[1])
          {
            std::cout << "skipping dry thermal tests with Lagrangian microphysics" << std::endl;
            continue; 
          }

          ostringstream cmd, opts;
          opts << opts_common << " " << opts_m << " " << opts_d << " " << opts_c << " " << opts_p << " " << opts_additional;
          // we want outdir=opts.str(), but that gives a long outdir that h5diff has trouble with reading and gives error when comparing const.h5 with refdata
          // hence we hash opts to get outdir
          auto outdir = std::hash<std::string>{}(opts.str());
          ofdict << outdir << " : " << opts.str() << std::endl;

          cmd << av[1] <<  "/../../build/uwlcm " << opts.str() << " --outdir=\"output/" << outdir << "\"";
 
          cerr << endl << "=========" << endl;
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

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
  if (ac != 3 && ac != 4) error_macro("expecting two or three arguments: 1. CMAKE_BINARY_DIR 2. should piggybacking be tested (bool) 3. additional command line options (optional)");
  string opts_additional = ac == 4 ? av[3] : "";
  bool run_piggy = std::stoi(av[2]);

  string opts_common = 
    "--outfreq=1000 --nt=2 --spinup=1 --dt=1 --serial=true --prs_tol=1e-3 --rng_seed=44 --case_n_stp_multiplier=1e-8"; 
  vector<string> opts_dim({
    "--nx=4 --nz=4",
    "--nx=4 --nz=4 --X=1000 --Z=-1",
    "--nx=4 --ny=4 --nz=4",
    "--nx=4 --ny=4 --nz=4 --X=1000 --Y=1000 --Z=-1"
  });
  vector<string> opts_micro({
    "--micro=blk_1m"  ,
    "--micro=blk_2m"  ,
    "--async=false --micro=lgrngn --backend=serial --sd_conc=8",
    "--async=false --micro=lgrngn --backend=serial --sd_conc=8 --gccn=1",
    "--async=false --micro=lgrngn --backend=serial --sd_conc=8 --relax_ccn=1",
    "--async=false --micro=lgrngn --backend=serial --sd_conc=8 --gccn=1 --relax_ccn=1"
  });
  vector<string> opts_case({
    "--case=moist_thermal",
    "--case=dry_thermal --cond=0 --coal=0",
    "--case=dycoms_rf02",
    "--case=dycoms_rf02 --out_dry_spec=1 --out_wet_spec=1",
    "--case=dycoms_rf02 --relax_th_rv=1",
    "--case=rico11",
    "--case=dycoms_rf01",
    "--case=cumulus_congestus"
  });
  vector<string> opts_piggy({
    "--piggy=0",
    "--piggy=0 --save_vel=1"
  });
  // run the piggybacker, if required
  if(run_piggy)
    opts_piggy.push_back("--piggy=1 --vel_in=velocity_out.dat");  // take vel file from blk, cause it's ran first

  system("mkdir output");

  // file with a dict translating hashed outdir into options
  std::ofstream ofdict("hash_dict"+opts_additional+".txt");

  for (auto &opts_d : opts_dim)
    for (auto &opts_m : opts_micro)
      for (auto &opts_c : opts_case)
        for (auto &opts_p : opts_piggy) // piggy has to be last to prevent overwriting of vel_out
        {
          if((opts_c == opts_case[1]) && (opts_d == opts_dim[2] || opts_d == opts_dim[3]))
          {
            std::cout << "skipping 3d dry thermal tests" << std::endl;
            continue;
          }
          if((opts_c == opts_case[1]) && opts_m != opts_micro[0])
          {
            std::cout << "skipping dry thermal tests with microphysics other than blk_1m" << std::endl;
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

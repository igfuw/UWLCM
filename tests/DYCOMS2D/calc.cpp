#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"
#include "bins.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string bins_wet_str, bins_dry_str, outdir;

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_dry();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_dry_str = tmp.str();
  }

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_wet();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_wet_str = tmp.str();
  }
/*  outdir=std::getenv("STORAGE");
  string tmp=std::getenv("NODE_CONF");
  outdir+="/"+tmp+"/out_lgrngn";
  std::cout << "output directory: " << outdir << std::endl;
*/
  string opts_common = 
    "--outfreq=200 --nt=12000 --spinup=9600 --nx=97 --nz=301 --dt=0.75 --relax_th_rv=false"; // DYCOMS: 128x300 ; 600 21600 3600
  set<string> opts_micro({
//    "--micro=blk_1m --outdir=out_blk_1m",
//    "--micro=lgrngn --outdir=out_lgrngn",
    "--micro=blk_2m --outdir=out_blk_2m \
     --adv_serial=false --async=true  --backend=OpenMP --case=dycoms_rf01 --piggy=false \
     --acnv_A=1350 --acnv_b=2.47 --acnv_c=-1.79 \
     --blk2m_mean_rd=0.05e-6 --blk2m_sdev_rd=1.5 --blk2m_N_stp=080e6 \
     --w_src=1 --uv_src=1 --rv_src=1 --th_src=1 --subsidence=1 "
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

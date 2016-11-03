#include <cstdlib> // system()
#include <unordered_set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"
#include "../../drawbicyc/src/PlotterMicro.hpp"
#include "../../drawbicyc/src/common_filters.hpp"

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
    "--micro=blk_1m --outdir=tmp_out --cond=true --cevp=true --revp=false --conv=false --accr=false --sedi=false"  
  });

  for (auto &opts_m : opts_micro)
  {
    // run the simulation
    ostringstream cmd;
    cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m;
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
 
    // instantiate plotter to read in output
    using Plotter_t = PlotterMicro_t<2>;
    Plotter_t plotter("tmp_out", "blk_1m");
    auto& n = plotter.map;
    float data_com[] = {0., 0., 915.282, 1006.67, 1101.82, 1183.7, 1245.3, 1292.47, 1326.17, 1340.71, 1346.41};
    float data_avg[] = { 0, 1.3249e-05, 0.00011178, 0.000275819, 0.000421079, 0.000552953, 0.000621904, 0.000583853, 0.00051231, 0.000448865, 0.000410936};

    blitz::Array<float, 1> result_com(n["t"]);
    blitz::Array<float, 1> result_avg(n["t"]);
    blitz::Array<double, 1> rel_err(n["t"]);
    double eps = 1e-5; // relative precision; mostly due to data_ number of digits

    // read statistics
    for (int at = 0; at < n["t"]; ++at)
    {
      {
        auto tmp = plotter.h5load_rc_timestep(plotter.file, at * 60);
        typename Plotter_t::arr_t snap(tmp);

        // center of mass of cloud droplets
        {
          typename Plotter_t::arr_t snap2(tmp);
          snap2 = snap2 * plotter.LastIndex * n["dz"];
          if(blitz::sum(snap) > 1e-3)
            result_com(at) = blitz::sum(snap2) / blitz::sum(snap);
          else
            result_com(at) = 0.;
        }
        // cloud water mixing ratio averaged over cells with r_c>1e-5
        {
          typename Plotter_t::arr_t snap2(tmp);
          snap2 = iscloudy_rc(snap2); // find cells with rc>1e-5
          snap *= snap2; // apply filter

          // mean only over updraught cells
          if(blitz::sum(snap2) > 0.)
            result_avg(at) = blitz::sum(snap) / blitz::sum(snap2);
          else
            result_avg(at) = 0.;
        }
      }
    }

    // compare center of mass
    blitz::Array<float, 1> expected_result(data_com, 11, blitz::neverDeleteData);
    rel_err = where(expected_result > 0, abs(result_com - expected_result) / expected_result, 0);
    if(any(rel_err > eps))
    {
      std::cerr << rel_err;
      error_macro("cloud droplets center of mass discrepancy");
    }

    // compare avg rc
    expected_result = blitz::Array<float, 1>(data_avg, 11, blitz::neverDeleteData);
    rel_err = where(expected_result > 0, abs(result_avg - expected_result) / expected_result, 0);
    if(any(rel_err > eps))
    {
      std::cerr << rel_err;
      error_macro("average cloud water mixing ratio discrepancy");
    }
  }
}

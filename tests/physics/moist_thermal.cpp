// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <cstdlib> // system()
#include <unordered_map>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"
#include "../../drawbicyc/PlotterMicro.hpp"
#include "../../drawbicyc/common_filters.hpp"

using std::ostringstream;
using std::unordered_map;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string outdir;
  string opts_common = 
    "--outfreq=60 --nt=600 --spinup=1000 --dt=1 --nx=181 --nz=121 --case=moist_thermal";
  unordered_map<string, string> opts_micro({
    {"blk_1m", "--micro=blk_1m --outdir=tmp_out_blk_1m --cond=true --cevp=true --revp=false --conv=false --accr=false --sedi=false"},
    {"lgrngn", "--micro=lgrngn --outdir=tmp_out_lgrngn --cond=true --adve=true --sedi=false --coal=false --backend=OpenMP --sd_conc=16"}
  });

  // expected results
  // center of mass of rc
  unordered_map<string, std::array<float, 11>> data_com = {
    {"blk_1m", {{0., 0., 915.282, 1006.67, 1101.82, 1183.7, 1245.3, 1292.47, 1326.17, 1340.71, 1346.41}}},
    {"lgrngn", {{0, 0, 914.761, 1006.64, 1100.94, 1179.43, 1235.54, 1282.71, 1337.99, 1370.28, 1344.92}}}
  };
  // average rc
  unordered_map<string, std::array<float, 11>> data_avg = {
    {"blk_1m", {{0, 1.3249e-05, 0.000111779, 0.000275816, 0.000421085, 0.000552938, 0.000621531, 0.000585304, 0.000513864, 0.000440379, 0.000406745}}},
    {"lgrngn", {{0, 0, 9.43111e-05, 0.000258181, 0.00041635, 0.000533337, 0.000590126, 0.000575933, 0.000496402, 0.000391189, 0.000295669}}}
  };
  // relative precision at given timestep
  unordered_map<string, std::array<float, 11>> eps = { 
    {"blk_1m", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 2.5e-2, 5e-2}} },  // why larger near the end?
    {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 15e-2, 5e-1, 5e-1, 1}} }   // during evaporation we get large fluctuations
  };
  // out dir
  unordered_map<string, string> tmp_out = { {"blk_1m", "tmp_out_blk_1m"}, {"lgrngn", "tmp_out_lgrngn"}};

  for (auto &opts_m : opts_micro)
  {
    // run the simulation
    ostringstream cmd;
    cmd << av[1] << "/src/bicycles " << opts_common << " " << opts_m.second;
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
 
    // instantiate plotter to read in output
    using Plotter_t = PlotterMicro_t<2>;
    Plotter_t plotter(tmp_out[opts_m.first], opts_m.first);
    auto& n = plotter.map;

    blitz::Array<float, 1> result_com(n["t"]);
    blitz::Array<float, 1> result_avg(n["t"]);
    blitz::Array<double, 1> rel_err(n["t"]);

    // read statistics
    for (int at = 0; at < n["t"]; ++at)
    {
      {
        auto tmp = plotter.h5load_ract_timestep(plotter.file, at * 60);
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
    blitz::Array<float, 1> expected_result(data_com[opts_m.first].data(), 11, blitz::neverDeleteData);
    blitz::Array<float, 1> epsilon(eps[opts_m.first].data(), 11, blitz::neverDeleteData);
    rel_err = where(expected_result > 0, abs(result_com - expected_result) / expected_result - epsilon, 0);
    if(any(rel_err > 0.))
    {
      std::cerr << result_com;
      std::cerr << rel_err;
      error_macro("cloud droplets center of mass discrepancy");
    }

    // compare avg rc
    expected_result = blitz::Array<float, 1>(data_avg[opts_m.first].data(), 11, blitz::neverDeleteData);
    rel_err = where(expected_result > 0, abs(result_avg - expected_result) / expected_result - epsilon, 0);
    if(any(rel_err > 0.))
    {
      std::cerr << result_avg;
      std::cerr << rel_err;
      error_macro("average cloud water mixing ratio discrepancy");
    }
  }
}

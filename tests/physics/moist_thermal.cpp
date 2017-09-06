// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <cstdlib> // system()
#include <unordered_map>
#include <string>
#include <sstream> // std::ostringstream

//#include "../common.hpp"
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
//    {"blk_1m", {{0., 0.,      915.282, 1006.67, 1101.82, 1183.7,  1245.3,  1292.47, 1326.17, 1340.71, 1346.41}}}, // old values from before the twomey SD bubble paper
    {"blk_1m", {{0., 841.758, 917.667, 1022.48, 1134.36, 1232.02, 1308.47, 1364.44, 1412.22, 1476.16, 1521.43}}},  // new values, i.e. for env profs from the twomey SD paper and for abs instead of iga&fct
    {"lgrngn", {{0, 0, 914.761, 1006.64, 1100.94, 1179.43, 1235.54, 1282.71, 1337.99, 1370.28, 1344.92}}} // old values are still good, because of the large error allowed here?
  };
  // average rc
  unordered_map<string, std::array<float, 11>> data_avg = {
//    {"blk_1m", {{0, 1.3249e-05, 0.000111779, 0.000275816, 0.000421085, 0.000552938, 0.000621531, 0.000585304, 0.000513864, 0.000440379, 0.000406745}}}, // old values from before the twomey SD bubble paper (ammonium sulphate aerosol + old env_profs + iga&fct)
    {"blk_1m", {{0, 0, 0.000185381, 0.000357674, 0.000518615, 0.000635203, 0.000707202, 0.000733656, 0.000686596, 0.000556995, 0.000425793 }}}, // new values, i.e. for env profs from the twomey SD paper and for abs instead of iga&fct
//    {"lgrngn", {{0, 0, 9.43111e-05, 0.000258181, 0.00041635, 0.000533337, 0.000590126, 0.000575933, 0.000496402, 0.000391189, 0.000295669}}} // old values from before the twomey SD bubble paper (ammonium sulphate aerosol + old env_profs + iga&fct)
    {"lgrngn", {{ 0, 0, 0.000179877, 0.000374551, 0.000552171, 0.000704754, 0.000790675, 0.000791144, 0.00073452, 0.000648194, 0.000548742}}}  // values for NaCl and env_profs used in the twomey SD bubble paper
  };
  // relative precision at given timestep
  unordered_map<string, std::array<float, 11>> eps = { 
    {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },  // why larger near the end?
    {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 10e-2, 15e-2, 5e-1, 5e-1, 1.5}} }   // during evaporation we get large fluctuations
  };
  // out dir
  unordered_map<string, string> tmp_out = { {"blk_1m", "tmp_out_blk_1m"}, {"lgrngn", "tmp_out_lgrngn"}};
  bool err_flag = 0;

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

    blitz::Array<float, 1> result(n["t"]);
    blitz::Array<double, 1> rel_err(n["t"]);

    // compare statistics
    // height of the center of mass of cloud droplets
    for (int at = 0; at < n["t"]; ++at)
      result(at) = plotter.act_com_z_timestep(at * 60); 

    blitz::Array<float, 1> expected_result(data_com[opts_m.first].data(), 11, blitz::neverDeleteData);
    blitz::Array<float, 1> epsilon(eps[opts_m.first].data(), 11, blitz::neverDeleteData);
    rel_err = where(expected_result > 0, abs(result - expected_result) / expected_result - epsilon, 0);
    std::cout << "height of the center of mass: " << result;
    if(any(rel_err > 0.))
    {
      std::cerr << "ERROR" << std::endl;
      std::cerr << "expected result: " << expected_result;
      std::cerr << "relative error minus precision: " << rel_err;
      err_flag = 1;
    }

    // average cloud water mixing ratio in cloudy cells
    for (int at = 0; at < n["t"]; ++at)
    {
      {
        auto tmp = plotter.h5load_ract_timestep(at * 60);
        typename Plotter_t::arr_t ract(tmp);
        typename Plotter_t::arr_t mask(tmp);
        mask = iscloudy_rc(mask);
        ract *= mask; // apply filter

        if(blitz::sum(mask) > 0.)
          result(at) = blitz::sum(ract) / blitz::sum(mask);
        else
          result(at) = 0.;
      }
    }
    expected_result = blitz::Array<float, 1>(data_avg[opts_m.first].data(), 11, blitz::neverDeleteData);
    rel_err = where(expected_result > 0, abs(result - expected_result) / expected_result - epsilon, 0);
    std::cout << "average cloud water mixing ratio in cloudy cells: " << result;
    if(any(rel_err > 0.))
    {
      std::cerr << "ERROR" << std::endl;
      std::cerr << "expected result: " << expected_result;
      std::cerr << "relative error minus precision: " << rel_err;
      err_flag = 1;
    }

    // average concentration of activated droplets in cloudy cells
    for (int at = 0; at < n["t"]; ++at)
    {
      {
        auto tmp = plotter.h5load_ract_timestep(at * 60);
        typename Plotter_t::arr_t ract(tmp);
        typename Plotter_t::arr_t mask(tmp);
        mask = iscloudy_rc(mask);
        ract *= mask; // apply filter

        if(blitz::sum(mask) > 0.)
          result(at) = blitz::sum(ract) / blitz::sum(mask);
        else
          result(at) = 0.;
      }
    }
  }

  if(err_flag)
    error_macro("error in one of the statistics");    
}

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
using barr1d = blitz::Array<double, 1>;

bool errcheck(barr1d result, barr1d expected_result, barr1d epsilon)
{
  barr1d rel_err(result.shape());
  rel_err = where(expected_result > 0, abs(result - expected_result) / expected_result - epsilon, 0);
  if(any(rel_err > 0.))
  {
    std::cerr << "ERROR" << std::endl;
    std::cerr << "expected result: " << expected_result;
    std::cerr << "relative error minus precision: " << rel_err;
    return 1;
  }
  else
    return 0;
}

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string outdir;
  string opts_common = 
    "--outfreq=60 --nt=600 --spinup=1000 --dt=1 --nx=181 --nz=121 --case=moist_thermal";
  unordered_map<string, string> opts_micro({
    {"blk_1m", "--micro=blk_1m --outdir=tmp_out_blk_1m --cond=true --cevp=true --revp=false --conv=false --accr=false --sedi=false"},
    {"lgrngn", "--micro=lgrngn --outdir=tmp_out_lgrngn --cond=true --adve=true --sedi=false --coal=false --backend=OpenMP --sd_conc=32"}
  });

  // expected results
  // center of mass of rc
  unordered_map<string, std::array<double, 11>> data_com = {
//    {"blk_1m", {{0., 0.,      915.282, 1006.67, 1101.82, 1183.7,  1245.3,  1292.47, 1326.17, 1340.71, 1346.41}}}, // old values from before the twomey SD bubble paper
    {"blk_1m", {{0., 841.758, 917.667, 1022.48, 1134.36, 1232.02, 1308.47, 1364.44, 1412.22, 1476.16, 1521.43}}},  // new values, i.e. for env profs from the twomey SD paper and for abs instead of iga&fct
    {"lgrngn", {{0, 0, 914.761, 1006.64, 1100.94, 1179.43, 1235.54, 1282.71, 1337.99, 1370.28, 1344.92}}} // old values are still good, because of the large error allowed here?
  };
  // average rc
  unordered_map<string, std::array<double, 11>> data_avg = {
//    {"blk_1m", {{0, 1.3249e-02, 0.111779, 0.275816, 0.421085, 0.552938, 0.621531, 0.585304, 0.513864, 0.440379, 0.406745}}}, // old values from before the twomey SD bubble paper (ammonium sulphate aerosol + old env_profs + iga&fct)
    {"blk_1m", {{0, 0, 0.185381, 0.357674, 0.518615, 0.635203, 0.707202, 0.733656, 0.686596, 0.556995, 0.425793 }}}, // new values, i.e. for env profs from the twomey SD paper and for abs instead of iga&fct
//    {"lgrngn", {{0, 0, 9.43111e-02, 0.258181, 0.41635, 0.533337, 0.590126, 0.575933, 0.496402, 0.391189, 0.295669}}} // old values from before the twomey SD bubble paper (ammonium sulphate aerosol + old env_profs + iga&fct)
    {"lgrngn", {{ 0, 0, 0.179877, 0.374551, 0.552171, 0.704754, 0.790675, 0.791144, 0.73452, 0.648194, 0.548742}}}  // values for NaCl and env_profs used in the twomey SD bubble paper
  };
  // std_dev_rc rc
  unordered_map<string, std::array<double, 11>> data_std_dev_rc = {
    {"blk_1m", {{0, 0, 0.0189641, 0.0669555, 0.137429, 0.210689, 0.266551, 0.30766, 0.333871, 0.310118, 0.232963}}},
    {"lgrngn", {{0, 0, 0.0551891, 0.116698, 0.202856, 0.290047, 0.375072, 0.442025, 0.485523, 0.525853, 0.523559}}} 
  };
  // relative precision at given timestep
  unordered_map<string, std::array<double, 11>> eps = { 
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

    blitz::Array<double, 1> result(n["t"]);

    // compare statistics
    // height of the center of mass of cloud droplets
    for (int at = 0; at < n["t"]; ++at)
      result(at) = plotter.act_com_z_timestep(at * 60); 
    std::cout << "height of the center of mass: " << result;
    blitz::Array<double, 1> expected_result(data_com[opts_m.first].data(), 11, blitz::neverDeleteData);
    blitz::Array<double, 1> epsilon(eps[opts_m.first].data(), 11, blitz::neverDeleteData);
    err_flag = errcheck(result, expected_result, epsilon) || err_flag;

    // average cloud water mixing ratio in cloudy cells
    for (int at = 0; at < n["t"]; ++at)
      result(at) = (plotter.cloud_ract_stats_timestep(at * 60)).first;
    std::cout << "average cloud water mixing ratio in cloudy cells: " << result;
    expected_result = blitz::Array<double, 1>(data_avg[opts_m.first].data(), 11, blitz::neverDeleteData);
    err_flag = errcheck(result, expected_result, epsilon) || err_flag;

    // std_dev of cloud water mixing ratio in cloudy cells
    for (int at = 0; at < n["t"]; ++at)
      result(at) = (plotter.cloud_ract_stats_timestep(at * 60)).second;
    std::cout << "std dev of cloud water mixing ratio in cloudy cells: " << result;
    expected_result = blitz::Array<double, 1>(data_std_dev_rc[opts_m.first].data(), 11, blitz::neverDeleteData);
    err_flag = errcheck(result, expected_result, epsilon) || err_flag;


    // average concentration of activated droplets in cloudy cells
    for (int at = 0; at < n["t"]; ++at)
    {
    }
  }

  if(err_flag)
    error_macro("error in one of the statistics");    
}

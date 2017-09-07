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
//  if(any(rel_err > 0.))
//  {
    std::cerr << "ERROR" << std::endl;
    std::cerr << "expected result: " << expected_result;
    std::cerr << "relative error minus precision: " << rel_err;
    return 1;
//  }
//  else
//    return 0;
}

struct test_data
{
  string name;
  unordered_map<string, std::array<double, 11>> expected,
                                                epsilon;
};

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

  // container for the expecged result and the epslion (req precision) for each tested statistic
  vector<test_data> tests;
  // populate it - height of the center of mass
  tests.push_back({
    // test name
    "com_z",
    // expected values map
    {
      {"blk_1m", {{0., 841.758, 917.667, 1022.48, 1134.36, 1232.02, 1308.47, 1364.44, 1412.22, 1476.16, 1521.43}}},
      {"lgrngn", {{0., 842.771, 916.496, 1020.23, 1129.24, 1226.24, 1299.9,  1351.78, 1398.83, 1443.72, 1467.06}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 1e-1}} } 
    }
  });

  // average mass mixing ratio of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "rc_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.185381, 0.357674, 0.518615, 0.635203, 0.707202, 0.733656, 0.686596, 0.556995, 0.425793 }}},
      {"lgrngn", {{ 0, 0, 0.179877, 0.374551, 0.552171, 0.704754, 0.790675, 0.83, 0.83, 0.8, 0.73}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-1}} }
    }
  });

  // standard deviation of the mass mixing ratio of activated droplets in cloudy cells
  tests.push_back({
    // test name
    "rc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.0189641, 0.0669555, 0.137429, 0.210689, 0.266551, 0.30766, 0.333871, 0.310118, 0.232963}}},
      {"lgrngn", {{0, 0, 0.0551891, 0.12, 0.202856, 0.290047, 0.375072, 0.442025, 0.485523, 0.525853, 0.523559}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 5e-2, 5e-2, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1}} } 
    }
  });

  // average concentration of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "actconc_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.185381, 0.357674, 0.518615, 0.635203, 0.707202, 0.733656, 0.686596, 0.556995, 0.425793 }}},
      {"lgrngn", {{ 0, 0, 0.179877, 0.374551, 0.552171, 0.704754, 0.790675, 0.791144, 0.73452, 0.648194, 0.548742}}}  
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 10e-2, 15e-2, 5e-1, 5e-1, 1.5}} } 
    }
  });

  // stqandard deviation of the concentration of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "actconc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.185381, 0.357674, 0.518615, 0.635203, 0.707202, 0.733656, 0.686596, 0.556995, 0.425793 }}},
      {"lgrngn", {{ 0, 0, 0.179877, 0.374551, 0.552171, 0.704754, 0.790675, 0.791144, 0.73452, 0.648194, 0.548742}}}  
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 10e-2, 15e-2, 5e-1, 5e-1, 1.5}} } 
    }
  });

  // out dir
  unordered_map<string, string> tmp_out = { {"blk_1m", "tmp_out_blk_1m"}, {"lgrngn", "tmp_out_lgrngn"}};

  // flag indicating if there was an error in any statistic
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
    for(auto &test: tests)
    {
      // calculate the statistic from the simulation at each output step
      for (int at = 0; at < n["t"]; ++at)
      {
        if(test.name == "com_z")
          result(at) = plotter.act_com_z_timestep(at * 60); 
        else if(test.name == "rc_avg")
          result(at) = (plotter.cloud_ract_stats_timestep(at * 60)).first;
        else if(test.name == "rc_std_dev")
          result(at) = (plotter.cloud_ract_stats_timestep(at * 60)).second;
        else if(test.name == "actconc_avg")
          result(at) = (plotter.cloud_actconc_stats_timestep(at * 60)).first;
        else if(test.name == "actconc_std_dev")
          result(at) = (plotter.cloud_actconc_stats_timestep(at * 60)).second;
      }

      // output the result
      cout << test.name<< " : " << result;
      // read in expected result and the epsilon
      blitz::Array<double, 1> expected_result(test.expected[opts_m.first].data(), 11, blitz::neverDeleteData);
      blitz::Array<double, 1> epsilon(test.epsilon[opts_m.first].data(), 11, blitz::neverDeleteData);
      // check if there is an agreement
      err_flag = errcheck(result, expected_result, epsilon) || err_flag;
    }
  }

  if(err_flag)
    error_macro("error in one of the statistics");    
}

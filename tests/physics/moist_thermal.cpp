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
    std::cerr << "expected result: " << expected_result;
    std::cerr << "relative error minus precision: " << rel_err;
  if(any(rel_err > 0.))
  {
    std::cerr << "ERROR" << std::endl;
    return 1;
  }
  else
    return 0;
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
    "--outfreq=60 --nt=600 --spinup=0 --dt=1 --nx=181 --nz=121 --case=moist_thermal";
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
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 2e-2, 3e-2, 6e-2, 10e-2}} } 
    }
  });

  // average mass mixing ratio of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "rc_avg",
    // expected values map
    {
      {"blk_1m", {{ 0, 0, 0.185215, 0.355712, 0.520675, 0.63674, 0.720535, 0.758259, 0.740406, 0.624473, 0.527215 }}},
      {"lgrngn", {{ 0, 0, 0.179877, 0.374551, 0.552171, 0.704754, 0.790675, 0.83, 0.83, 0.8, 0.73}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 3e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{7e-2, 7e-2, 7e-2, 7e-2, 7e-2, 7e-2, 7e-2, 7e-2, 7e-2, 15e-2, 20e-2}} }
    }
  });

  // standard deviation of the mass mixing ratio of activated droplets in cloudy cells
  tests.push_back({
    // test name
    "rc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.0188977, 0.0701371, 0.135743, 0.211222, 0.258296, 0.314703, 0.34613, 0.35574, 0.327544}}},
      {"lgrngn", {{0, 0, 0.0551891, 0.12, 0.202856, 0.290047, 0.375072, 0.442025, 0.485523, 0.525853, 0.523559}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 3e-2, 5e-2, 7e-2, 10e-2}} },
      {"lgrngn", {{1e-2, 25e-2, 25e-2, 25e-2, 25e-2, 25e-2, 25e-2, 25e-2, 25e-2, 25e-2, 25e-2}} } 
    }
  });

  // average concentration of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "actconc_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0, 89.2595, 87.5162, 86.9183, 86.0386, 83.6924, 80.9611, 78.7904, 76.5399, 73.252}}}  

    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 7e-2, 10e-2, 10e-2}} } 
    }
  });

  // stqandard deviation of the concentration of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "actconc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0, 27.902, 30.7292, 31.9791, 30.947, 32.4312, 34.3208, 35.8673, 35.6399, 34.2308}}}  

    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{10e-2, 10e-2, 20e-2, 15e-2, 15e-2, 20e-2, 15e-2, 10e-2, 10e-2, 15e-2, 20e-2}} } 
    }
  });

  // average supersaturation in supersaturated cells
  tests.push_back({
    // test name
    "supersat_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0.476906, 0.428138, 0.374729, 0.333207, 0.277702, 0.224544, 0.213603, 0.251973, 0.298691, 0.351449}}}

    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 2e-2, 4e-2, 5e-2, 9e-2, 10e-2, 15e-2, 30e-2, 35e-2, 35e-2, 50e-2}} } 
    }
  });

  // std dev of supersaturation in supersaturated cells
  // disabled due to too low precision - would require higher sd_conc?
  /*
  tests.push_back({
    // test name
    "supersat_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0.121373, 0.212959, 0.204964, 0.267225, 0.274048, 0.272687, 0.367227, 0.394839, 0.391756, 0.469211}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{40e-2, 40e-2, 40e-2, 40e-2, 40e-2, 40e-2, 40e-2, 40e-2, 40e-2, 50e-2, 60e-2}} } // high, makes any sense? 
    }
  });
  */

  // average number of SDs in cloudy cells
  tests.push_back({
    // test name
    "sdconc_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0, 16.138, 15.7586, 15.5091, 15.3118, 15.1084, 15.0345, 15.098, 15.1314, 14.9388}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2}} } 
    }
  });

  // std_dev of number of SDs in cloudy cells
  tests.push_back({
    // test name
    "sdconc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0, 2.82919, 2.9137, 3.16832, 3.27304, 3.21091, 3.56438, 3.4866, 3.60563, 3.73938}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{16e-2, 16e-2, 16e-2, 16e-2, 16e-2, 16e-2, 16e-2, 16e-2, 16e-2, 16e-2, 20e-2}} } 
    }
  });

  // average mean radius of activated droplets
  tests.push_back({
    // test name
    "meanr_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0, 7.63714, 9.70857, 11.1745, 12.1136, 12.4604, 12.7156, 12.6461, 12.3354, 11.9144}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 5e-2, 8e-2}} } 
    }
  });

  // average std dev of radius of activated droplets
  tests.push_back({
    // test name
    "stddevr_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0, 0.966154, 1.22183, 1.22192, 1.18918, 1.296, 1.43948, 1.74801, 2.06894, 2.31939}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{15e-2, 15e-2, 15e-2, 15e-2, 15e-2, 15e-2, 15e-2, 15e-2, 15e-2, 17e-2, 20e-2}} } 
    }
  });

  // cloud fraction
  tests.push_back({
    // test name
    "clfrac",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.0232409, 0.024839, 0.0257066, 0.0267568, 0.0264828, 0.0253413, 0.0230583, 0.0200904, 0.0155244}}},
      {"lgrngn", {{0, 0, 0.0210493, 0.0241998, 0.024976, 0.0257522, 0.0259349, 0.0247477, 0.0239258, 0.0208666, 0.0168029}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 10e-2}} },
      {"lgrngn", {{5e-2, 5e-2, 7e-2, 5e-2, 5e-2, 5e-2, 5e-2, 10e-2, 15e-2, 15e-2, 15e-2}} } 
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
        else if(test.name == "supersat_avg")
          result(at) = (plotter.positive_supersat_stats_timestep(at * 60)).first;
        else if(test.name == "supersat_std_dev")
          result(at) = (plotter.positive_supersat_stats_timestep(at * 60)).second;
        else if(test.name == "sdconc_avg")
          result(at) = (plotter.cloud_sdconc_stats_timestep(at * 60)).first;
        else if(test.name == "sdconc_std_dev")
          result(at) = (plotter.cloud_sdconc_stats_timestep(at * 60)).second;
        else if(test.name == "meanr_avg")
          result(at) = (plotter.cloud_meanr_stats_timestep(at * 60)).first;
        else if(test.name == "stddevr_avg")
          result(at) = (plotter.cloud_stddevr_stats_timestep(at * 60)).first;
        else if(test.name == "clfrac")
        {
          Plotter_t::arr_t ract(plotter.h5load_ract_timestep(at * 60));
          ract = iscloudy_rc(ract);
          result(at) = blitz::mean(ract); 
        }
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

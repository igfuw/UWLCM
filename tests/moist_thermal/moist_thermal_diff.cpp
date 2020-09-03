// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <sstream> // std::ostringstream

#include "common.hpp"
#include <drawbicyc/PlotterMicro.hpp>
#include <drawbicyc/common_filters.hpp>

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
  // container for the expecged result and the epslion (req precision) for each tested statistic
  vector<test_data> tests;
  // populate it - height of the center of mass
  tests.push_back({
    // test name
    "com_z",
    // expected values map
    {
      {"blk_1m", {{0., 841.758, 917.667, 1022.48, 1134.36, 1232.02, 1308.47, 1364.44, 1412.22, 1476.16, 1521.43}}},
      {"lgrngn", {{0, 844.156, 919.85, 1024.4, 1135.51, 1233.24, 1307.44, 1354.7, 1386.39, 1405.33, 1417.56}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2}} } 
    }
  });

  // average mass mixing ratio of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "rc_avg",
    // expected values map
    {
      {"blk_1m", {{ 0, 0    , 0.19, 0.37    , 0.533   , 0.675   , 0.74    , 0.74, 0.66, 0.5, 0.39 }}},
      {"lgrngn", {{0, 0.108071, 0.192615, 0.391706, 0.58865, 0.734293, 0.814996, 0.8343, 0.801823, 0.762012, 0.716567}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 3e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 3e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 2e-2, 2e-2, 4e-2, 6e-2, 8e-2}} }
    }
  });

  // standard deviation of the mass mixing ratio of activated droplets in cloudy cells
  tests.push_back({
    // test name
    "rc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.022    , 0.076, 0.154   , 0.206, 0.271, 0.332, 0.326, 0.262, 0.228}}},
      {"lgrngn", {{0, 0.0058853, 0.0623114, 0.139442, 0.22448, 0.310132, 0.407476, 0.466262, 0.508034, 0.520107, 0.496754}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 1e-2, 1e-1, 2e-2, 4e-2, 1e-2, 2e-2, 3e-2, 5e-2, 7e-2, 10e-2}} },
      {"lgrngn", {{1e-2, 5e-2, 3e-2, 3e-2, 3e-2, 3e-2, 4e-2, 5e-2, 5e-2, 7e-2, 11e-2}} } 
    }
  });

  // average concentration of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "actconc_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 115.683, 88.5124, 86.5786, 85.1575, 84.6427, 80.8379, 78.769, 72.4859, 71.4538, 68.7489}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 2e-2, 2e-2, 2e-2, 4e-2, 5e-2, 6e-2}} } 
    }
  });

  // stqandard deviation of the concentration of activated dropletes in cloudy cells
  tests.push_back({
    // test name
    "actconc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0.313923, 28.015, 29.9676, 31.7979, 33.1697, 35.4694, 36.5521, 36.0021, 36.4161, 34.7416}}} 
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 2e-2, 2e-2, 2e-2, 3e-2, 4e-2, 6e-2, 8e-2, 8e-2, 10e-2, 15e-2}} } 
    }
  });

  // average supersaturation in supersaturated cells
  tests.push_back({
    // test name
    "supersat_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0., 0.513721, 0.455369, 0.408629, 0.367086, 0.305751, 0.269151, 0.252091, 0.266137, 0.32337, 0.435048}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 2e-2, 4e-2, 4e-2, 4e-2, 4e-2, 9e-2, 10e-2, 20e-2, 35e-2}} } 
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
      {"lgrngn", {{10e-2, 10e-2, 10e-2, 10e-2, 10e-2, 10e-2, 10e-2, 10e-2, 10e-2, 10e-2, 10e-2}} } // high, makes any sense? 
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
      {"lgrngn", {{0, 14, 16.0447, 15.6962, 15.4722, 15.2388, 15.115, 15.0581, 15.0314, 14.9429, 15.0769}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 2e-2, 2e-2, 2e-2}} } 
    }
  });

  // std_dev of number of SDs in cloudy cells
  tests.push_back({
    // test name
    "sdconc_std_dev",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 1, 2.93024, 2.96181, 3.06923, 3.18366, 3.29947, 3.3961, 3.82987, 3.61984, 3.78123}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1.5e-2, 3e-2, 4e-2, 7e-2, 9e-2, 10e-2, 12e-2, 12e-2}} } 
    }
  });

  // average mean radius of activated droplets
  tests.push_back({
    // test name
    "meanr_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 5.77758, 7.78731, 9.90543, 11.4234, 12.1887, 12.6271, 12.6339, 12.6476, 12.3377, 12.1446}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1.5e-2, 2e-2, 3e-2, 3e-2}} } 
    }
  });

  // average std dev of radius of activated droplets
  tests.push_back({
    // test name
    "stddevr_avg",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
      {"lgrngn", {{0, 0.991351, 0.909609, 1.18174, 1.18827, 1.27298, 1.38948, 1.70126, 1.87223, 2.12477, 2.23955}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 2e-1}} },
      {"lgrngn", {{1e-2, 3e-2, 7e-2, 6e-2, 6e-2, 6e-2, 8e-2, 8e-2, 8e-2, 8e-2, 12e-2}} } 
    }
  });

  // cloud fraction
  tests.push_back({
    // test name
    "clfrac",
    // expected values map
    {
      {"blk_1m", {{0, 0, 0.0232409, 0.0246, 0.0257979, 0.0256, 0.0255, 0.0239, 0.0214, 0.0191, 0.0123}}},
      {"lgrngn", {{0, 9.132e-05, 0.0214145, 0.0237432, 0.0246564, 0.0254326, 0.02525, 0.0236519, 0.0215515, 0.0178531, 0.0142916}}}
    },
    // epsilons map
    {
      {"blk_1m", {{1e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 2e-2, 5e-2, 10e-2}} },
      {"lgrngn", {{1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 2e-2, 3e-2, 4e-2, 6e-2, 8e-2}} } 
    }
  });

  // flag indicating if there was an error in any statistic
  bool err_flag = 0;

  for (auto &opts_m : opts_micro)
  {
    // instantiate plotter to read in output
    using Plotter_t = PlotterMicro_t<2>;
    Plotter_t plotter(outdir.at(opts_m.first), opts_m.first);
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

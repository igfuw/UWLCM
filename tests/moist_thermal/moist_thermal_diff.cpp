// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <sstream> // std::ostringstream

#include "common.hpp"
#include <UWLCM_plotters/PlotterMicro.hpp>
#include <UWLCM_plotters/common_filters.hpp>

using barr1d = blitz::Array<double, 1>;
using namespace std;

bool errcheck(barr1d result, barr1d expected_result, barr1d epsilon)
{
  barr1d rel_err(result.shape());
  rel_err = where(expected_result > 0, abs(result - expected_result) / expected_result - epsilon, 0);
    std::cerr << "expected result: " << expected_result;
    std::cerr << "error tolerance: " << epsilon;
    std::cerr << "relative error minus tolerance: " << rel_err;
  if(any(rel_err > 0.))
  {
    std::cerr << "ERROR" << std::endl;
    return 1;
  }
  else
    return 0;
}

int main(int ac, char** av)
{

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
    for(auto stat_name: stat_names)
    {
      // calculate the statistic from the simulation at each output step
      for (int at = 0; at < n["t"]; ++at)
      {
        if(stat_name == "com_z")
          result(at) = plotter.act_com_z_timestep(at * 60); 
        else if(stat_name == "rc_avg")
          result(at) = (plotter.cloud_ract_stats_timestep(at * 60)).first;
        else if(stat_name == "rc_std_dev")
          result(at) = (plotter.cloud_ract_stats_timestep(at * 60)).second;
        else if(stat_name == "actconc_avg")
          result(at) = (plotter.cloud_actconc_stats_timestep(at * 60)).first;
        else if(stat_name == "actconc_std_dev")
          result(at) = (plotter.cloud_actconc_stats_timestep(at * 60)).second;
        else if(stat_name == "supersat_avg")
          result(at) = (plotter.positive_supersat_stats_timestep(at * 60)).first;
        else if(stat_name == "supersat_std_dev")
          result(at) = (plotter.positive_supersat_stats_timestep(at * 60)).second;
        else if(stat_name == "sdconc_avg")
          result(at) = (plotter.cloud_sdconc_stats_timestep(at * 60)).first;
        else if(stat_name == "sdconc_std_dev")
          result(at) = (plotter.cloud_sdconc_stats_timestep(at * 60)).second;
        else if(stat_name == "meanr_avg")
          result(at) = (plotter.cloud_meanr_stats_timestep(at * 60)).first;
        else if(stat_name == "stddevr_avg")
          result(at) = (plotter.cloud_stddevr_stats_timestep(at * 60)).first;
        else if(stat_name == "clfrac")
        {
          Plotter_t::arr_t ract(plotter.h5load_ract_timestep(at * 60));
          ract = iscloudy_rc(ract);
          result(at) = blitz::mean(ract); 
        }
      }

      // output the result
      cout << stat_name<< " : " << result;

      // read reference data
      std::ifstream fref("../../moist_thermal/refdata/stats.txt");

      std::string micro;
      barr1d mean, std_dev;
      // find the line with the stat
      bool found=0;
      std::string line;
      while(getline( fref, line ) ) 
      {
        if( line.find( stat_name ) != string::npos ) // data found
        {
          // find the line with the micro
          for(int i=0; i<2; ++i) // done twice: two type of micro
          {
            fref >> micro;
            if(micro == opts_m.first)
            {
              fref >> mean;
              getline(fref, line); // blitz reading arrays doesnt move to the next line, need to do it manually here
              fref >> std_dev;
              getline(fref, line); // blitz reading arrays doesnt move to the next line, need to do it manually here
              found=1;
            }
          }
          break;
        }
      }
      if(!found)
        throw runtime_error("One of the stats not found in the refrence data file.");

      // check if there is an agreement
      cerr << result;
      cerr << mean;
      cerr << std_dev;
//      err_flag = errcheck(result, test.mean[opts_m.first], test.std_dev[opts_m.first]) || err_flag;
    }
  }

  if(err_flag)
    error_macro("error in one of the statistics");    
}

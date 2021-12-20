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

const int tolerance = 3; // throw an error if the result is outside of: expected_value +/- tolerance * standard_deviation

bool errcheck(barr1d result, barr1d mean, barr1d std_dev)
{
  cout << "result: " << result;
  cout << "reference mean: " << mean;
  cout << "reference standard deviation: " << std_dev;
  cout << "reference relative standard deviation [%]: " << barr1d(where(mean > 0, std_dev / mean * 100., 0));
  barr1d diff(result.shape());
  diff = result - mean;
  cout << "(result - mean) / mean [%]: " << barr1d(where(mean > 0, diff / mean * 100., 0));
  cout << "abs(result - mean) / std_dev: " << barr1d(where(std_dev > 0, abs(diff) / std_dev, 0));
  diff = (abs(diff) - tolerance * std_dev) / std_dev;

  if(any(diff > 0.))
  {
    cerr << "error flag: " << diff;
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
      cout << stat_name<< endl;

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
            fref >> mean;
            getline(fref, line); // blitz reading arrays doesnt move to the next line, need to do it manually here
            fref >> std_dev;
            getline(fref, line); // blitz reading arrays doesnt move to the next line, need to do it manually here
            found=1;
            if(micro == opts_m.first)
              break;
          }
          break;
        }
      }
      if(!found)
        throw runtime_error("One of the stats not found in the refrence data file.");

      // check if there is an agreement
      err_flag = errcheck(result, mean, std_dev) || err_flag;
    }
  }

  if(err_flag)
    error_macro("error in one of the statistics");    
}

// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <sstream> // std::ostringstream
#include <fstream>
#include <cmath>

#include "common.hpp"
#include <UWLCM_plotters/PlotterMicro.hpp>
#include <UWLCM_plotters/common_filters.hpp>

using barr1d = blitz::Array<double, 1>;
using namespace std;


// https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison?page=1&tab=votes#tab-top
//const double dbl_epsilon = 1000*std::numeric_limits<double>::epsilon();
const double dbl_epsilon = 1e-5; 
bool less_than(double a, double b)
{
  return (b - a) > ( (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * dbl_epsilon);
}
bool greater_than(double a, double b)
{
  return (a - b) > ( (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * dbl_epsilon);
}
BZ_DECLARE_FUNCTION2_RET(less_than, bool);
BZ_DECLARE_FUNCTION2_RET(greater_than, bool);

// check if modeled result is within [min,max] range from 1000 test runs
// note: there is no randomness with blk_1m microphysics
// note2: plenty of precision lost when storing refdata in files, so a relative precision of 1e-5 is required
bool errcheck_minmax(barr1d result, barr1d min, barr1d max)
{
  barr1d err_flag(result.shape());

  barr1d tolerance(result.shape());
  tolerance = (max - min)*0.02;

  err_flag = where(
    less_than(result,min-tolerance), 
      1,
      where(greater_than(result, max+tolerance), 1 , 0)
  ); 

  if(any(err_flag > 0.))
  {
    cerr << "min: " << min;
    cerr << "result: " << result;
    cerr << "max: " << max;
    cerr << "error flag: " << err_flag;
    return 1;
  }
  else
    return 0;
}

// check if the result is within mean +/- tolerance * std_dev, where mean and std_dev come from 1000 test runs
// note: currently not used
bool errcheck_stddev(barr1d result, barr1d mean, barr1d std_dev)
{
  const int tolerance = 3;
  barr1d diff(result.shape());
  diff = (abs(diff) - tolerance * std_dev);

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
    // instantiate plotter to read output
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

      //different reference file for bulk micro compiled with MPI, see refdata/readme.md for details
      string reffile_name = 
        opts_m.first == "lgrngn" ? "../../moist_thermal/refdata/stats_lgrngn_ens_1000.txt" :  // lgrngn 
          plotter.map["MPI_compiler"] ? "../../moist_thermal/refdata/stats_mpi_blk_ens_1.txt": // bulk with mpi
            "../../moist_thermal/refdata/stats_blk_ens_1.txt"; // bulk without mpi

      cout << "checking " << stat_name << " reference file: " << reffile_name << endl;
      std::ifstream fref(reffile_name);

      std::string micro;
      barr1d mean, std_dev, min, max;
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
            getline(fref, line); 
            fref >> min;
            getline(fref, line); 
            fref >> max;
            getline(fref, line); 
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
//      err_flag = (errcheck_stddev(result, mean, std_dev) && errcheck_minmax(result, min, max)) || err_flag;
      err_flag = errcheck_minmax(result, min, max) || err_flag;
    }
  }

  if(err_flag)
    error_macro("Error in one of the statistics. Make sure that libcloudph++ and UWLCM are compiled with the same flags as in the moist_thermal job in https://github.com/igfuw/UWLCM/blob/master/.github/workflows/test_uwlcm_hlpr.yml");
}

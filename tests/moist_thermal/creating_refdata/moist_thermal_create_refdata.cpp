// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <sstream> // std::ostringstream
#include <fstream>

#include "../common.hpp"
#include <UWLCM_plotters/PlotterMicro.hpp>
#include <UWLCM_plotters/common_filters.hpp>

using arr1d = blitz::Array<double, 1>;

int main(int ac, char** av)
{
  const int ensemble_size = ac-1;

  unordered_map<string, unordered_map<string, unordered_map<string, arr1d>>> micro_stat_dir_data_map; // microname -> statname -> dirname -> data
  unordered_map<string, unordered_map<string, arr1d>> micro_stat_mean_map; // microname -> statname -> mean
  unordered_map<string, unordered_map<string, arr1d>> micro_stat_StdDev_map; // microname -> statname -> standard deviation

  int nt = -1; // number of output timesteps

  // loop over input dirs, each with tmp_out_lgrngn and tmp_out_blk_1m subdirectories
  for (char **a = av+1; a != av+ac ; a++) 
  {
    const string dirname = *a;

    // loop over blk_1m and lgrngn micro
    for (auto &opts_m : opts_micro)
    {
      // instantiate plotter to read in output
      using Plotter_t = PlotterMicro_t<2>;
      Plotter_t plotter(dirname+"/"+outdir.at(opts_m.first), opts_m.first);
  
      // read const.h5
      auto& n = plotter.map;
  
      if(nt < 0) nt = n["t"];
      else if (n["t"] != nt) throw runtime_error("Input directories have different numbers of output timesteps");

      arr1d result(nt);

  
      for(auto &stat_name: stat_names)
      {
        // calculate the statistic from the simulation at each output step
        for (int at = 0; at < nt; ++at)
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

        pair<string, arr1d> dirname_data(dirname, result.copy());
        micro_stat_dir_data_map[opts_m.first][stat_name].insert(dirname_data);
      }
    }
  }

  arr1d stat_mean(nt),
        stat_std_dev(nt),
        stat_min(nt),
        stat_max(nt);

  // calculate min/max, mean and std dev of stats from the ensemble
  // and store in a file
  // note: it is assumed that the micro/stats/etc. order does not change
  ofstream outf("stats_ens_"+to_string(ensemble_size)+".txt");

  for(auto &stat_name: stat_names)
  {
    outf << stat_name << endl;
    for (auto &opts_m : opts_micro)
    {
      outf << opts_m.first << endl;

      stat_mean=0;
      for (char **a = av+1; a != av+ac ; a++) 
        stat_mean += micro_stat_dir_data_map[opts_m.first][stat_name][string(*a)];
      stat_mean /= ensemble_size;
      outf << stat_mean;

      stat_std_dev=0;
      for (char **a = av+1; a != av+ac ; a++) 
        stat_std_dev += blitz::pow(micro_stat_dir_data_map[opts_m.first][stat_name][string(*a)] - stat_mean, 2.);
      stat_std_dev = blitz::sqrt(stat_std_dev/ensemble_size);
      outf << stat_std_dev;
      
      stat_min =  1e20; // crap
      for (char **a = av+1; a != av+ac ; a++) 
        stat_min = blitz::where(micro_stat_dir_data_map[opts_m.first][stat_name][string(*a)] < stat_min, micro_stat_dir_data_map[opts_m.first][stat_name][string(*a)], stat_min);
      outf << stat_min;
      
      stat_max =  -1e20; // crap
      for (char **a = av+1; a != av+ac ; a++) 
        stat_max = blitz::where(micro_stat_dir_data_map[opts_m.first][stat_name][string(*a)] > stat_max, micro_stat_dir_data_map[opts_m.first][stat_name][string(*a)], stat_max);
      outf << stat_max;

    }
  }
}

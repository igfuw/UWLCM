#include "common.hpp"
#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include "plots.hpp"

const double D = 3.75e-6; //[1/s], ugly, large-scale horizontal wind divergence TODO: read from model output

template<class Plotter_t>
void plot_series(Plotter_t plotter, Plots plots)
{

  auto& n = plotter.map;
  for(auto elem : n)
  {
     std::cout << elem.first << " " << elem.second << std::endl;
  }
  Gnuplot gp;
  string file = plotter.file + "_series.svg";
  int hor = min<int>(plots.series.size(), 4);
  int ver = double(plots.series.size()) / 4. + 0.99999;
  init_prof(gp, file, ver, hor); 

  string prof_file = plotter.file + "_series.dat";
  std::ofstream oprof_file(prof_file);

  // read in density
  auto tmp = plotter.h5load(plotter.file + "/const.h5", "G");
  typename Plotter_t::arr_t rhod(tmp);
  typename Plotter_t::arr_t rtot(rhod.shape());

  typename Plotter_t::arr_t res_tmp(rhod.shape());


  // read opts
  po::options_description opts("profile plotting options");
  opts.add_options()
    ("series_start", po::value<int>()->default_value(0) , "time in sec when we start drawin series")
    ("series_end", po::value<int>()->default_value(0) , "time in sec when we end drawing series")
  ;
  po::variables_map vm; 
  handle_opts(opts, vm);

  int first_timestep =  vm["series_start"].as<int>() / n["dt"] / n["outfreq"];
  int last_timestep =  vm["series_end"].as<int>() / n["dt"] / n["outfreq"];
  if(last_timestep == 0) last_timestep = n["t"]-1;

  Array<double, 1> res_prof(last_timestep - first_timestep + 1);
  Array<double, 1> res_prof_std_dev(last_timestep - first_timestep + 1);
  Array<double, 1> res_pos(last_timestep - first_timestep + 1),
    com_N_c(last_timestep - first_timestep + 1), // particles concentration at the center of mass
    com_miu(last_timestep - first_timestep + 1); // to keep mean particle radius at the center of mass
  Array<int, 1> com_z_idx(last_timestep - first_timestep + 1), 
    com_x_idx(last_timestep - first_timestep + 1); // index of the center of mass cell

  // save time steps to the series file
  oprof_file << plotter.timesteps;


  for (auto &plt : plots.series)
  {
    bool plot_std_dev = 0;
    res_prof_std_dev = 0;
    res_prof = 0;
    res_pos = 0;

    std::ifstream f_precip(plotter.file + "/prec_vol.dat");
    std::string row;
    double prec_vol = 0.;
    double prec_vol_prev;

    for (int at = first_timestep; at <= last_timestep; ++at) // TODO: mark what time does it actually mean!
    {
      res_pos(at) = at * n["outfreq"] * n["dt"] / 3600.;
      // store accumulated precip volume
      prec_vol_prev = prec_vol;
      // read in accumulated precipitation volume
      // wet precip vol is the 9th value, followed by druy vol and an empty line
      for(int i=0; i<9; ++i)
        std::getline(f_precip, row);
      sscanf(row.c_str(), "%*d %lf", &prec_vol);
      std::getline(f_precip, row);
      std::getline(f_precip, row);

      if (plt == "clfrac")
      {
/*
        try
        {
          // cloud fraction (cloudy if q_c > 0.1 g/kg)
          // read activated droplets mixing ratio to res_tmp 
          auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = iscloudy_rc(snap); // find cells with rc>1e-5
          res_prof(at) = blitz::mean(res_tmp); 
        }
*/
        try
        {
          auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; //g/kg
          typename Plotter_t::arr_t snap(tmp); 
          snap += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; //g/kg
          snap *= rhod; // water per cubic metre (should be wet density...)
          plotter.k_i = blitz::sum(snap, plotter.LastIndex) * n["dz"]; // LWP [g/m2] in the column 
          plotter.k_i = where(plotter.k_i > 20 , 1 , 0); // cloudiness as in Ackermann et al. 
          res_prof(at) = blitz::mean(plotter.k_i);
        }
        catch(...){;}
      }
      // max RH in the domain
      else if (plt == "RH_max")
      {
        try
        {
          // read RH 
          auto tmp = plotter.h5load_timestep("RH", at * n["outfreq"]);

          typename Plotter_t::arr_t snap(tmp);
          res_prof(at) = blitz::max(snap);
        }
        catch(...) {;}
      }
      // r_act averaged over cloudy cells
      else if (plt == "ract_avg")
      {
        try
        {
          auto stats = plotter.cloud_ract_stats_timestep(at * n["outfreq"]);
          res_prof(at) = stats.first;
          res_prof_std_dev(at) = stats.second;
        }
        catch(...) {;}
      }
      // averagew sd_conc in clloudy cells
      else if (plt == "sd_conc_avg")
      {
        try
        {
          auto stats = plotter.cloud_sdconc_stats_timestep(at * n["outfreq"]);
          res_prof(at) = stats.first;
          res_prof_std_dev(at) = stats.second;
        }
        catch(...) {;}
      }
      // average activated sd_conc in clloudy cells
      else if (plt == "sd_conc_act_avg")
      {
        try
        {
          auto stats = plotter.cloud_sdconc_act_stats_timestep(at * n["outfreq"]);
          res_prof(at) = stats.first;
          res_prof_std_dev(at) = stats.second;
        }
        catch(...) {;}
      }
      else if (plt == "tot_water")
      {
        try
        {
/*
          {
            auto tmp = plotter.h5load_timestep("aerosol_rw_mom3", at * n["outfreq"]) * 4./3. * 3.1416 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_prof(at) = blitz::mean(snap);
          }
*/
          {
            auto tmp = plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]) * 4./3. * 3.1416 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_prof(at) += blitz::mean(snap);
          }
          {
            auto tmp = plotter.h5load_timestep("rain_rw_mom3", at * n["outfreq"]) * 4./3. * 3.1416 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_prof(at) += blitz::mean(snap);
          }
          {
            auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]);
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_prof(at) += blitz::mean(snap);
          } 
        }
        catch(...) {;}
      }
      else if (plt == "ract_com")
      {
	// center of mass of activated droplets
        try
        {
          res_prof(at) = plotter.act_com_z_timestep(at * n["outfreq"]);
        }        
        catch(...) {;}
      }
      else if (plt == "com_vel")
      {
	// vertical velocity at the center of mass of activated droplets
        try
        {
          auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          typename Plotter_t::arr_t snap2(tmp);
          typename Plotter_t::arr_t snap3(tmp);
          
          snap2 = snap2 * plotter.LastIndex;
          snap3 = snap3 * blitz::tensor::i;
          if(blitz::sum(snap) > 1e-3)
          {
            int z_idx = blitz::sum(snap2) / blitz::sum(snap); 
            int x_idx = blitz::sum(snap3) / blitz::sum(snap); 
            auto tmp2 = plotter.h5load_timestep("w", at * n["outfreq"]);
            typename Plotter_t::arr_t snap_mom(tmp2);
            res_prof(at) = snap_mom(x_idx, z_idx);
          } 
          else 
            res_prof(at) = 0.;
        }
        catch(...) {;}
      }
      else if (plt == "com_supersat")
      {
	// supersaturation at the center of mass of activated droplets
        try
        {
          auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          typename Plotter_t::arr_t snap2(tmp);
          typename Plotter_t::arr_t snap3(tmp);
          
          snap2 = snap2 * plotter.LastIndex;
          snap3 = snap3 * blitz::tensor::i;
          if(blitz::sum(snap) > 1e-3)
          {
            int z_idx = blitz::sum(snap2) / blitz::sum(snap); 
            int x_idx = blitz::sum(snap3) / blitz::sum(snap); 
            auto tmp2 = plotter.h5load_timestep("RH", at * n["outfreq"]);
            typename Plotter_t::arr_t snap_mom(tmp2);
            res_prof(at) = snap_mom(x_idx, z_idx) - 1;
          } 
          else 
            res_prof(at) = 0.;
        }
        catch(...) {;}
      }
      else if (plt == "com_mom0")
      {
	// 0th moment of rw distribution at the center of mass of activated droplets (particles concentration), 2D only
        try
        {
          auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          typename Plotter_t::arr_t snap2(tmp);
          typename Plotter_t::arr_t snap3(tmp);
          
          snap2 = snap2 * plotter.LastIndex;
          snap3 = snap3 * blitz::tensor::i;
          if(blitz::sum(snap) > 1e-3)
          {
            com_z_idx(at) = blitz::sum(snap2) / blitz::sum(snap); 
            com_x_idx(at) = blitz::sum(snap3) / blitz::sum(snap); 
            std::cout << at << ": (" << com_x_idx(at) << "," << com_z_idx(at) << ")" << std::endl;
            auto tmp2 = plotter.h5load_timestep("actrw_rw_mom0", at * n["outfreq"]);
            typename Plotter_t::arr_t snap_mom(tmp2);
            com_N_c(at) = snap_mom(com_x_idx(at), com_z_idx(at)); // 0th raw moment / mass [1/kg]
            snap_mom *= rhod; // now per m^3
            res_prof(at) = snap_mom(com_x_idx(at), com_z_idx(at));
          }
          else 
          {
            com_N_c(at) = 0.;
            res_prof(at) = 0.;
          }
        }
        catch(...) {;}
      }
      else if (plt == "com_mom1")
      {
        // mean droplet radius at the center of mass
        try
        {
          auto tmp = plotter.h5load_timestep("actrw_rw_mom1", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp); // 1st raw moment / mass [m / kg]
          if(com_N_c(at) > 0)
            res_prof(at) = snap(com_x_idx(at), com_z_idx(at)) / com_N_c(at);
          else
            res_prof(at) = 0.;
          com_miu(at) = res_prof(at); // mean radius [m]
        }
        catch(...) {;}
      }
      else if (plt == "com_mom2")
      {
        // std deviation of distribution of radius at center of mass
        try
        {
          auto tmp = plotter.h5load_timestep("actrw_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t zeroth_raw_mom(tmp); // 0th raw moment / mass [1 / kg]
          tmp = plotter.h5load_timestep("actrw_rw_mom1", at * n["outfreq"]);
          typename Plotter_t::arr_t first_raw_mom(tmp); // 1st raw moment / mass [m / kg]
          tmp = plotter.h5load_timestep("actrw_rw_mom2", at * n["outfreq"]);
          typename Plotter_t::arr_t second_raw_mom(tmp); // 2nd raw moment / mass [m^2 / kg]
          tmp = plotter.h5load_timestep("sd_conc", at * n["outfreq"]);
          typename Plotter_t::arr_t sd_conc(tmp); // number of SDs
          if(com_N_c(at) > 0)
          {
            double SD_no = sd_conc(com_x_idx(at), com_z_idx(at));
            if(SD_no > 1 && com_miu(at) > 0)
            {
              res_prof(at) = ( 
                SD_no / (SD_no - 1) /
                com_N_c(at) * (
                  second_raw_mom(com_x_idx(at), com_z_idx(at)) - 
                  2. * com_miu(at) * first_raw_mom(com_x_idx(at), com_z_idx(at)) + 
                  com_miu(at) * com_miu(at) * zeroth_raw_mom(com_x_idx(at), com_z_idx(at))
                )
              );
              
              // could not be true due to numerics?
              if(res_prof(at) > 0.) 
                res_prof(at) = sqrt(res_prof(at));
              else 
                res_prof(at) = 0.;
            }
          }
          else
            res_prof(at) = 0.;
        }
        catch(...) {;}
      }
      else if (plt == "com_sd_conc")
      {
        // number of SDs at the center of mass
        try
        {
          tmp = plotter.h5load_timestep("sd_conc", at * n["outfreq"]);
          typename Plotter_t::arr_t sd_conc(tmp); // number of SDs
          res_prof(at) = sd_conc(com_x_idx(at), com_z_idx(at));
        }
        catch(...) {;}
      }
      else if (plt == "th_com")
      {
	// center of mass of temp perturb
        try
        {
          auto tmp = plotter.h5load_timestep("th", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          
          res_tmp = is_th_prtrb(snap); // find cells with th>300.1
          snap *= res_tmp; // apply filter
          res_tmp = snap * plotter.LastIndex * n["dz"];
          
          if(blitz::sum(res_tmp) > 0.)
            res_prof(at) = blitz::sum(res_tmp) / blitz::sum(snap); 
          else
            res_prof(at) = 0.;
        }
        catch(...) {;}
      }
      else if (plt == "nc")
      {
	// cloud droplet (0.5um < r < 25 um) concentration
        try
        {
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap /= 1e6; // per cm^3
          snap *= rhod; // b4 it was per milligram
          res_prof(at) = blitz::mean(snap); 
        }
        catch(...) {;}
      }
      else if (plt == "cl_nc")
      {
	// cloud droplet (0.5um < r < 25 um) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          typename Plotter_t::arr_t snap2;
          snap2.resize(snap.shape());
          snap2=snap;
          snap = iscloudy(snap); // cloudiness mask
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_prof(at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_prof(at) = 0;
        }
        catch(...){;}
      }
      else if (plt == "cloud_base")
      {
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          snap(plotter.hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
          plotter.k_i = blitz::first((snap == 1), plotter.LastIndex); 
          auto cloudy_column = plotter.k_i.copy();
          cloudy_column = blitz::sum(snap, plotter.LastIndex);
          cloudy_column = where(cloudy_column > 0, 1, 0);
          plotter.k_i = where(cloudy_column == 0, 0, plotter.k_i);
          if(blitz::sum(cloudy_column) > 0)
            res_prof(at) = double(blitz::sum(plotter.k_i)) / blitz::sum(cloudy_column) * n["dz"];
          else
            res_prof(at) = 0;
        }
        catch(...){;}
      }
      // average concentration of activated droplets in cloudy cells
      else if (plt == "cloud_avg_act_conc")
      {
        try
        {
          auto stats = plotter.cloud_actconc_stats_timestep(at * n["outfreq"]);
          res_prof(at) = stats.first;
          res_prof_std_dev(at) = stats.second;
        }
        catch(...){;}
      }
      // average supersaturation in cells with S>0
      else if (plt == "cloud_avg_supersat")
      {
        try
        {
          auto stats = plotter.positive_supersat_stats_timestep(at * n["outfreq"]);
          res_prof(at) = stats.first;
          res_prof_std_dev(at) = stats.second;
        }
        catch(...){;}
      }
      // average radius of activated droplets in cloudy cells
      else if (plt == "cl_avg_cloud_rad")
      {
        try
        {
          auto stats = plotter.cloud_meanr_stats_timestep(at * n["outfreq"]);
          res_prof(at) = stats.first;
        }
        catch(...){;}
      }

      // average standard deviation of the acttivated droplets radius distribution in cloudy cells
      else if (plt == "cloud_avg_std_dev_act_rad")
      {
        try
        {
          auto stats = plotter.cloud_stddevr_stats_timestep(at * n["outfreq"]);
          res_prof(at) = stats.first;
        }
        catch(...){;}
      }

      else if (plt == "mass_dry")
      {
	// total dry mass
        double rho_dry = 1769; //[kg/m^3] - density of ammonium sulfate from wikipedia
        try
        {
          auto tmp = plotter.h5load_timestep("rd_rng000_mom3", at * n["outfreq"]) * 4./3. * 3.14 * rho_dry * 1e3;
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod * plotter.CellVol; // turn mixing ratio in g/kg to total mass in g
          res_prof(at) = blitz::sum(snap); 
        }
        catch(...) {;}
      }
      else if (plt == "surf_precip")
      {
        // surface precipitation [mm/day]
        try
        {
          res_prof(at) = plotter.calc_surf_precip(prec_vol - prec_vol_prev);
        }
        catch(...) {;}
      }
      else if (plt == "acc_precip")
      {
        // accumulated surface precipitation [mm]
        try
        {
          res_prof(at) = plotter.calc_acc_surf_precip(prec_vol);
        }
        catch(...) {;}
      }
      else if (plt == "lwp")
      {   
        // liquid water path
        try
        {
          {
            auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; //g/kg
            typename Plotter_t::arr_t snap(tmp); 
            snap += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; //g/kg
            snap *= rhod; // water per cubic metre (should be wet density...)
            res_prof(at) = blitz::mean(snap); 
          }
        }
        catch(...) {;}
      }   
      else if (plt == "surf_flux_latent")
      {   
        try
        {
          {
            typename Plotter_t::arr_t snap(plotter.h5load_timestep("latent surface flux", at * n["outfreq"], true)); 
            res_prof(at) = blitz::mean(snap); 
          }
        }
        catch(...) {;}
      }   
      else if (plt == "surf_flux_sensible")
      {   
        try
        {
          {
            typename Plotter_t::arr_t snap(plotter.h5load_timestep("sensible surface flux", at * n["outfreq"], true)); 
            res_prof(at) = blitz::mean(snap); 
          }
        }
        catch(...) {;}
      }   
      else if (plt == "er")
      {   
        //entrainment rate as in the 2009 Ackerman paper
        // to store total mixingg ratio
        try
        {
          {
            auto tmp = plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; //g/kg
            typename Plotter_t::arr_t snap(tmp); 
            snap += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; //g/kg
            rtot = snap;
          }
          {
            auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]) * 1e3;
            typename Plotter_t::arr_t snap(tmp); // vapor mixing ratio [g/kg]
            rtot += snap;
          }
          plotter.k_i = 0;
          plotter.k_i = blitz::first((rtot < 8.), plotter.LastIndex); 
          res_prof(at) = blitz::mean(plotter.k_i);
        }
        catch (...) {;}
      }
      else if (plt == "wvarmax")
      {
        // maximum variance of vertical velocity
        try
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
    //      Array<double, 1> mean(n["z"]);
          snap = snap * snap; // 2nd power, w_mean = 0
          // mean variance of w in horizontal
//          mean = blitz::mean(snap(tensor::j, tensor::i), tensor::j); // mean over x and y
          auto mean = plotter.horizontal_mean(snap);
          res_prof(at) = blitz::max(mean); // the max value
        }
        catch(...) {;}
      }
      else if (plt == "tot_tke")
      {
        try
        {
          auto u = plotter.h5load_timestep("u", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(u);
          plotter.subtract_horizontal_mean(snap);
          snap = snap * snap;
          auto mean = plotter.horizontal_mean(snap);
          res_prof(at) = blitz::sum(mean);

          {
            auto w = plotter.h5load_timestep("w", at * n["outfreq"]);
            snap = w;
            plotter.subtract_horizontal_mean(snap);
            snap = snap * snap;
            auto mean = plotter.horizontal_mean(snap);
            res_prof(at) += blitz::sum(mean);
          }
        
          if (Plotter_t::n_dims > 2)
          {
            auto v = plotter.h5load_timestep("v", at * n["outfreq"]);
            snap = v;
            plotter.subtract_horizontal_mean(snap);
            snap = snap * snap;
            auto mean = plotter.horizontal_mean(snap);
            res_prof(at) += blitz::sum(mean);
          }
          
          res_prof(at) *= 0.5 * n["dz"];

          {
            auto tke = plotter.h5load_timestep("tke", at * n["outfreq"]);
            typename Plotter_t::arr_t snap(tke);
            res_prof(at) += blitz::sum(plotter.horizontal_mean(snap));
          }
        }
        catch(...) {;}
      }
      else if (plt == "cl_gccn_conc")
      {
	// gccn (r_d > 2 um) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
          snap2 *= rhod; // b4 it was specific moment
          snap2 /= 1e6; // per cm^3
          snap2 *= snap;
          if(blitz::sum(snap) > 0)
            res_prof(at) = blitz::sum(snap2) / blitz::sum(snap); 
          else
            res_prof(at) = 0;
        }
        catch(...){;}
      }
      else if (plt == "cl_gccn_meanr")
      {
	// gccn (r_d > 2 um) mean radius in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap_m0(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
          snap *= ispositive(snap_m0); // is cloudy and has gccn
          typename Plotter_t::arr_t snap_m1(plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]));
          snap_m1 = where(snap > 0, 1e6 * snap_m1 / snap_m0, 0); // in microns
          //snap_m1 *= snap;
          if(blitz::sum(snap) > 0)
            res_prof(at) = blitz::sum(snap_m1) / blitz::sum(snap); 
          else
            res_prof(at) = 0;
        }
        catch(...){;}
      }
      else if (plt == "gccn_conc")
      {
	// gccn (r_d > 2 um) concentration
        try
        {
          typename Plotter_t::arr_t snap2(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
          snap2 /= 1e6; // per cm^3
          snap2 *= rhod; // b4 it was per milligram
          res_prof(at) = blitz::mean(snap2); 
        }
        catch(...){;}
      }
      else if (plt == "cl_meanr")
      {
	// cloud droplets mean radius in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap); // cloudiness mask
          typename Plotter_t::arr_t snap_m0(plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]));
          typename Plotter_t::arr_t snap_m1(plotter.h5load_timestep("cloud_rw_mom1", at * n["outfreq"]));
          snap_m1 = where(snap > 0, 1e6 * snap_m1 / snap_m0, 0); // in microns
          //snap_m1 *= snap;
          if(blitz::sum(snap) > 0)
            res_prof(at) = blitz::sum(snap_m1) / blitz::sum(snap); 
          else
            res_prof(at) = 0;
        }
        catch(...){;}
      }
      else assert(false);
    } // time loop

    gp << "set yrange[*:*]\n";
    gp << "set xrange[*:*]\n";

    if (plt == "clfrac")
    {
    //  res_pos *= 60.;
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
      gp << "set title 'cloud fraction'\n";
    }
    else if (plt == "ract_com")
    {
      res_prof /= 1000.;
      res_pos *= 60.;
      gp << "set ylabel 'q_c center of mass [km]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'center of mass of cloud water'\n";
    }
    else if (plt == "th_com")
    {
      res_prof /= 1000.;
      res_pos *= 60.;
      gp << "set ylabel 'th prtrb center of mass [km]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'center of mass height'\n";
    }
    else if (plt == "ract_avg")
    {
      plot_std_dev = true;
      res_pos *= 60.;
      gp << "set ylabel '<q_c>  [g/kg]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'avg cloud water mixing ratio'\n";
    }
    else if (plt == "ract_std_dev")
    {
      res_pos *= 60.;
      gp << "set ylabel 'sigma(q_c) / <q_c>'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'relative std dev of q_c'\n";
    }
    else if (plt == "cloud_avg_act_conc")
    {
      plot_std_dev = true;
      res_pos *= 60.;
      gp << "set ylabel '<N_c>  [1/cm^3]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'avg concentration of activated droplets'\n";
    }
    else if (plt == "cloud_std_dev_act_conc")
    {
      res_pos *= 60.;
      gp << "set ylabel 'sigma(N_c) / <N_c>'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'relative std dev N_c'\n";
    }
    else if (plt == "cloud_avg_supersat")
    {
      plot_std_dev = true;
      res_pos *= 60.;
      gp << "set ylabel '<S> [%]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'avg supersaturation'\n";
    }
    else if (plt == "cloud_std_dev_supersat")
    {
      res_pos *= 60.;
      gp << "set ylabel 'sigma(S) / <S>'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'relative std dev S'\n";
    }
    else if (plt == "cloud_avg_act_rad")
    {
      res_pos *= 60.;
      gp << "set ylabel '<r_{mean}> [um]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'average mean radius of activated droplets'\n";
    }
    else if (plt == "cloud_avg_std_dev_act_rad")
    {
      res_pos *= 60.;
 ///     res_prof *= 1e6;
      gp << "set ylabel '<sigma(r)> [um]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'average std dev of radius of activated droplets'\n";
    }
    else if (plt == "sd_conc_avg")
    {
      plot_std_dev = true;
      res_pos *= 60.;
      gp << "set ylabel '<N_{SD}>'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'average SD number'\n";
    }
/*
    else if (plt == "sd_conc_std_dev")
    {
      res_pos *= 60.;
      gp << "set ylabel 'sigma(N_{SD}) / <N_{SD}>'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'relative std dev of N_{SD}'\n";
    }
*/
    else if (plt == "sd_conc_act_avg")
    {
      plot_std_dev = true;
      res_pos *= 60.;
      gp << "set ylabel '<N_{SD}^{act}>'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'average activated SD number'\n";
    }
    else if (plt == "tot_water")
    {
      gp << "set title 'mean water mass density'\n";
      res_pos *= 60.;
      res_prof *= 1e3;
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel 'mean (hydrometeors+vapor) water mass density [g/m^3]'\n";
    }
    else if (plt == "com_vel")
    {
      gp << "set title 'vertical velocity of the COM\n";
      res_pos *= 60.;
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel 'w [m/s]'\n";
    }
    else if (plt == "com_supersat")
    {
      gp << "set title 'supersaturation at the COM\n";
      res_pos *= 60.;
      res_prof *= 100.; // to get %
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel 'S [%]'\n";
    }
    else if (plt == "com_mom0")
    {
      gp << "set title 'cloud drops concentration at the center of mass\n";
      res_pos *= 60.;
      res_prof /= 1e6;
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel 'N_c [cm^{-3}]'\n";
    }
    else if (plt == "com_mom1")
    {
      gp << "set title 'mean wet radius at the center of mass\n";
      res_pos *= 60.;
      res_prof *= 1e6;
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel '<r> [micron]'\n";
    }
    else if (plt == "com_mom2")
    {
      gp << "set title 'std dev of r at the center of mass\n";
      res_pos *= 60.;
      res_prof *= 1e6;
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel 'std dev [micron]'\n";
    }
    else if (plt == "com_sd_conc")
    {
      gp << "set title 'number of SDs at the center of mass\n";
      res_pos *= 60.;
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel 'SD no'\n";
    }
    else if (plt == "nc")
      gp << "set title 'average cloud drop conc [1/cm^3]'\n";
    else if (plt == "cl_nc")
    {
      gp << "set title 'average cloud drop conc [1/cm^3] in cloudy cells'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "wvarmax")
    {
      gp << "set title 'max variance of w [m^2 / s^2]'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "surf_precip")
    {
      gp << "set title 'surface precipitation [mm/d]'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "acc_precip")
    {
      gp << "set title 'accumulated surface precipitation [mm]'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "mass_dry")
      gp << "set title 'total dry mass [g]'\n";
    else if (plt == "RH_max")
    {
      gp << "set title 'max RH in the domain'\n";
      res_pos *= 60.;
      gp << "set xlabel 'time [min]'\n";
      gp << "set ylabel 'RH'\n";
    }
    else if (plt == "lwp")
    {
      gp << "set title 'liquid water path [g / m^2]'\n";
      res_prof *= (n["z"] - 1) * n["dz"]; // top and bottom cells are smaller
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "cloud_base")
    {
      gp << "set title 'cloud base [m]'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "cl_gccn_conc")
    {
      gp << "set title 'average gccn conc [1/cm^3] in cloudy cells'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "cl_gccn_meanr")
    {
      gp << "set title 'average wet radius [um] of GCCNs in cloudy cells'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "cl_meanr")
    {
      gp << "set title 'average wet radius [um] of cloud droplets in cloudy cells'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "cl_avg_cloud_rad")
    {
      gp << "set title 'average wet radius [um] of cloud droplets in cloudy cells'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "gccn_conc")
    {
      gp << "set title 'average gccn conc [1/cm^3]'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "surf_flux_sensible")
    {
      gp << "set title 'sensible surf flux [W/m^2]'\n";
    }
    else if (plt == "surf_flux_latent")
    {
      gp << "set title 'latent surf flux [W/m^2]'\n";
    }
    else if (plt == "er")
    {
      // central difference, in cm
      Range nofirstlast = Range(1, last_timestep-1);
      auto res_prof_tmp = res_prof.copy();
      res_prof(nofirstlast) = where(res_prof_tmp(nofirstlast+1) > 0., (res_prof_tmp(nofirstlast+1) - res_prof_tmp(nofirstlast-1)) * n["dz"] * 1e2 / (2 * n["dt"] * n["outfreq"])  + D * (res_prof_tmp(nofirstlast) - 0.5) * n["dz"] * 1e2, 0.);

      // larger stencil
//      Range notwo = Range(2, last_timestep-2);
   //   res_prof(notwo) = where(res_prof_tmp(notwo+1) > 0., ( 2. / 3. * (res_prof_tmp(notwo+1) - res_prof_tmp(notwo-1)) + 1. / 12. * (res_prof_tmp(notwo+2) - res_prof_tmp(notwo-2)) ) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"])  + D * (res_prof_tmp(notwo) - 0.5) * n["dz"] * 1e2, 0.);

      //res_prof(0) = 0.;
      res_prof(0) = (res_prof_tmp(1) - res_prof_tmp(0)) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"])  + D * (res_prof_tmp(0) - 0.5) * n["dz"] * 1e2;
      res_prof(last_timestep) = (res_prof_tmp(last_timestep) - res_prof_tmp(last_timestep-1)) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"])  + D * (res_prof_tmp(last_timestep) - 0.5) * n["dz"] * 1e2;
      gp << "set title 'entrainment rate [cm / s]'\n";
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
    }
    else if (plt == "tot_tke")
    {
      gp << "set xlabel ''\n";
      gp << "set ylabel ''\n";
      gp << "set title 'turbulent kinetic energy (resolved + sgs) [m^3 / s^2]'\n";
    }

    gp << "plot '-' with l";
    if(plot_std_dev)
      gp << ", '-' w l, '-' w l";
    gp << " \n";

    std::cout << plt << " " << res_pos << res_prof << res_prof_std_dev << std::endl;
    gp.send1d(boost::make_tuple(res_pos, res_prof));
    oprof_file << res_prof ;
    if(plot_std_dev)
    {
      oprof_file << res_prof_std_dev ;
      res_prof = res_prof + res_prof_std_dev;
      gp.send1d(boost::make_tuple(res_pos, res_prof));
      res_prof = res_prof - 2*res_prof_std_dev;
      gp.send1d(boost::make_tuple(res_pos, res_prof));
    }
   // plot(gp, res);
  } // var loop
}

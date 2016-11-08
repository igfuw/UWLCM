#include "common.hpp"
#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>

const double D = 3.75e-6; //[1/s], ugly, large-scale horizontal wind divergence TODO: read from model output

const std::set<std::string> plots({"wvarmax", "clfrac", "lwp", "er", "surf_precip"/*, "mass_dry"*/, "acc_precip", "cl_nc", "rc_com", "rc_avg", "th_com", "tot_water"});

template<class Plotter_t>
void plot_series(Plotter_t plotter)
{
  auto& n = plotter.map;
  for(auto elem : n)
  {
     std::cout << elem.first << " " << elem.second << std::endl;
  }
  Gnuplot gp;
  string file = plotter.file + "_series.svg";
  int hor = min<int>(plots.size(), 4);
  int ver = double(plots.size()) / 4. + 0.99999;
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

  Array<float, 1> res_prof(last_timestep - first_timestep + 1);
  Array<float, 1> res_pos(last_timestep - first_timestep + 1);

  for (auto &plt : plots)
  {
    res_prof = 0;
    res_pos = 0;

    std::ifstream f_precip(plotter.file + "/prec_vol.dat");
    std::string row;
    double prec_vol;

    for (int at = first_timestep; at <= last_timestep; ++at) // TODO: mark what time does it actually mean!
    {
      res_pos(at) = at * n["outfreq"] * n["dt"] / 3600.;
      // read in precipitation volume
      std::getline(f_precip, row);
      sscanf(row.c_str(), "%*d %lf", &prec_vol);
      if (plt == "clfrac")
      {
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng000_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap);
          res_prof(at) = blitz::mean(snap); 
        }
        catch(...){;}
      }
      // rc averaged over cells with rc > 1.e-5
      else if (plt == "rc_avg")
      {
        try
        {
          // read rc to res_tmp 
          auto tmp = plotter.h5load_rc_timestep(plotter.file, at * n["outfreq"]);

          typename Plotter_t::arr_t snap(tmp);
          
          res_tmp = iscloudy_rc(snap); // find cells with rc>1e-5
          snap *= res_tmp; // apply filter
          
          // mean only over updraught cells
          if(blitz::sum(res_tmp) > 0.)
            res_prof(at) = blitz::sum(snap) / blitz::sum(res_tmp); 
          else
            res_prof(at) = 0.;
        }
        catch(...) {;}
      }
      else if (plt == "tot_water")
      {
        try
        {
          {
            auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng002_mom3", at * n["outfreq"]) * 4./3. *3.14159265359 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_prof(at) = blitz::sum(snap);
          } 
          {
            auto tmp = plotter.h5load_timestep(plotter.file, "rv", at * n["outfreq"]);
            typename Plotter_t::arr_t snap(tmp);
            snap *= rhod;
            res_prof(at) += blitz::sum(snap);
          } 
        }
        catch(...) {;}
      }
      else if (plt == "rc_com")
      {
	// center of mass of cloud droplets
        try
        {
          auto tmp = plotter.h5load_rc_timestep(plotter.file, at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          typename Plotter_t::arr_t snap2(tmp);
          
          snap2 = snap2 * plotter.LastIndex * n["dz"];
          if(blitz::sum(snap) > 1e-3)
            res_prof(at) = blitz::sum(snap2) / blitz::sum(snap); 
          else 
            res_prof(at) = 0.;
        }
        catch(...) {;}
      }
      else if (plt == "th_com")
      {
	// center of mass of temp perturb
        try
        {
          auto tmp = plotter.h5load_timestep(plotter.file, "th", at * n["outfreq"]);
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
          auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng000_mom0", at * n["outfreq"]);
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
          auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng000_mom0", at * n["outfreq"]);
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
      else if (plt == "mass_dry")
      {
	// total dry mass
        double rho_dry = 1769; //[kg/m^3] - density of ammonium sulfate from wikipedia
        try
        {
          auto tmp = plotter.h5load_timestep(plotter.file, "rd_rng000_mom3", at * n["outfreq"]) * 4./3. * 3.14 * rho_dry * 1e3;
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
          res_prof(at) = prec_vol / (double(n["dx"]) * rhod.extent(0)) / (double(n["outfreq"]) * n["dt"] / 3600. / 24.) * 1e3;
        }
        catch(...) {;}
      }
      else if (plt == "acc_precip")
      {
        // accumulated surface precipitation [mm]
        try
        {
          if(at==0)
            res_prof(at) = prec_vol / plotter.DomainSurf * 1e3; 
          else
            res_prof(at) = res_prof(at-1) + prec_vol / plotter.DomainSurf * 1e3; 
        }
        catch(...) {;}
      }
      else if (plt == "lwp")
      {   
        // liquid water path
        try
        {
          {
            auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            typename Plotter_t::arr_t snap(tmp); // cloud water mixing ratio [g/kg]
            snap *= rhod; // cloud water per cubic metre (should be wet density...)
            res_prof(at) = blitz::mean(snap); 
          }
          {
            auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            typename Plotter_t::arr_t snap(tmp); // rain water mixing ratio [g/kg]
            snap *= rhod; // rain water per cubic metre (should be wet density...)
            res_prof(at) += blitz::mean(snap); 
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
            auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            typename Plotter_t::arr_t snap(tmp); // cloud water mixing ratio [g/kg]
            rtot = snap;
          }
          {
            auto tmp = plotter.h5load_timestep(plotter.file, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            typename Plotter_t::arr_t snap(tmp); // rain water mixing ratio [g/kg]
            rtot += snap;
          }
          {
            auto tmp = plotter.h5load_timestep(plotter.file, "rv", at * n["outfreq"]) * 1e3;
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
          auto tmp = plotter.h5load_timestep(plotter.file, "w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
    //      Array<float, 1> mean(n["z"]);
          snap = snap * snap; // 2nd power, w_mean = 0
          // mean variance of w in horizontal
//          mean = blitz::mean(snap(tensor::j, tensor::i), tensor::j); // mean over x and y
          auto mean = plotter.horizontal_mean(snap);
          res_prof(at) = blitz::max(mean); // the max value
        }
        catch(...) {;}
      }
      else assert(false);
    } // time loop

    gp << "set yrange[*:*]\n";
    gp << "set xrange[*:*]\n";

    if (plt == "clfrac")
      gp << "set title 'average cloud fraction'\n";
    else if (plt == "rc_com")
    {
      res_prof /= 1000.;
      res_pos *= 60.;
      gp << "set ylabel 'r_c center of mass [km]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'center of mass'\n";
    }
    else if (plt == "th_com")
    {
      res_prof /= 1000.;
      res_pos *= 60.;
      gp << "set ylabel 'th prtrb center of mass [km]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'center of mass'\n";
    }
    else if (plt == "rc_avg")
    {
      res_prof *= 1000.;
      res_pos *= 60.;
      gp << "set ylabel 'average r_c  [g/kg]'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'average rc'\n";
    }
    else if (plt == "tot_water")
      gp << "set title 'total water'\n";
    else if (plt == "nc")
      gp << "set title 'average cloud drop conc [1/cm^3]'\n";
    else if (plt == "cl_nc")
      gp << "set title 'average cloud drop conc [1/cm^3] in cloudy cells (should be ca. 55)'\n";
    else if (plt == "wvarmax")
      gp << "set title 'max variance of w [m^2 / s^2]'\n";
    else if (plt == "surf_precip")
      gp << "set title 'surface precipitation [mm/d]'\n";
    else if (plt == "acc_precip")
      gp << "set title 'accumulated surface precipitation [mm]'\n";
    else if (plt == "mass_dry")
      gp << "set title 'total dry mass [g]'\n";
    else if (plt == "lwp")
    {
      gp << "set title 'liquid water path [g / m^2]'\n";
      res_prof *= (n["dz"] - 1) * n["z"]; // top and bottom cells are smaller
    }
    else if (plt == "er")
    {
      // forward difference, in cm
      Range nolast = Range(0, last_timestep-1);
      res_prof(nolast) = (res_prof(nolast+1) - res_prof(nolast)) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"]) + D * (res_prof(nolast) - 0.5) * n["dz"] * 1e2;
      res_prof(last_timestep) = 0.;
      gp << "set title 'entrainment rate [cm / s]'\n";
    }

    gp << "plot '-' with line\n";
    std::cout << plt << " " << res_pos << res_prof << std::endl;
    gp.send1d(boost::make_tuple(res_pos, res_prof));

    oprof_file << res_prof ;
   // plot(gp, res);
  } // var loop
}

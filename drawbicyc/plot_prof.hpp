#include "common.hpp"
#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include <libcloudph++/common/const_cp.hpp>
#include "plots.hpp"

template<class Plotter_t>
void plot_profiles(Plotter_t plotter, Plots plots)
{

  // read opts
  po::options_description opts("profile plotting options");
  opts.add_options()
    ("prof_start", po::value<int>()->required() , "time in sec when we start collecting profiles")
    ("prof_end", po::value<int>()->required() , "time in sec when we end collecting profiles")
  ;
  po::variables_map vm;
  handle_opts(opts, vm);
  std::string prof_start_s = std::to_string(vm["prof_start"].as<int>());
  std::string prof_end_s = std::to_string(vm["prof_end"].as<int>());

  auto& n = plotter.map;
  for(auto elem : n)
  {
     std::cout << elem.first << " " << elem.second << std::endl;
  }
  Gnuplot gp; 
  string file = plotter.file + "_profiles_" + prof_start_s + "_" + prof_end_s+".svg";
  init_prof(gp, file, 3, 5); 

  string prof_file = plotter.file + "_profiles_" + prof_start_s + "_" + prof_end_s+".dat";
  std::ofstream oprof_file(prof_file);

  // read in density
  typename Plotter_t::arr_t rhod(plotter.h5load(plotter.file + "/const.h5", "G"));
  typename Plotter_t::arr_t rtot(rhod.shape());

  int k_i = 0; // inversion cell

  int first_timestep =  vm["prof_start"].as<int>() / n["dt"] / n["outfreq"];
  int last_timestep =  vm["prof_end"].as<int>() / n["dt"] / n["outfreq"];

  // some ugly constants
  const double p_1000 = 1000.;
  const double L = 2.5e6;
  const double R_d = 287.;
  const double c_p = 1004;
  const double c_pd = c_p;

  double z_i;

  for (auto &plt : plots.profs)
  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    typename Plotter_t::arr_t res(rhod.shape());
    typename Plotter_t::arr_t res_tmp(rhod.shape());
    typename Plotter_t::arr_t res_tmp2(rhod.shape());
    blitz::Array<float, 1> res_prof(n["z"]);
    blitz::Array<float, 1> res_prof2(n["z"]);
    blitz::Array<float, 1> res_pos(n["z"]);
    blitz::Range all = blitz::Range::all();
    res = 0;
    res_prof = 0;
    res_prof2 = 0;

    for (int at = first_timestep; at <= last_timestep; ++at) // TODO: mark what time does it actually mean!
    {
      std::cout << at * n["outfreq"] << std::endl;
      if (plt == "rliq")
      {
	// liquid water content (cloud + rain, missing droplets with r<0.5um!)
        //res += plotter.h5load_timestep("aerosol_rw_mom3", at * n["outfreq"]) * 4./3 * 3.1416 * 1e3 * 1e3;
        res += plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]) * 4./3 * 3.1416 * 1e3 * 1e3;
        res += plotter.h5load_timestep("rain_rw_mom3", at * n["outfreq"]) * 4./3 * 3.1416 * 1e3 * 1e3;
        gp << "set title 'liquid water [g/kg]'\n";
      }
      if (plt == "gccn_rw")
      {
	// gccn (rd>1um) droplets dry radius
        res_tmp = plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]) * 1e6;
        typename Plotter_t::arr_t snap(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
        res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        res += res_tmp;
        gp << "set title 'gccn-based droplets mean wet radius'\n";
      }
      if (plt == "non_gccn_rw")
      {
	//non gccn (rd<1um) droplets dry radius
        res_tmp = plotter.h5load_timestep("non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
        typename Plotter_t::arr_t snap(plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]));
        res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        res += res_tmp;
        gp << "set title 'non-gccn-based droplets mean wet radius'\n";
      }
      if (plt == "non_gccn_rw_down")
      {
	// non-gccn (rd<2um) droplets dry radius in downdraughts
        res_tmp = plotter.h5load_timestep("non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
        {
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]));
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp); // mean radius
        }
        {
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("w", at * n["outfreq"]));
          res_tmp2 = isdowndraught(snap); // downdraft mask
          res_tmp *= res_tmp2; // apply the mask
        }
        // mean only over downdraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);
        gp << "set title 'non-gccn-based droplets mean wet radius (downdraughts only)'\n";
      }
      if (plt == "gccn_rw_down")
      {
	// gccn (rd>2um) droplets dry radius in downdraughts
        {
          auto tmp = plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isdowndraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over downdraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);
        gp << "set title 'gccn-based droplets mean wet radius (downdraughts only)'\n";
      }
      if (plt == "non_gccn_rw_up")
      {
	// non-gccn (rd<2um) droplets dry radius in updraughts
        {
          auto tmp = plotter.h5load_timestep("non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);
        gp << "set title 'non-gccn-based droplets mean wet radius (updraughts only)'\n";
      }
      if (plt == "gccn_rw_up")
      {
	// gccn (rd>2um) droplets dry radius in updraughts
        {
          auto tmp = plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);
        gp << "set title 'gccn-based droplets mean wet radius (updraughts only)'\n";
      }
      if (plt == "ugccn_rw_down")
      {
	// ultra-gccn (rd>5um) droplets dry radius in downdraughts
        {
          auto tmp = plotter.h5load_timestep("ugccn_rw_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("ugccn_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isdowndraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over downdraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);
        gp << "set title 'ultra-gccn-based droplets mean wet radius (downdraughts only)'\n";
      }
      if (plt == "act_conc")
      {
        // 0th-mom (concentration) of droplets with RH > Sc
        {
          auto tmp = plotter.h5load_timestep("actRH_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
        }
        // updraft only
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);

        // 0th-mom of droplets with rw>rc
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
        }
        // updraft only
        res_tmp *= res_tmp2;
        res_prof2 += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);
        gp << "set title 'activated droplets concentation [1/cm^3] (updrafts)'\n";
      }
      if (plt == "act_rd")
      {
	// RH > Sc droplets first dry mom
        {
          auto tmp = plotter.h5load_timestep("actRH_rd_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        // divide by 0th-mom (number of droplets)
        {
          auto tmp = plotter.h5load_timestep("actRH_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(snap > 0 , res_tmp / snap, res_tmp);
        }
        // updraft only
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);

	// rw > rc droplets first dry mom
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        // divide by 0th-mom (number of droplets)
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(snap > 0 , res_tmp / snap, res_tmp);
        }
        // updraft only
        res_tmp *= res_tmp2;
        res_prof2 += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);
        gp << "set title 'activated droplets mean dry radius (updrafts)'\n";
      }
      if (plt == "actRH_rd")
      {
	// activated droplets dry radius
        {
          auto tmp = plotter.h5load_timestep("actRH_rd_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("actRH_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(snap > 0 , res_tmp / snap, res_tmp);
          res_tmp = snap;
        }
        res += res_tmp;
        gp << "set title 'activated (RH>Sc) droplets mean dry radius'\n";
      }
      if (plt == "actrw_rd")
      {
	// activated droplets dry radius
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = where(snap > 0 , res_tmp / snap, res_tmp);
          res_tmp = snap;
        }
        res += res_tmp;
        gp << "set title 'activated (rw>rc) droplets mean dry radius'\n";
      }
      else if (plt == "rv")
      {
        {
          auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]) * 1e3;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
        }
        res += res_tmp;
        gp << "set title 'rv [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
        gp << "set yrange [0.:0.6]\n";
        gp << "set xrange [9.2:9.5]\n";
      }
      else if (plt == "sat_RH")
      {
        {
          auto tmp = plotter.h5load_timestep("RH", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap -1;
        }
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof += where(res_pos > 0 , plotter.horizontal_sum(res_tmp) / res_pos, 0);

        res += res_tmp;
        gp << "set title 'supersaturation RH-based in updrafts only'\n";
        gp << "set yrange [0.45:1.]\n";
        gp << "set xrange [0.000:*]\n";
      }
      else if (plt == "00rtot")
      {
	// total water content (vapor + cloud + rain, missing droplets with r<0.5um!)
	 /*
        {
          auto tmp = plotter.h5load_timestep("aerosol_rw_mom3", at * n["outfreq"]) * 4./3 * 3.1416 * 1e3 * 1e3;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }*/
        {
          auto tmp = plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]) * 4./3 * 3.1416 * 1e3 * 1e3;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("rain_rw_mom3", at * n["outfreq"]) * 4./3 * 3.1416 * 1e3 * 1e3;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp += snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]) * 1e3;
          typename Plotter_t::arr_t snap(tmp);
          res_tmp += snap;
        }
        res += res_tmp;
        res_prof = plotter.horizontal_mean(res_tmp); // average in x
        // find instantaneous inversion height
        k_i +=  blitz::first((res_prof < 8.));
        gp << "set title 'total water [g/kg]'\n";
      }
      else if (plt == "N_c")
      {
	// cloud drops concentration [1/cm^3]
        {
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          res += snap; 
        }
        gp << "set title 'cloud droplets ( 0.5um < r < 25um) concentration [1/cm^3]'\n";
      }
      else if (plt == "thl")
      {
	// liquid potential temp [K]
        {
          {
            auto tmp = plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            res_tmp2 = snap; 
          }
          {
            auto tmp = plotter.h5load_timestep("rain_rw_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3;
            typename Plotter_t::arr_t snap(tmp);
            res_tmp2 += snap; 
          }
          // res_tmp2 is now q_l (liq water content)
          auto tmp = plotter.h5load_timestep("th", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp); // snap is theta_dry
          res_tmp = pow(snap * pow(rhod * R_d / (p_1000 * 100), R_d / c_pd), c_pd / (c_pd - R_d)); // res_tmp is now temperature; 1 bar = 100 000Pa
          snap *= (res_tmp - res_tmp2 * L / c_p) / res_tmp; 
          res += snap; 
//          res += res_tmp2;
        }
        gp << "set title 'liquid potential temp [K]'\n";
      }
      else if (plt == "clfrac")
      {
	// cloud fraction (cloudy if N_c > 20/cm^3)
        {
          auto tmp = plotter.h5load_timestep("cloud_rw_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap);
          res += snap; 
        }
        gp << "set title 'cloud fraction'\n";
      }
      else if (plt == "prflux")
      {
	// precipitation flux(doesnt include vertical volicty w!)
        { 
          auto tmp = plotter.h5load_timestep("precip_rate", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap = snap *  4./3 * 3.14 * 1e3 // to get mass
                     / plotter.CellVol    // averaged over cell volume, TODO: make precip rate return specific moment? wouldnt need the dx and dy
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
	// add vertical velocity to precipitation flux (3rd mom of cloud drops * w)
/*
        { 
          auto tmp = plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]); // this time its a specific moment
          typename Plotter_t::arr_t snap(tmp);
	  auto tmp2 = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap2(tmp2);
          snap = - (snap * snap2) *  4./3 * 3.14 * 1e3 // to get mass
                     * rhod           // dry air density
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
	// add vertical velocity to precipitation flux (3rd mom of rain drops * w)
        { 
          auto tmp = plotter.h5load_timestep("rain_rw_mom3", at * n["outfreq"]); // this time its a specific moment
          typename Plotter_t::arr_t snap(tmp);
	  auto tmp2 = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap2(tmp2);
          snap = - (snap * snap2) *  4./3 * 3.14 * 1e3 // to get mass
                     * rhod           // dry air density
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
*/
        // turn 3rd mom * velocity into flux in [W/m^2]
        gp << "set title 'precipitation flux [W/m^2]'\n";
      }
      else if (plt == "wvar")
      {
	// variance of vertical velocity, w_mean=0
	auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
        typename Plotter_t::arr_t snap(tmp);
        snap = snap * snap; // 2nd power
        res += snap;
        gp << "set title 'variance of w [m^2 / s^2]'\n";
      }
      else if (plt == "w3rd")
      {
	// 3rd mom of vertical velocity, w_mean=0
	auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
        typename Plotter_t::arr_t snap(tmp);
        snap = snap * snap * snap; // 3rd power
        res += snap;
        gp << "set title '3rd mom of w [m^3 / s^3]'\n";
      }
//      else assert(false);
    } // time loop
    res /= last_timestep - first_timestep + 1;
    
    z_i = (double(k_i)-0.5) / (last_timestep - first_timestep + 1) * n["dz"];
    std::cout << "average inversion height " << z_i;
    res_pos = i * n["dz"] / z_i; 

    if (plt != "act_rd" && plt != "act_conc")
    {
      if (plt == "ugccn_rw_down" || plt == "sat_RH" || plt=="gccn_rw_down" || plt=="non_gccn_rw_down" || plt=="gccn_rw_up" || plt=="non_gccn_rw_up")
        res_prof /= last_timestep - first_timestep + 1;
      else
        res_prof = plotter.horizontal_mean(res); // average in x


      // res_prof(0) is ground level - we dont know what is there? surf fluxes shouldnt be added to it?! anyway, set res_prof(0)=res_prof(1) for plotting purposes
      res_prof(0) = res_prof(1);

      gp << "plot '-' with line\n";
      gp.send1d(boost::make_tuple(res_prof, res_pos));

      oprof_file << res_prof ;

      if(plt == "rv" || plt == "sat" || plt == "sat_RH")
      {
        gp << "set yrange [0.:1.2]\n";
        gp << "set xrange [*:*]\n";
      }
    }
    else 
    {
      gp << "plot '-' with line title 'RH > Sc', '-' w l t 'rw > rc'\n";
      res_prof /= last_timestep - first_timestep + 1;
      res_prof2 /= last_timestep - first_timestep + 1;
      // res_prof(0) is ground level - we dont know what is there? surf fluxes shouldnt be added to it?! anyway, set res_prof(0)=res_prof(1) for plotting purposes
      res_prof(0) = res_prof(1);
      res_prof2(0) = res_prof2(1);
      gp.send1d(boost::make_tuple(res_prof, res_pos));
      gp.send1d(boost::make_tuple(res_prof2, res_pos));
      oprof_file << res_prof ;
      oprof_file << res_prof2 ;
    }

//    plot(gp, res);
  } // var loop
  oprof_file << z_i << std::endl;
}


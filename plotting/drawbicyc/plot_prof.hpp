#include "common.hpp"
#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include <libcloudph++/common/const_cp.hpp>
#include "plots.hpp"
#include "gnuplot_profs_set_labels.hpp"

template<class Plotter_t>
void plot_profiles(Plotter_t plotter, Plots plots, std::string type, const bool normalize)
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
  string file = plotter.file + "_" + type + "_profiles_" + prof_start_s + "_" + prof_end_s+".svg";
  init_prof(gp, file, 3, 5); 

  string prof_file = plotter.file + "_" + type + "_profiles_" + prof_start_s + "_" + prof_end_s+".dat";
  std::ofstream oprof_file(prof_file);

  // read in density
  typename Plotter_t::arr_t rhod(plotter.h5load(plotter.file + "/const.h5", "G"));
  typename Plotter_t::arr_t rtot(rhod.shape());

//  int k_i = 0; // inversion cell

  int first_timestep =  vm["prof_start"].as<int>() / int(n["dt"] * n["outfreq"]);
  int last_timestep =  vm["prof_end"].as<int>() / int(n["dt"] * n["outfreq"]);

  // some ugly constants
  const double p_1000 = 100000.;
  const double L = 2.5e6;
  const double R_d = 287.0024888;
  const double c_p = 1005;
  const double c_pd = c_p;

//  double z_i;

  blitz::firstIndex i;

  bool res_pos_out_done = false;

  blitz::Array<float, 1> res_pos(n["z"]);      // uniform vertical axis, height normalized by initial inversion height
  blitz::Array<float, 1> res_pos_hlpr(n["z"]); // actual vertical axis, height normalized by current inversion height
  if(normalize)
    res_pos = i * n["dz"] / 795; // TODO: using hardcoded DYCOMS initial inversion height just to get some uniform grid
  else
    res_pos = i * n["dz"];

  for (auto &plt : plots.profs)
  {
    blitz::secondIndex j;
    typename Plotter_t::arr_t res(rhod.shape());
    typename Plotter_t::arr_t res_tmp(rhod.shape());
    typename Plotter_t::arr_t res_tmp2(rhod.shape());
    blitz::Array<float, 1> res_prof_sum(n["z"]);  // profile interpolate to the uniform grid summed over timesteps
    blitz::Array<float, 1> res_prof(n["z"]);      // profile interpolate to the uniform grid
    blitz::Array<float, 1> res_prof_hlpr(n["z"]); // actual profile
    blitz::Array<float, 1> prof_tmp(n["z"]);
    blitz::Array<int, 1>   occur_no(n["z"]);      // number of occurances - for unusual profiles like base_prflux_vs_clhght
    blitz::Range all = blitz::Range::all();

    res_prof_sum = 0;
    occur_no = 0;

    for (int at = first_timestep; at <= last_timestep; ++at) // TODO: mark what time does it actually mean!
    {
      std::cout << at * n["outfreq"] << std::endl;
      res = 0;

      if (plt == "rliq")
      {
	// liquid water content
        res += plotter.h5load_ra_timestep(at * n["outfreq"]) * 1e3; // aerosol
        res += plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; // cloud
        res += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; // rain
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      if (plt == "gccn_rw")
      {
	// gccn (rd>1um) droplets dry radius
        res = plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]) * 1e6;
        res_tmp = plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]);
        res = where(res > 0 , res / res_tmp, res);
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      if (plt == "non_gccn_rw")
      {
	//non gccn (rd<1um) droplets dry radius
        res = plotter.h5load_timestep("non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
        res_tmp = plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]);
        res = where(res > 0 , res / res_tmp, res);
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      if (plt == "gccn_rw_cl")
      {
	// gccn (rd>2um) droplets dry radius in cloudy cells
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
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          res_tmp2 = iscloudy_rc_rico(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over downdraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
      }
      if (plt == "non_gccn_rw_cl")
      {
	// non-gccn (rd<2um) droplets dry radius in cloudy cells
        res_tmp = plotter.h5load_timestep("non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
        {
          typename Plotter_t::arr_t snap(plotter.h5load_timestep("non_gccn_rw_mom0", at * n["outfreq"]));
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp); // mean radius
        }
        {
          typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
          res_tmp2 = iscloudy_rc_rico(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over downdraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
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
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
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
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
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
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
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
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
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
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
      }
      if (plt == "act_conc_up")
      {
        // 0th-mom (concentration) of droplets with RH > Sc
        {
          // disabled, since UWLCM doesnt store actRH anymore
/*
          auto tmp = plotter.h5load_timestep("actRH_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
*/
        }
        // updraft only
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
     //     res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
       // res_prof2 += where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);

        // 0th-mom of droplets with rw>rc
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
        }
        // updraft only
        res_tmp *= res_tmp2;
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
      }
      if (plt == "nc_up")
      {
        // updraft only
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
        }
        // mean only over updraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level

        {
          auto tmp = plotter.h5load_nc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
        }
        // updraft only
        res_tmp *= res_tmp2;
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
      }
      if (plt == "cl_nc_up")
      {
        // updraft only
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
        }

        {
          auto tmp = plotter.h5load_nc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          res_tmp = snap;
          snap = iscloudy(snap); // cloudiness mask
          res_tmp2 *= snap; // cloudy updrafts only
        }

        // mean only over cloudy updraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level

        // updraft only
        res_tmp *= res_tmp2;
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
      }
      if (plt == "nc_down")
      {
        // updraft only
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isdowndraught(snap);
        }
        // mean only over updraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level

        {
          auto tmp = plotter.h5load_nc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
        }
        // updraft only
        res_tmp *= res_tmp2;
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
      }
      if (plt == "act_rd_up")
      {
	// RH > Sc droplets first dry mom
	/*
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
        */
        // updraft only
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
 //         res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
//        res_prof2 += where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);

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
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);
      }
      if (plt == "actRH_rd")
      {
	// activated droplets dry radius
        {
          auto tmp = plotter.h5load_timestep("actRH_rd_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("actRH_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res = where(snap > 0 , res / snap, res);
        }
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      if (plt == "actrw_rd")
      {
	// activated droplets dry radius
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom1", at * n["outfreq"]) * 1e6;
          typename Plotter_t::arr_t snap(tmp);
          res = snap; 
        }
        {
          auto tmp = plotter.h5load_timestep("actrw_rd_mom0", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res = where(snap > 0 , res / snap, res);
          res = snap;
        }
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "rv")
      {
        res = plotter.h5load_timestep("rv", at * n["outfreq"]) * 1e3;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "u")
      {
        res = plotter.h5load_timestep("u", at * n["outfreq"]);
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "v")
      {
        res = plotter.h5load_timestep("v", at * n["outfreq"]);
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "w")
      {
        res = plotter.h5load_timestep("w", at * n["outfreq"]);
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sd_conc")
      {
        res = plotter.h5load_timestep("sd_conc", at * n["outfreq"]);
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "vel_div")
      {
        res = plotter.h5load_timestep("vel_div", at * n["outfreq"]);
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sat_RH")
      {
        res = plotter.h5load_RH_timestep(at * n["outfreq"]);
        res = (res -1) * 100;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sat_RH_up")
      {
        {
          auto tmp = plotter.h5load_RH_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = (snap -1) * 100;
        }
 // for mean over updraught cells
        {
          auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        prof_tmp = plotter.horizontal_sum(res_tmp2); // number of downdraft cells on a given level
        res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(res_tmp) / prof_tmp, 0);

//        mean over all cells
//        res_prof += plotter.horizontal_sum(res_tmp);
      }
      else if (plt == "00rtot")
      {
        res = plotter.h5load_ra_timestep(at * n["outfreq"]) * 1e3; // aerosol
        res += plotter.h5load_rc_timestep(at * n["outfreq"]) * 1e3; // cloud
        res += plotter.h5load_rr_timestep(at * n["outfreq"]) * 1e3; // rain
        res += plotter.h5load_timestep("rv", at * n["outfreq"]) * 1e3; // vapour

        res_prof_hlpr = plotter.horizontal_mean(res); // average in x

        if(normalize)
        {
          // find instantaneous inversion height
          int k_i =  blitz::first((res_prof_hlpr < 8.)); // inversion cell index

          // normalize vertical axis
          res_pos_hlpr = i / float(k_i); 
        }

//    z_i = (double(k_i)-0.5) / (last_timestep - first_timestep + 1) * n["dz"];
//    std::cout << "average inversion height " << z_i;
      }
      else if (plt == "N_c")
      {
	// cloud drops concentration [1/cm^3]
        res = plotter.h5load_nc_timestep(at * n["outfreq"]) * rhod / 1e6; // from sepcific to normal moment + per cm^3
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "cl_nc")
      {
	// cloud droplet (0.5um < r < 25 um) concentration in cloudy grid cells
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = plotter.h5load_nc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          typename Plotter_t::arr_t snap2;
          snap2.resize(snap.shape());
          snap2=snap;
          snap = iscloudy(snap); // cloudiness mask
          snap2 *= snap;

          // mean only over cloudy cells
          prof_tmp = plotter.horizontal_sum(snap); // number of cloudy cells on a given level
          res_prof_hlpr = where(prof_tmp > 0 , plotter.horizontal_sum(snap2) / prof_tmp, 0);
        }
        catch(...){;}
      }
      else if (plt == "thl")
      {
	// liquid potential temp [K]
        {
          auto &ql(res_tmp2);
          ql  = plotter.h5load_ra_timestep(at * n["outfreq"]); // aerosol
          ql  += plotter.h5load_rc_timestep(at * n["outfreq"]); // cloud
          ql  += plotter.h5load_rr_timestep(at * n["outfreq"]); // rain
          // ql is now q_l (liq water content)
//          auto tmp = plotter.h5load_timestep("th", at * n["outfreq"]);
  //        typename Plotter_t::arr_t th_d(tmp); 
          typename Plotter_t::arr_t th(plotter.h5load_timestep("th", at * n["outfreq"]));
          ///auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]);
//          typename Plotter_t::arr_t rv(plotter.h5load_timestep("rv", at * n["outfreq"]));

          typename Plotter_t::arr_t T = th.copy();
          T *= pow(plotter.p_e(plotter.LastIndex) / p_1000, R_d / c_pd);
// (plotter.h5load_timestep("libcloud_temperature", at * n["outfreq"]));
          // init pressure, from rv just to get correct size
//          typename Plotter_t::arr_t p(rv); 
  //        T = pow(th_d * pow(rhod * R_d / (p_1000), R_d / c_pd), c_pd / (c_pd - R_d)); 
          // TODO: env pressure should be used below!
    //      p = rhod * R_d * (1 + 29./18. * rv) * T;  // Rv/Rd = 29/18
          res = th / T * (T - ql * L / c_p); 
//          res += ql;
        }
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "clfrac")
      {
	// cloud fraction (cloudy if N_c > 20/cm^3)
        {
          auto tmp = plotter.h5load_nc_timestep(at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap);
          res += snap; 
        }
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "base_prflux_vs_clhght")
      // TODO: 'normalize' messes with this plot
      {
        res_prof_hlpr = 0;
        // cloudy cells (cloudy if q_c > 0.01 g/kg as in RICO paper. NOTE: also add q_r ?)
        typename Plotter_t::arr_t snap(plotter.h5load_rc_timestep(at * n["outfreq"]));
        snap = iscloudy_rc_rico(snap); // cloudiness mask
        plotter.k_i = blitz::sum(snap, plotter.LastIndex); // sum in the vertical, assumes that all cloudy cells in a column belong to the same cloud

        plotter.tmp_int_hrzntl_slice = blitz::first(snap > 0, plotter.LastIndex); // cloud base hgt over dz

        // precipitation flux(doesnt include vertical velocity w!)
        res = plotter.h5load_prflux_timestep(at * n["outfreq"]);
        plotter.tmp_float_hrzntl_slice = 0; // just to get zero in columns without clouds
        plotter.tmp_float_hrzntl_slice = plotter.get_value_at_hgt(res, plotter.tmp_int_hrzntl_slice); // precip flux at cloud base

        // NOTE: we assume that k_i and tmp_float_hr... is contiguous in memory
        for(int i=0; i<plotter.k_i.size(); ++i)
        {
          const int cl_hgt_over_dz = *(plotter.k_i.data() + i);
          if(cl_hgt_over_dz > 0)
          {
            occur_no(cl_hgt_over_dz)+=1;
            res_prof_hlpr(cl_hgt_over_dz) += *(plotter.tmp_float_hrzntl_slice.data() + i);
          }
        }
      }
      else if (plt == "prflux")
      {
	// precipitation flux(doesnt include vertical volicty w!)
        { 
          res = plotter.h5load_prflux_timestep(at * n["outfreq"]);
          res_prof_hlpr = plotter.horizontal_mean(res); // average in x
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
      }
      else if (plt == "rad_flx")
      {
        auto tmp = plotter.h5load_timestep("radiative_flux", at * n["outfreq"]);
        typename Plotter_t::arr_t snap(tmp);
        res = snap; 
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "wvar")
      {
	// variance of vertical velocity, w_mean=0
	auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
        typename Plotter_t::arr_t snap(tmp);
        snap = snap * snap; // 2nd power
        res = snap;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "w3rd")
      {
	// 3rd mom of vertical velocity, w_mean=0
	auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
        typename Plotter_t::arr_t snap(tmp);
        snap = snap * snap * snap; // 3rd power
        res = snap;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sgs_tke")
      {
        {
          auto tmp = plotter.h5load_timestep("tke", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
        }
        res = res_tmp;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "k_m")
      {
        {
          auto tmp = plotter.h5load_timestep("k_m", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
        }
        res = res_tmp;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sgs_tht_flux")
      {
        {
          auto tmp = plotter.h5load_timestep("sgs_tht_flux", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
        }
        res = res_tmp;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sgs_rv_flux")
      {
        {
          auto tmp = plotter.h5load_timestep("sgs_rv_flux", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
        }
        res = res_tmp;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sgs_rc_flux")
      {
        {
          auto tmp = plotter.h5load_timestep("sgs_rc_flux", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
        }
        res = res_tmp;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }
      else if (plt == "sgs_u_flux")
      {
        {
          auto tmp = plotter.h5load_timestep("sgs_u_flux", at * n["outfreq"]);
          typename Plotter_t::arr_t snap(tmp);
          res_tmp = snap;
        }
        res = res_tmp;
        res_prof_hlpr = plotter.horizontal_mean(res); // average in x
      }



// ======================================================================================

      if(normalize)
      {
        // interpolate profile to uniform vertical grid (height normalized by inversion height)
        // lowest level the same, no need to interpolate
        res_prof(0) = res_prof_hlpr(0);
        for(int k = 1; k < n["z"]; ++k)
        {
          auto last = blitz::last(res_pos_hlpr < res_pos(k));
          if(last == n["z"]-1) // need to extrapolate up
          {
            res_prof(k) = res_prof_hlpr(last) + 
                         (res_prof_hlpr(last) - res_prof_hlpr(last-1)) / (res_pos_hlpr(last) - res_pos_hlpr(last-1)) *
                         (res_pos(k) - res_pos_hlpr(last));
          }
          else // interpolate
          {
            res_prof(k) = res_prof_hlpr(last) + 
                         (res_prof_hlpr(last+1) - res_prof_hlpr(last)) / (res_pos_hlpr(last+1) - res_pos_hlpr(last)) *
                         (res_pos(k) - res_pos_hlpr(last));
          }
        }
      }
      else
        res_prof = res_prof_hlpr;

      res_prof_sum += res_prof;

//      else assert(false);
    } // time loop

//    res /= last_timestep - first_timestep + 1;
    
//    z_i = (double(k_i)-0.5) / (last_timestep - first_timestep + 1) * n["dz"];
//    std::cout << "average inversion height " << z_i;
//    res_pos = i * n["dz"] / z_i; 
    if(!res_pos_out_done)
    {
      oprof_file << "position" << endl;
      oprof_file << res_pos;
      res_pos_out_done = true;
    }

    if (plt != "base_prflux_vs_clhght")
      res_prof_sum /= last_timestep - first_timestep + 1;
    else
      res_prof_sum = where(occur_no > 0, res_prof_sum/occur_no, 0);
 
    // set labels for the gnuplot plot
    gnuplot_profs_set_labels(gp, plt, normalize);

    // do the plotting
    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(res_prof_sum, res_pos));

    if (plt == "base_prflux_vs_clhght")
    {
      oprof_file << plt << " number of occurances" << endl;
      oprof_file << occur_no;
    }
    oprof_file << plt << endl;
    oprof_file << res_prof_sum ;

    // plotting two lines, eg. to compare different activation conditions (RH- vs rw- based)
/*
    {
      gp << "plot '-' with line title 'RH > Sc', '-' w l t 'rw > rc'\n";
      res_prof /= last_timestep - first_timestep + 1;
      res_prof2 /= last_timestep - first_timestep + 1;
      gp.send1d(boost::make_tuple(res_prof, res_pos));
      gp.send1d(boost::make_tuple(res_prof2, res_pos));
      oprof_file << res_prof ;
      oprof_file << res_prof2 ;
    }
*/

//    plot(gp, res);
  } // var loop
//  oprof_file << z_i << std::endl;
}


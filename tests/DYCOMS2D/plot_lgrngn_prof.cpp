#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <boost/tuple/tuple.hpp>
#include <libcloudph++/common/const_cp.hpp>

using namespace blitz;

double iscloudy(double x)
{
  return x > 20. ? 1. : 0.;
}
BZ_DECLARE_FUNCTION(iscloudy)

double isdowndraught(double x)
{
  return  x < -0.2 ? 1. : 0.;
}
BZ_DECLARE_FUNCTION(isdowndraught)

double isupdraught(double x)
{
  return  x > 0.2 ? 1. : 0.;
}
BZ_DECLARE_FUNCTION(isupdraught)

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: dir containing out_lgrngn")

  std::string
    dir = string(av[1]),
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);
  int k_i = 0; // inversion cell

  // average profile between 2h and 6h (in paper its between 2h and 6h! - do longer sims)
  int first_timestep = 3000. / n["dt"] / n["outfreq"];
  int last_timestep  = 3600. / n["dt"] / n["outfreq"];

  const double p_1000 = 1000.;
  const double L = 2.5e6;
  const double R_d = 287.;
  const double c_p = 1004;
  const double c_pd = c_p;

  double z_i;

  Gnuplot gp;
  string file = h5 + "_profiles.svg";
  init_prof(gp, file, 3, 5); 

  string prof_file = h5 + "_profiles.dat";
  std::ofstream oprof_file(prof_file);

  blitz::Array<float, 2> rhod;
  rhod.resize(n["x"], n["z"]);
  // read density
  {
    notice_macro("about to open file: " << h5)
    H5::H5File h5f(h5 + "/const.h5", H5F_ACC_RDONLY);
  
    notice_macro("about to read dataset: G")
    H5::DataSet h5d = h5f.openDataSet("G");
    H5::DataSpace h5s = h5d.getSpace();
  
    if (h5s.getSimpleExtentNdims() != 2)  
      error_macro("need 2 dimensions")
  
    hsize_t n[2];
    enum {x, z}; 
    h5s.getSimpleExtentDims(n, NULL);
  
    rhod.resize(n[x], n[z]);
  
    hsize_t 
      cnt[2] = { n[x], n[z] },  
      off[2] = { 0,    0    };  
    h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[2] = { 
      hsize_t(rhod.extent(0)), 
      hsize_t(rhod.extent(1))
    };  
    h5d.read(rhod.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(2, ext), h5s);
  }

  for (auto &plt : std::set<std::string>({"00rtot", "rliq", "thl", "wvar", "w3rd", "prflux", "act_conc", "clfrac", "N_c", "non_gccn_rw_up", "gccn_rw_up", "non_gccn_rw_down", "gccn_rw_down", "sat_RH"})) // rtot has to be first
  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::Array<float, 2> res(n["x"], n["z"]);
    blitz::Array<float, 2> res_tmp(n["x"], n["z"]);
    blitz::Array<float, 2> res_tmp2(n["x"], n["z"]);
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
        {
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res += snap; 
        }
        {
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res += snap; 
        }
        gp << "set title 'liquid water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
      }
      if (plt == "gccn_rw")
      {
	// gccn (rd>1um) droplets dry radius
        {
          auto tmp = h5load(h5, "gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "gccn_rw_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        res += res_tmp;
        gp << "set title 'gccn-based droplets mean wet radius'\n";
      }
      if (plt == "non_gccn_rw")
      {
	//non gccn (rd<1um) droplets dry radius
        {
          auto tmp = h5load(h5, "non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "non_gccn_rw_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        res += res_tmp;
        gp << "set title 'non-gccn-based droplets mean wet radius'\n";
      }
      if (plt == "non_gccn_rw_down")
      {
	// non-gccn (rd<2um) droplets dry radius in downdraughts
        {
          auto tmp = h5load(h5, "non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "non_gccn_rw_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isdowndraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over downdraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);
        gp << "set title 'non-gccn-based droplets mean wet radius (downdraughts only)'\n";
      }
      if (plt == "gccn_rw_down")
      {
	// gccn (rd>2um) droplets dry radius in downdraughts
        {
          auto tmp = h5load(h5, "gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "gccn_rw_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isdowndraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over downdraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);
        gp << "set title 'gccn-based droplets mean wet radius (downdraughts only)'\n";
      }
      if (plt == "non_gccn_rw_up")
      {
	// non-gccn (rd<2um) droplets dry radius in updraughts
        {
          auto tmp = h5load(h5, "non_gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "non_gccn_rw_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);
        gp << "set title 'non-gccn-based droplets mean wet radius (updraughts only)'\n";
      }
      if (plt == "gccn_rw_up")
      {
	// gccn (rd>2um) droplets dry radius in updraughts
        {
          auto tmp = h5load(h5, "gccn_rw_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "gccn_rw_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);
        gp << "set title 'gccn-based droplets mean wet radius (updraughts only)'\n";
      }
      if (plt == "ugccn_rw_down")
      {
	// ultra-gccn (rd>5um) droplets dry radius in downdraughts
        {
          auto tmp = h5load(h5, "ugccn_rw_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "ugccn_rw_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(res_tmp > 0 , res_tmp / snap, res_tmp);
        }
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isdowndraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over downdraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);
        gp << "set title 'ultra-gccn-based droplets mean wet radius (downdraughts only)'\n";
      }
      if (plt == "act_conc")
      {
        // 0th-mom (concentration) of droplets with RH > Sc
        {
          auto tmp = h5load(h5, "actRH_rd_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
        }
        // updraft only
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);

        // 0th-mom of droplets with rw>rc
        {
          auto tmp = h5load(h5, "actrw_rd_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap;
          res_tmp *= rhod / 1e6; // per cm^3
        }
        // updraft only
        res_tmp *= res_tmp2;
        res_prof2 += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);
        gp << "set title 'activated droplets concentation [1/cm^3] (updrafts)'\n";
      }
      if (plt == "act_rd")
      {
	// RH > Sc droplets first dry mom
        {
          auto tmp = h5load(h5, "actRH_rd_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        // divide by 0th-mom (number of droplets)
        {
          auto tmp = h5load(h5, "actRH_rd_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(snap > 0 , res_tmp / snap, res_tmp);
        }
        // updraft only
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);

	// rw > rc droplets first dry mom
        {
          auto tmp = h5load(h5, "actrw_rd_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        // divide by 0th-mom (number of droplets)
        {
          auto tmp = h5load(h5, "actrw_rd_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(snap > 0 , res_tmp / snap, res_tmp);
        }
        // updraft only
        res_tmp *= res_tmp2;
        res_prof2 += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);
        gp << "set title 'activated droplets mean dry radius (updrafts)'\n";
      }
      if (plt == "actRH_rd")
      {
	// activated droplets dry radius
        {
          auto tmp = h5load(h5, "actRH_rd_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "actRH_rd_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
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
          auto tmp = h5load(h5, "actrw_rd_mom1", at * n["outfreq"]) * 1e6;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "actrw_rd_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = where(snap > 0 , res_tmp / snap, res_tmp);
          res_tmp = snap;
        }
        res += res_tmp;
        gp << "set title 'activated (rw>rc) droplets mean dry radius'\n";
      }
      else if (plt == "rv")
      {
        {
          auto tmp = h5load(h5, "rv", at * n["outfreq"]) * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap;
        }
        res += res_tmp;
        gp << "set title 'rv [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
        gp << "set yrange [0.:0.6]\n";
        gp << "set xrange [9.2:9.5]\n";
      }
      else if (plt == "sat")
      {
        {
          auto tmp = h5load(h5, "th", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = snap;
        }
        // dry air pressure
        res_tmp = pow(res_tmp2 * rhod * R_d / pow(100000., R_d / c_pd), 1 / (1 - R_d / c_pd)); 
        // temperature
        res_tmp2 = res_tmp / (rhod * R_d);
        // saturated vapor pressure
        for(int x=0; x<n["x"]; ++x)
          for(int y=0; y<n["z"]; ++y)
            res_tmp2(x, y) = libcloudphxx::common::const_cp::p_vs<double>(res_tmp2(x, y) * si::kelvins) / si::pascals;
        
        {
          auto tmp = h5load(h5, "rv", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          // water vapor partial pressure
          res_tmp = snap * res_tmp / libcloudphxx::common::moist_air::eps<double>();
        }
        // supersaturation
        res_tmp = res_tmp / res_tmp2 -1;
        res += res_tmp;
        gp << "set title 'supersaturation'\n";
        gp << "set yrange [0.45:1.]\n";
        gp << "set xrange [0.000:0.003]\n";
      }
      else if (plt == "sat_RH")
      {
        {
          auto tmp = h5load(h5, "RH", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap -1;
        }
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          res_tmp2 = isupdraught(snap);
          res_tmp *= res_tmp2;
        }
        // mean only over updraught cells
        res_pos = blitz::sum(res_tmp2(j, i), j);
        res_prof += where(res_pos > 0 , blitz::sum(res_tmp(j, i), j) / res_pos, 0);

        res += res_tmp;
        gp << "set title 'supersaturation RH-based in updrafts only'\n";
        gp << "set yrange [0.45:1.]\n";
        gp << "set xrange [0.000:*]\n";
      }
      else if (plt == "00rtot")
      {
	// total water content (vapor + cloud + rain, missing droplets with r<0.5um!)
        {
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res_tmp += snap; 
        }
        {
          auto tmp = h5load(h5, "rv", at * n["outfreq"]) * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res_tmp += snap;
        }
        res += res_tmp;
        res_prof = blitz::mean(res_tmp(j, i), j); // average in x
        // find instantaneous inversion height
        k_i +=  blitz::first((res_prof < 8.));
        gp << "set title 'total water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
      }
      else if (plt == "N_c")
      {
	// cloud drops concentration [1/cm^3]
        {
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
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
            auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3;
            blitz::Array<float, 2> snap(tmp);
            res_tmp2 = snap; 
          }
          {
            auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3;
            blitz::Array<float, 2> snap(tmp);
            res_tmp2 += snap; 
          }
          // res_tmp2 is now q_l (liq water content)
          auto tmp = h5load(h5, "th", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp); // snap is theta_dry
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
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
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
          auto tmp = h5load(h5, "precip_rate", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          snap = snap *  4./3 * 3.14 * 1e3 // to get mass
                     / n["dx"] / n["dz"]    // averaged over cell volume, TODO: make precip rate return specific moment? wouldnt need the dx and dy
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
	// add vertical velocity to precipitation flux (3rd mom of cloud drops * w)
/*
        { 
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]); // this time its a specific moment
          blitz::Array<float, 2> snap(tmp);
	  auto tmp2 = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap2(tmp2);
          snap = - (snap * snap2) *  4./3 * 3.14 * 1e3 // to get mass
                     * rhod           // dry air density
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
	// add vertical velocity to precipitation flux (3rd mom of rain drops * w)
        { 
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]); // this time its a specific moment
          blitz::Array<float, 2> snap(tmp);
	  auto tmp2 = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap2(tmp2);
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
	auto tmp = h5load(h5, "w", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        snap = snap * snap; // 2nd power
        res += snap;
        gp << "set title 'variance of w [m^2 / s^2]'\n";
      }
      else if (plt == "w3rd")
      {
	// 3rd mom of vertical velocity, w_mean=0
	auto tmp = h5load(h5, "w", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
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
        res_prof = blitz::mean(res(j, i), j); // average in x

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
      gp.send1d(boost::make_tuple(res_prof, res_pos));
      gp.send1d(boost::make_tuple(res_prof2, res_pos));
      oprof_file << res_prof ;
      oprof_file << res_prof2 ;
    }

//    plot(gp, res);
  } // var loop
  oprof_file << z_i << std::endl;
} // main

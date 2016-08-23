#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <boost/tuple/tuple.hpp>

using namespace blitz;

double iscloudy(double x)
{
  return x > 20. ? 1. : 0.;
}
BZ_DECLARE_FUNCTION(iscloudy)

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: dir containing out_lgrngn")

  std::string
    dir = string(av[1]),
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);

  const double z_i = 795; // [m]
  const double D = 3.75e-6; //[1/s], ugly, large-scale horizontal wind divergence

  Gnuplot gp;
  string file = h5 + "_series.svg";
  init_prof(gp, file, 3, 3); 

  string prof_file = h5 + "_series.dat";
  std::ofstream oprof_file(prof_file);

  // read density
  blitz::Array<float, 2> rhod;
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
      cnt[3] = { n[x], n[z] },  
      off[3] = { 0,    0    };  
    h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[2] = { 
      hsize_t(rhod.extent(0)), 
      hsize_t(rhod.extent(1)) 
    };  
    h5d.read(rhod.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(2, ext), h5s);
  }

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::Array<float, 1> res_prof(n["t"]);
  blitz::Array<float, 1> res_pos(n["t"]);
  blitz::Array<int, 1> k_i(n["x"]); // index of the inversion cell
  blitz::Array<float, 2> rtot(n["x"],  n["z"]); 
  blitz::Range all = blitz::Range::all();

  for (auto &plt : std::set<std::string>({"wvarmax", "clfrac", "lwp", "er", "surf_precip", "mass_dry", "acc_precip", "cl_nc"}))
  {
    res_prof = 0;
    res_pos = 0;

    std::ifstream f_precip(h5 + "/prec_vol.dat");
    std::string row;
    double prec_vol;

    for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
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
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap);
          res_prof(at) = blitz::mean(snap); 
        }
        catch(...){;}
      }
      else if (plt == "nc")
      {
	// cloud droplet (0.5um < r < 25 um) concentration
        try
        {
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
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
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          blitz::Array<float, 2> snap2;
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
          auto tmp = h5load(h5, "rd_rng000_mom3", at * n["outfreq"]) * 4./3. * 3.14 * rho_dry * 1e3;
          blitz::Array<float, 2> snap(tmp);
          snap *= rhod * n["dx"] * n["dz"]; // turn mixing ratio in g/kg to total mass in g
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
            res_prof(at) = prec_vol / (double(n["dx"]) * rhod.extent(0)) * 1e3; 
          else
            res_prof(at) = res_prof(at-1) + prec_vol / (double(n["dx"]) * rhod.extent(0)) * 1e3; 
        }
        catch(...) {;}
      }
      else if (plt == "lwp")
      {   
        // liquid water path
        try
        {
          {
            auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            blitz::Array<float, 2> snap(tmp); // cloud water mixing ratio [g/kg]
            snap *= rhod; // cloud water per cubic metre (should be wet density...)
            res_prof(at) = blitz::mean(snap); 
          }
          {
            auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            blitz::Array<float, 2> snap(tmp); // rain water mixing ratio [g/kg]
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
            auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            blitz::Array<float, 2> snap(tmp); // cloud water mixing ratio [g/kg]
            rtot = snap;
          }
          {
            auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            blitz::Array<float, 2> snap(tmp); // rain water mixing ratio [g/kg]
            rtot += snap;
          }
          {
            auto tmp = h5load(h5, "rv", at * n["outfreq"]) * 1e3;
            blitz::Array<float, 2> snap(tmp); // vapor mixing ratio [g/kg]
            rtot += snap;
          }
          k_i = 0;
          k_i = blitz::first((rtot < 8.), j); 
          res_prof(at) = blitz::mean(k_i);
        }
        catch (...) {;}
      }
      else if (plt == "wvarmax")
      {
        // maximum variance of vertical velocity
        try
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          blitz::Array<float, 1> mean(n["z"]);
          snap = snap * snap; // 2nd power, w_mean = 0
          // mean variance of w in horizontal
          mean = blitz::mean(snap(j, i), j); // mean over x and y
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
      blitz::Range nolast = blitz::Range(0, n["t"]-2);
      res_prof(nolast) = (res_prof(nolast+1) - res_prof(nolast)) * n["dz"] * 1e2 / (n["dt"] * n["outfreq"]) + D * (res_prof(nolast) - 0.5) * n["dz"] * 1e2;
      res_prof(n["t"]-1) = 0.;
      gp << "set title 'entrainment rate [cm / s]'\n";
    }

    gp << "plot '-' with line\n";
    std::cout << plt << " " << res_pos << res_prof << std::endl;
    gp.send1d(boost::make_tuple(res_pos, res_prof));

    oprof_file << res_prof ;
//    plot(gp, res);
  } // var loop
} // main

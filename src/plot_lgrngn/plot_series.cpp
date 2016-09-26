#include "common.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <boost/tuple/tuple.hpp>

// helper type displayer
template<typename T>
class TD;

using namespace blitz;

double iscloudy(double x)
{
  return x > 20. ? 1. : 0.;
}
BZ_DECLARE_FUNCTION(iscloudy)

template<class Plotter_t>
void plot_series(const Plotter_t &plotter)
{
  auto& n = plotter.map;
  for(auto elem : n)
  {
     std::cout << elem.first << " " << elem.second << std::endl;
  }
}

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: dir containing out_lgrngn")

  std::string
    dir = string(av[1]),
    h5  = dir + "out_lgrngn";

  int n_dims = 2;

  if(n_dims == 2)
    plot_series(Plotter_t<2>(h5));
//  else if(n_dims == 3)
//    plot_series(Plotter_t<3>(h5));

return 0;
/*
  const double D = 3.75e-6; //[1/s], ugly, large-scale horizontal wind divergence TODO: read from model output

  Gnuplot gp;
  string file = h5 + "_series.svg";
  init_prof(gp, file, 3, 3); 

  string prof_file = h5 + "_series.dat";
  std::ofstream oprof_file(prof_file);

  // read density
  Array<float, 2> rhod(n["x"], n["z"]);
  rhod = h5load(h5 + "/const.h5", "G");

  Array<float, 1> res_prof(n["t"]);
  Array<float, 1> res_pos(n["t"]);
  Array<int, 1> k_i(n["x"]); // index of the inversion cell
  Array<float, 2> rtot(n["x"],  n["z"]); 

  std::set<std::string> plots({"wvarmax", "clfrac", "lwp", "er", "surf_precip", "mass_dry", "acc_precip", "cl_nc"});
  for (auto &plt : plots)
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
          auto tmp = h5load_timestep(h5, "rw_rng000_mom0", at * n["outfreq"]);
          Array<float, 2> snap(tmp);
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
          auto tmp = h5load_timestep(h5, "rw_rng000_mom0", at * n["outfreq"]);
          Array<float, 2> snap(tmp);
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
          auto tmp = h5load_timestep(h5, "rw_rng000_mom0", at * n["outfreq"]);
          Array<float, 2> snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          Array<float, 2> snap2;
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
          auto tmp = h5load_timestep(h5, "rd_rng000_mom3", at * n["outfreq"]) * 4./3. * 3.14 * rho_dry * 1e3;
          Array<float, 2> snap(tmp);
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
            auto tmp = h5load_timestep(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            Array<float, 2> snap(tmp); // cloud water mixing ratio [g/kg]
            snap *= rhod; // cloud water per cubic metre (should be wet density...)
            res_prof(at) = blitz::mean(snap); 
          }
          {
            auto tmp = h5load_timestep(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            Array<float, 2> snap(tmp); // rain water mixing ratio [g/kg]
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
            auto tmp = h5load_timestep(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            Array<float, 2> snap(tmp); // cloud water mixing ratio [g/kg]
            rtot = snap;
          }
          {
            auto tmp = h5load_timestep(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
            Array<float, 2> snap(tmp); // rain water mixing ratio [g/kg]
            rtot += snap;
          }
          {
            auto tmp = h5load_timestep(h5, "rv", at * n["outfreq"]) * 1e3;
            Array<float, 2> snap(tmp); // vapor mixing ratio [g/kg]
            rtot += snap;
          }
          k_i = 0;
          k_i = blitz::first((rtot < 8.), tensor::j); 
          res_prof(at) = blitz::mean(k_i);
        }
        catch (...) {;}
      }
      else if (plt == "wvarmax")
      {
        // maximum variance of vertical velocity
        try
        {
          auto tmp = h5load_timestep(h5, "w", at * n["outfreq"]);
          Array<float, 2> snap(tmp);
          Array<float, 1> mean(n["z"]);
          snap = snap * snap; // 2nd power, w_mean = 0
          // mean variance of w in horizontal
          mean = blitz::mean(snap(tensor::j, tensor::i), tensor::j); // mean over x and y
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
      Range nolast = Range(0, n["t"]-2);
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
*/
} // main

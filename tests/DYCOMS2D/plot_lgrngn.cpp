#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <unordered_set>
#include <iomanip> 

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: out_lgrngn parent dir")

  std::string
    dir = string(av[1]) + "/tests/fig_a/", 
    h5  = dir + "out_lgrngn";

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::Range all = blitz::Range::all();
  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::unordered_set<std::string>({"rl", "rr", "nc", "nr", "ef", "na", "th", "rv", "u", "w", "sd_conc", "r_dry", "RH"}))
    {
      std::cout << at * n["outfreq"] << " : " << plt << std::endl;
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

      if (plt == "rl")
      {
	// cloud water content
	//                                                         rho_w  kg2g
	auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
        std::string title = "cloud (0.5um < r < 25um) water mixing ratio [g/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
	plot(gp, tmp);
      }
      else if (plt == "rr")
      {
	// rain water content
	//                                                         rho_w  kg2g
	auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
	gp << "set logscale cb\n";
 	std::string title = "rain (r > 25um) water mixing ratio [g/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [1e-2:1]\n";
	plot(gp, tmp);
	gp << "unset logscale cb\n";
      }
      else if (plt == "nc")
      {
	// cloud particle concentration
	auto tmp = 1e-6 * h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
	std::string title ="cloud (0.5um < r < 25um) droplet spec. conc. [mg^{-1}]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [0:150]\n";
	plot(gp, tmp);
      }
      else if (plt == "r_dry")
      {
        // dry mass content
        // assume ammonium sulfate density of 1769 kg / m^3 (c.f. wikipedia)
        double rho_dry = 1769;
	auto tmp = h5load(h5, "rd_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e9 * rho_dry;
	gp << "set logscale cb\n";
	std::string title ="dry mass mixing ratio [ug/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [1e-2:1]\n";
	plot(gp, tmp);
	gp << "unset logscale cb\n";
      }
      else if (plt == "nr")
      {
	// rain particle concentration
	auto tmp = 1e-6 * h5load(h5, "rw_rng001_mom0", at * n["outfreq"]);
	std::string title = "rain (r > 25um) drop spec. conc. [mg^{-1}]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [.01:10]\n";
	gp << "set logscale cb\n";
	plot(gp, tmp);
	gp << "unset logscale cb\n";
      }
      else if (plt == "ef")
      {
	// effective radius
	auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) / h5load(h5, "rw_rng000_mom2", at * n["outfreq"]) * 1e6;
	std::string title = "cloud (0.5um < r < 25um) droplet effective radius [μm]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [1:20]\n";
	plot(gp, tmp);
      }
/*
      else if (plt == "na")
      {
	// aerosol concentration
	blitz::Array<float, 3> tmp(h5load(h5, "rw_rng002_mom0", at * n["outfreq"]));
	vector<quantity<si::length>> left_edges = bins_wet();
	for (int i = 1; i < left_edges.size()-1; ++i)
	{
	  if (left_edges[i + 1] > 1e-6 * si::metres) break;
	  ostringstream str;
	  str << "rw_rng" << std::setw(3) << std::setfill('0') << i + 2  << "_mom0";
	  tmp = tmp + h5load(h5, str.str(), at * n["outfreq"]);
	}
        blitz::Array<float, 2> mean(n["x"], n["z"]);
        mean = blitz::mean(tmp(i, k, j), k); // average over 2nd dim
//	gp << "set cbrange [" << 0 << ":" << 150 << "]\n";
	gp << "set title 'aerosol concentration [mg^{-1}]'\n";
	tmp /= 1e6;
	plot(gp, tmp);
      }
*/
      else if (plt == "rv")
      {   
	std::string title = "water vapour mixing ratio [g/kg]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "rv", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "th")
      {   
	std::string title = "dry air potential temperature [K]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "th", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "u")
      {   
	std::string title = "velocity in x [m/s]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "u", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "w")
      {   
	std::string title = "velocity in z [m/s]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "w", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "RH")
      {   
	std::string title = "relative humidity"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "RH", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "sd_conc")
      {   
	std::string title = "number of super-droplets"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "sd_conc", at * n["outfreq"]);
        plot(gp, tmp);
      }   
    } // var loop
  } // time loop
} // main

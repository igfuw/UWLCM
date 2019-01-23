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
    dir = string(av[1]),
    h5  = "out_blk_2m";

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::Range all = blitz::Range::all();
  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::unordered_set<std::string>({"rl", "rr", "nc", "nr", "th", "rv", "u", "w"}))
    {
      std::cout << at * n["outfreq"] << " : " << plt << std::endl;
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".png", 1, 1, n);

      if (plt == "rl")
      {
	// cloud water content
	auto tmp = h5load(h5, "rc", at * n["outfreq"]) * 1e3;
        std::string title = "cloud (0.5um < r < 25um) water mixing ratio [g/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
	gp << "set cbrange [0:1.2]\n";
	plot(gp, tmp);
      }
      else if (plt == "rr")
      {
	// rain water content
	auto tmp = h5load(h5, "rr", at * n["outfreq"]) * 1e3;
	gp << "set logscale cb\n";
 	std::string title = "rain (r > 25um) water mixing ratio [g/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
	gp << "set cbrange [1e-3:1]\n";
	plot(gp, tmp);
	gp << "unset logscale cb\n";
      }
      else if (plt == "nc")
      {
	// cloud particle concentration
	auto tmp = 1e-6 * h5load(h5, "nc", at * n["outfreq"]);
	std::string title ="cloud (0.5um < r < 25um) droplet spec. conc. [mg^{-1}]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [0:150]\n";
	plot(gp, tmp);
      }
      else if (plt == "nr")
      {
	// rain particle concentration
	auto tmp = 1e-6 * h5load(h5, "nr", at * n["outfreq"]);
	std::string title = "rain (r > 25um) drop spec. conc. [mg^{-1}]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
	gp << "set cbrange [.01:20]\n";
	gp << "set logscale cb\n";
	plot(gp, tmp);
	gp << "unset logscale cb\n";
      }
      else if (plt == "rv")
      {
	std::string title = "water vapour mixing ratio [g/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "rv", at * n["outfreq"]) * 1e3;
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
    } // var loop
  } // time loop
} // main

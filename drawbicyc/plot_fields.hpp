#include <unordered_set>
#include <iomanip> 

#include "common.hpp"
#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include "plots.hpp"

template<class Plotter_t>
void plot_fields(Plotter_t plotter, Plots plots, std::string type)
{
  // read opts
  po::options_description opts("fields plotting options");
  opts.add_options()
    ("field_plotfreq", po::value<int>()->required() , "interval in sec with which fields are plotted")
  ;
  po::variables_map vm;
  handle_opts(opts, vm);
  const int plotfreq = vm["field_plotfreq"].as<int>();

  auto& n = plotter.map;

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::Range all = blitz::Range::all();

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    if(int(at * n["outfreq"]) % plotfreq != 0) continue;
    for (auto &plt : plots.fields)
    {
      std::cout << at * n["outfreq"] << " : " << plt << std::endl;
      Gnuplot gp;
      init(gp, plotter.file + "_" + type + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n, 1, 0.666666); 

      if (plt == "rl")
      {
        try{
	// cloud water content
        auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]) * 1e3;

        std::string title = "cloud water mixing ratio [g/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
	plotter.plot(gp, tmp);
        }
        catch(...){}
      }
      else if (plt == "rr")
      {
        try{
	// rain water content
	//                                                         rho_w  kg2g
	auto tmp = plotter.h5load_timestep("rain_rw_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
	gp << "set logscale cb\n";
 	std::string title = "rain (r > 25um) water mixing ratio [g/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [1e-2:1]\n";
	plotter.plot(gp, tmp);
	gp << "unset logscale cb\n";
        }
        catch(...){}
      }
      else if (plt == "nc")
      {
	// cloud particle concentration
        try{
	auto tmp = plotter.h5load_timestep("actrw_rw_mom0", at * n["outfreq"]) * 1e-6;
	std::string title ="activated droplet spec. conc. [mg^{-1}]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [0:150]\n";
	plotter.plot(gp, tmp);
        }
        catch(...){}
      }
      else if (plt == "r_dry")
      {
        try{
        // dry mass content
        // assume ammonium sulfate density of 1769 kg / m^3 (c.f. wikipedia)
        double rho_dry = 1769;
	auto tmp = plotter.h5load_timestep("rd_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e9 * rho_dry;
	gp << "set logscale cb\n";
	std::string title ="dry mass mixing ratio [ug/kg]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [1e-2:1]\n";
	plotter.plot(gp, tmp);
	gp << "unset logscale cb\n";
        }
        catch(...){}
      }
      else if (plt == "nr")
      {
        try{
	// rain particle concentration
	auto tmp = plotter.h5load_timestep("rain_rw_mom0", at * n["outfreq"]) * 1e-6;
	std::string title = "rain (r > 25um) drop spec. conc. [mg^{-1}]";
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [.01:10]\n";
	gp << "set logscale cb\n";
	plotter.plot(gp, tmp);
	gp << "unset logscale cb\n";
        }
        catch(...){}
      }
      else if (plt == "ef")
      {
        try{
	// effective radius
	auto tmp = plotter.h5load_timestep("cloud_rw_mom3", at * n["outfreq"]) / plotter.h5load_timestep("cloud_rw_mom2", at * n["outfreq"]) * 1e6;
	std::string title = "cloud (0.5um < r < 25um) droplet effective radius [μm]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
//	gp << "set cbrange [1:20]\n";
	plotter.plot(gp, tmp);
        }
        catch(...){}
      }
/*
      else if (plt == "na")
      {
	// aerosol concentration
	blitz::Array<float, 3> tmp(h5load(h5, "rw_rng002_mom0", at * n["outfreq"]));
	auto left_edges = bins_wet();
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
	plotter.plot(gp, tmp);
      }
*/
      else if (plt == "mrk")
      {   
        try{
	std::string title = "marker"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("mrk", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "rv")
      {   
        try{
	std::string title = "water vapour mixing ratio [g/kg]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("rv", at * n["outfreq"]) * 1e3;
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "th")
      {   
        try{
	std::string title = "dry air potential temperature [K]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("th", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "u")
      {   
        try{
	std::string title = "velocity in x [m/s]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("u", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "w")
      {   
        try{
	std::string title = "velocity in z [m/s]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("w", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "vel_div")
      {   
        try{
	std::string title = "velocity field divergence"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("vel_div", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "RH")
      {   
        try{
	std::string title = "relative humidity"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("RH", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "lib_pres")
      {   
        try{
	std::string title = "libcloud pressure [Pa]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("libcloud_pressure", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "lib_temp")
      {   
        try{
	std::string title = "libcloud temperature [K]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("libcloud_temperature", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "supersat")
      {   
        try{
	std::string title = "supersaturation [%]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        typename Plotter_t::arr_t snap(plotter.h5load_timestep("RH", at * n["outfreq"]));
        snap -= 1.;
        snap *= 100;
        gp << "set cbrange[0:3]\n";
        plotter.plot(gp, snap);
        }
        catch(...){}
      }   
      else if (plt == "sd_conc")
      {   
        try{
	std::string title = "number of super-droplets"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = plotter.h5load_timestep("sd_conc", at * n["outfreq"]);
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "gccn_conc")
      {   
        try{
	std::string title = "gccn(rd>2um) concentration [1/kg]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        typename Plotter_t::arr_t tmp(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
        plotter.plot(gp, tmp);
        }
        catch(...){}
      }   
      else if (plt == "gccn_mean_rw")
      {   
        try{
	std::string title = "gccn(rd>2um) mean rw [um]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        //auto mom1 = plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]);
        //auto mom0 = plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]);
        //auto meanr = mom1 / mom0;
        typename Plotter_t::arr_t tmp1(plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
        typename Plotter_t::arr_t tmp2(plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"])*1e6);
        typename Plotter_t::arr_t tmp3((plotter.h5load_timestep("gccn_rw_mom1", at * n["outfreq"]) * 1e6) / plotter.h5load_timestep("gccn_rw_mom0", at * n["outfreq"]));
        auto tmp4 = tmp2 / tmp1;
        plotter.plot(gp, tmp4);
        }
        catch(...){}
      }   
    } // var loop
  } // time loop
}

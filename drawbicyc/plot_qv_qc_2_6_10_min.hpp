#include <unordered_set>
#include <iomanip> 

#include "common.hpp"
#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include "plots.hpp"

template<class Plotter_t>
void plot_qv_qc_2_6_10_min(Plotter_t plotter, Plots plots)
{
  using std::string;
  auto& n = plotter.map;

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::Range all = blitz::Range::all();

  Gnuplot gp;
  init(gp, "plots/out_lgrngn_qv_qc_2_6_10_min.tex", 2, 3, n); 
  gp << "unset multiplot\n";
  gp << "set pm3d map\n";
  gp << "set terminal epslatex color size 5, 2.8\n";
  gp << "set multiplot layout 2,3\n";

  std::set<float> times = {120, 360, 600};

  string XTICS = "set xtics scale .5 rotate by 60 ('' 0.0, '0.6' 30, '1.2' 60, '1.8' 90, '2.4' 120, '3.0' 150, '' 180); set xlabel 'x [km]'\n";
  string NOXTICS = "set xtics scale .5 ('' 0.0, '' 30, '' 60, '' 90, '' 120, '' 150, '' 180); unset xlabel\n";
  string YTICS = "set xtics scale .5 rotate by 60 ('' 0.0, '0.6' 30, '1.2' 60, '1.8' 90, '2.4' 120); set ylabel 'y [km]'\n";
  string NOYTICS = "set xtics scale .5 rotate by 60 ('' 0.0, '' 30, '' 60, '' 90, '' 120); unset ylabel\n";
  // Margins for each row resp. column
  string TMARGIN = "set tmargin at screen .9; set bmargin at screen 0.5\n";
  string BMARGIN = "set tmargin at screen 0.5; set bmargin at screen 0.1\n";
  
  string LMARGIN = "set lmargin at screen 0.1; set rmargin at screen 0.36\n";
  string LRMARGIN = "set lmargin at screen 0.36; set rmargin at screen 0.62\n";
  string RMARGIN = "set lmargin at screen 0.62; set rmargin at screen 0.88\n";
  //color palette
  string PALETTE = "set colorbox\n";
  string NOPALETTE = "unset colorbox\n";
  //titles
  string TITLEQV = "set title offset 0,-0.8  '$q_v$ [g/kg]'\n";
  string TITLEQC = "set title offset 0,-0.8 '$q_c$ [g/kg]'\n";
  string NOTITLE = "unset title\n";
  
  
  gp << "set cbrange [0:8]\n";
  gp << "set cbtics (0, 2, 4, 6, 8)\n";

  // --- GRAPH a
  gp << TMARGIN;
  gp << LMARGIN;
  gp << NOXTICS;
  gp << YTICS;
  gp << NOPALETTE;
  gp << "set title offset 0, -0.8 '$q_v$ [g/kg], t = 2 min'\n";
  try{
  typename Plotter_t::arr_t tmp(plotter.h5load_timestep(plotter.file, "rv", 120) * 1e3);
  std::cout << tmp;
  plotter.plot(gp, tmp);
  }
  catch(...){}

  // --- GRAPH b
  gp << TMARGIN;
  gp << LRMARGIN;
  gp << NOXTICS;
  gp << NOYTICS;
  gp << NOPALETTE;
  gp << "set title offset 0, -0.8 '$q_v$ [g/kg], t = 6 min'\n";
  try{
  typename Plotter_t::arr_t tmp(plotter.h5load_timestep(plotter.file, "rv", 320) * 1e3);
  std::cout << tmp;
  plotter.plot(gp, tmp);
  }
  catch(...){}

  // --- GRAPH c
  gp << TMARGIN;
  gp << RMARGIN;
  gp << NOXTICS;
  gp << NOYTICS;
  gp << PALETTE;
  gp << "set title offset 0, -0.8 '$q_v$ [g/kg], t = 10 min'\n";
  try{
  typename Plotter_t::arr_t tmp(plotter.h5load_timestep(plotter.file, "rv", 600) * 1e3);
  std::cout << tmp;
  plotter.plot(gp, tmp);
  }
  catch(...){}

  
  gp << "set cbrange [0:2]\n";
  gp << "set cbtics (0, 0.5, 1, 1.5, 2)\n";

  // --- GRAPH d
  gp << BMARGIN;
  gp << LMARGIN;
  gp << XTICS;
  gp << YTICS;
  gp << NOPALETTE;
  gp << "set title offset 0, -0.8 '$q_c$ [g/kg], t = 2 min'\n";
  try{
  // cloud water content
  typename Plotter_t::arr_t tmp(plotter.h5load_ract_timestep(plotter.file, 60) * 1e3);
  std::cout << tmp;
  plotter.plot(gp, tmp);
  }
  catch(...){}

  // --- GRAPH e
  gp << BMARGIN;
  gp << LRMARGIN;
  gp << XTICS;
  gp << NOYTICS;
  gp << NOPALETTE;
  gp << "set title offset 0, -0.8 '$q_c$ [g/kg], t = 6 min'\n";
  try{
  // cloud water content
  typename Plotter_t::arr_t tmp(plotter.h5load_ract_timestep(plotter.file, 360) * 1e3);
  std::cout << tmp;
  plotter.plot(gp, tmp);
  }
  catch(...){}

  // --- GRAPH f
  gp << BMARGIN;
  gp << RMARGIN;
  gp << XTICS;
  gp << NOYTICS;
  gp << PALETTE;
  gp << "set title offset 0, -0.8 '$q_c$ [g/kg], t = 10 min'\n";
  try{
  // cloud water content
  typename Plotter_t::arr_t tmp (plotter.h5load_ract_timestep(plotter.file, 600) * 1e3);
  std::cout << tmp;
  plotter.plot(gp, tmp);
  }
  catch(...){}

}

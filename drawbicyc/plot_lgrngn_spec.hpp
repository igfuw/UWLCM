#include "common.hpp"
#include "PlotterMicro.hpp"
#include "plots.hpp"

#include "bins.hpp"
//#include "gnuplot.hpp"
//#include "hdf5.hpp"

#include <map>

/*#include <unordered_set>
#include <iomanip> 

#include "common.hpp"
#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include "plots.hpp"
*/

//plot spectra positions
template<class Plotter_t>
void plot_lgrngn_spec_positions(Plotter_t plotter, Plots plots, int at)
{
  auto& n = plotter.map;
  Gnuplot gp;

  try{
  // cloud water content
  auto tmp = plotter.h5load_ract_timestep(at * n["outfreq"]) * 1e3;

  std::string title = "cloud water mixing ratio [g/kg]";
  gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
  init(gp, plotter.file + "_spectra_positions_at" + zeropad(at) + ".svg", 1, 1, n, 2., 1); 
  plotter.plot(gp, tmp, blitz::Range(yslice_idx, yslice_idx));
  }
  catch(...){}
}



//plot spectra
template<class Plotter_t>
void plot_lgrngn_spec(Plotter_t plotter, Plots plots, int at)
{
  auto& n = plotter.map;
  for(auto elem : n)
  {
     std::cout << elem.first << " " << elem.second << std::endl;
  }
  Gnuplot gp;
//  int off = 2; // TODO!!!
  int off = 0; // TODO!!!
  string file = plotter.file + "_spectra_at" + zeropad(at) + ".svg";

  int hor = min<int>(focus_3d.size(), 2);
  int ver = double(focus_3d.size()) / 2. + 0.99999;

  init_prof(gp, file, ver, hor);

  gp << "set xrange [0:30]\n";
  gp << "set yrange [0:50]\n";
  gp << "set xlabel 'r [μm]'\n"; 
  gp << "set ylabel 'n [cm^{-3} μm^{-1}]'\n"; 

  // read in density
  auto tmp = plotter.h5load(plotter.file + "/const.h5", "G");
  typename Plotter_t::arr_t rhod(tmp);

  // focus to the gridbox from where the size distribution is plotted
  char lbl = 'i';
  for (auto &fcs : focus_3d)
  {
    const int &x = fcs[0], &y = fcs[1], &z = fcs[2];

    //gp << "set label 1 '(" << lbl << ")' at graph -.15, 1.02 font ',20'\n";
    //gp << "set title 'x=" << x << " y=" << y << "'\n";

//    std::map<float, float> focus_d;
    std::map<float, float> focus_w;

    //info on the number and location of histogram edges
//    vector<quantity<si::length>> left_edges_rd = bins_dry();
  //    int nsd = left_edges_rd.size() - 1;
    vector<quantity<si::length>> left_edges_rw = bins_wet();
    int nsw = left_edges_rw.size() - 1;

/*
    for (int i = 0; i < nsd; ++i)
    {
      const string name = "rd_rng" + zeropad(i) + "_mom0";
      blitz::Array<float, 2> tmp_d(1e-6 * h5load(h5, name, at));

      focus_d[left_edges_rd[i] / 1e-6 / si::metres] = sum(tmp_d(
        blitz::Range(x-1, x+1),
        blitz::Range(y-1, y+1)
      )) 
      / 9  // mean over 9 gridpoints
      / ((left_edges_rd[i+1] - left_edges_rd[i]) / 1e-6 / si::metres); // per micrometre
    }
*/

    for (int i = 0; i < nsw; ++i)
    {
      const string name = "rw_rng" + zeropad(i + off) + "_mom0";
      auto tmp_w = plotter.h5load_timestep(name, at * n["outfreq"]);
//      blitz::Array<float, 2> tmp_w(1e-6 * h5load(h5, name, at));

      focus_w[left_edges_rw[i] / 1e-6 / si::metres] =
/*
 (tmp_w * rhod)(x,y,z) * 1e-6 // per cm^{-3}
        / ((left_edges_rw[i+1] - left_edges_rw[i]) / 1e-6 / si::metres); // per micrometre
*/
      sum((tmp_w*rhod)(
        blitz::Range(x-1, x+1),
        blitz::Range(y-1, y+1),
        blitz::Range(z-1, z+1)
      )) * 1e-6 // per cm^{-3}
      / 27 
      / ((left_edges_rw[i+1] - left_edges_rw[i]) / 1e-6 / si::metres); // per micrometre
    }
    const string name = "rw_rng" + zeropad(nsw + off) + "_mom0";
    auto tmp_w = plotter.h5load_timestep(name, at * n["outfreq"]);

    notice_macro("setting-up plot parameters");
    gp << "set title 'larger drops conc: " << (tmp_w * rhod)(x,y,z) * 1e-6 <<"'" << endl;
    gp << "plot"
       << "'-' with line title 'wet radius' lw 3 lc rgb 'blue'," << endl;
//     << "'-' with histeps title 'dry radius' lw 1 lc rgb 'red' " << endl;
    gp.send(focus_w);
//    gp.send(focus_d);

  }
}



/*

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/paper_GMD_2015/fig_a/",
    h5  = dir + "out_lgrngn",
    svg = dir + "out_lgrngn_spec.svg";

  Gnuplot gp;

  float ymin = .4 * .01, ymax = .9 * 10000;
  const int at = 9000;

  gp << "set term svg dynamic enhanced fsize 15 size 900, 1500 \n";
  gp << "set output '" << svg << "'\n";
  gp << "set logscale xy\n";
  gp << "set xrange [.002:100]\n";
  gp << "set yrange [" << ymin << ":" << ymax << "]\n";
  gp << "set ylabel '[mg^{-1} μm^{-1}]'\n"; // TODO: add textual description (PDF?)
  gp << "set grid\n";
  gp << "set nokey\n";

  // FSSP range
  gp << "set arrow from .5," << ymin << " to .5," << ymax << " nohead\n";
  gp << "set arrow from 25," << ymin << " to 25," << ymax << " nohead\n";

  gp << "set xlabel offset 0,1.5 'particle radius [μm]'\n";
  gp << "set key samplen 1.2\n";
  gp << "set xtics rotate by 65 right (.01, .1, 1, 10, 100) \n";

// TODO: use dashed lines to allow printing in black and white... same in image plots

  assert(focus.first.size() == focus.second.size());
  gp << "set multiplot layout " << focus.first.size() << ",2 columnsfirst upwards\n";

}
*/

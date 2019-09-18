#include "common.hpp"
#include "PlotterMicro.hpp"
#include "plots.hpp"
#include "focus.hpp"
#include <map>

namespace
{
  po::variables_map vm;
};

//plot spectra positions
template<class Plotter_t>
void plot_lgrngn_spec_positions(Plotter_t plotter)
{
  // read opts
  po::options_description opts("spectra plotting options");
  opts.add_options()
    ("spectra_step", po::value<int>()->required() , "time step number at which spectra are plotted")
  ;
  handle_opts(opts, vm);
  const int spectra_step = vm["spectra_step"].as<int>();

  auto& n = plotter.map;
  Gnuplot gp;

  try{
  // cloud water content
  auto tmp = plotter.h5load_ract_timestep(spectra_step * n["outfreq"]) * 1e3;

  std::string title = "cloud water mixing ratio [g/kg]";
  gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(spectra_step) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
  init(gp, plotter.file + "_spectra_positions_at" + zeropad(spectra_step) + ".svg", 1, 1, n, 2.); 

  {   
    for (auto &fcs : focus_3d)
    {   
      auto &x = fcs[0];
      auto &y = fcs[2];

      // black square
      gp << "set arrow from " << x-box_size << "," << y-box_size << " to " << x+(box_size+1) << "," << y-box_size << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
      gp << "set arrow from " << x-box_size << "," << y+(box_size+1) << " to " << x+(box_size+1) << "," << y+(box_size+1) << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
      gp << "set arrow from " << x-box_size << "," << y-box_size << " to " << x-box_size << "," << y+(box_size+1) << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
      gp << "set arrow from " << x+(box_size+1) << "," << y-box_size << " to " << x+(box_size+1) << "," << y+(box_size+1) << " nohead lw 4 lc rgbcolor '#ffffff' front\n";
      // white square
      gp << "set arrow from " << x-box_size << "," << y-box_size << " to " << x+(box_size+1) << "," << y-box_size << " nohead lw 2 front\n";
      gp << "set arrow from " << x-box_size << "," << y+(box_size+1) << " to " << x+(box_size+1) << "," << y+(box_size+1) << " nohead lw 2 front\n";
      gp << "set arrow from " << x-box_size << "," << y-box_size << " to " << x-box_size << "," << y+(box_size+1) << " nohead lw 2 front\n";
      gp << "set arrow from " << x+(box_size+1) << "," << y-box_size << " to " << x+(box_size+1) << "," << y+(box_size+1) << " nohead lw 2 front\n";
    }   
  }   

  // labels
  {   
    char lbl = 'a';
    for (auto &fcs : focus_3d)
    {   
      auto &x = fcs[0];
      auto &y = fcs[2];
      gp << "set label " << int(lbl) << " '" << lbl << "' at " << x-1 << "," << y-7 << " front font \",20\"\n";
      lbl += 1;
    }
  }   

  plotter.plot(gp, tmp, blitz::Range(yslice_idx, yslice_idx));
  }
  catch(...){}
}



//plot spectra
template<class Plotter_t>
void plot_lgrngn_spec(Plotter_t plotter)
{
  const int spectra_step = vm["spectra_step"].as<int>();

  const double r_c_adiab = 7.8; // adiabatic water content [g/kg] coming from Wojtek
  auto& n = plotter.map;
  for(auto elem : n)
  {
     std::cout << elem.first << " " << elem.second << std::endl;
  }
  Gnuplot gp;
//  int off = 2; // TODO!!!
  int off = 0; // TODO!!!
  string file = plotter.file + "_spectra_at" + zeropad(spectra_step) + ".svg";

  int hor = min<int>(focus_3d.size(), 5);
  int ver = double(focus_3d.size()) / 5. + 0.99999;

  init_prof(gp, file, ver, hor);


  gp << "set xrange [0:30]\n";
  gp << "set yrange [0:50]\n";
  gp << "set xlabel 'r [μm]'\n"; 

  // read in density
  auto tmp = plotter.h5load(plotter.file + "/const.h5", "G");
  typename Plotter_t::arr_t rhod(tmp);

  // focus to the gridbox from where the size distribution is plotted
  char lbl = 'a';
  bool setylabel = 1;
  for (auto &fcs : focus_3d)
  {
    if(setylabel)
    {
      gp << "set ylabel 'n [cm^{-3} μm^{-1}]'\n"; 
      setylabel = 0;
    }
    else
      gp << "unset ylabel\n"; 
    const int &x = fcs[0], &y = fcs[1], &z = fcs[2];
    const blitz::RectDomain<3> focusBox({x-box_size, y-box_size, z-box_size}, {x+box_size, y+box_size, z+box_size});

    gp << "set label 1 '(" << lbl << ")' at graph .1, .93 font \",15\"\n";
    lbl += 1;

    // calc ratio of water content to adiabatic water content
    {
      auto tmp = plotter.h5load_ract_timestep(spectra_step * n["outfreq"]) * 1e3;
      double ratio = mean(tmp(focusBox)) / r_c_adiab;
      gp << "set label 4 'AF = " << std::fixed << std::setprecision(2) << ratio << "' at graph .2, .63 font \",15\"\n";
    }

    // calc mean and std dev of radius of acivated droplets in the box
    double act_conc = 0.; // concentration of activated droplets [1/kg]
    double act_rw_mean = 0.; 
    double act_rw_std_dev = 0.; 
    try
    {
      auto tmp = plotter.h5load_timestep("actrw_rw_mom0", spectra_step * n["outfreq"]);
      typename Plotter_t::arr_t snap(tmp); 
      act_conc = mean(snap(focusBox));
    }
    catch(...) {;}
    // mean droplet radius in the box
    try
    {
      auto tmp = plotter.h5load_timestep("actrw_rw_mom1", spectra_step * n["outfreq"]);
      typename Plotter_t::arr_t snap(tmp); // 1st raw moment / mass [m / kg]
      if(act_conc > 0)
        act_rw_mean = mean(snap(focusBox)) / act_conc;
      else
        act_rw_mean = 0.;
    }
    catch(...) {;}
    // std deviation of distribution of radius in the box
    try
    {
      auto tmp = plotter.h5load_timestep("actrw_rw_mom0", spectra_step * n["outfreq"]);
      typename Plotter_t::arr_t zeroth_raw_mom(tmp); // 0th raw moment / mass [1 / kg]
      tmp = plotter.h5load_timestep("actrw_rw_mom1", spectra_step * n["outfreq"]);
      typename Plotter_t::arr_t first_raw_mom(tmp); // 1st raw moment / mass [m / kg]
      tmp = plotter.h5load_timestep("actrw_rw_mom2", spectra_step * n["outfreq"]);
      typename Plotter_t::arr_t second_raw_mom(tmp); // 2nd raw moment / mass [m^2 / kg]
      tmp = plotter.h5load_timestep("sd_conc", spectra_step * n["outfreq"]);
      typename Plotter_t::arr_t sd_conc(tmp); // number of SDs
      if(act_conc > 0)
      {
        double SD_no = sum(sd_conc(focusBox));
        if(SD_no > 1 && act_rw_mean > 0)
        {
          act_rw_std_dev = (
            SD_no / (SD_no - 1) /
            act_conc * (
              mean(second_raw_mom(focusBox)) -
              2. * act_rw_mean * mean(first_raw_mom(focusBox)) +
              act_rw_mean * act_rw_mean * mean(zeroth_raw_mom(focusBox))
            )
          );

          // could not be true due to numerics?
          if(act_rw_std_dev > 0.)
            act_rw_std_dev = sqrt(act_rw_std_dev);
          else
            act_rw_std_dev = 0.;
        }
      }
      else
        act_rw_std_dev = 0.;
    }
    catch(...) {;}


    gp << "set label 2 '<r> = " << std::fixed << std::setprecision(2) << act_rw_mean * 1e6 <<" µm' at graph .2, .83 font \",15\"\n";
    gp << "set label 3 'σ = " << std::fixed << std::setprecision(2) << act_rw_std_dev * 1e6 <<" µm' at graph .2, .73 font \",15\"\n";


    std::map<float, float> focus_w;

    //info on the number and location of histogram edges
    auto left_edges_rw = bins_wet();
    int nsw = left_edges_rw.size() - 1;

    for (int i = 0; i < nsw; ++i)
    {
      const string name = "rw_rng" + zeropad(i + off) + "_mom0";
      auto tmp_w = plotter.h5load_timestep(name, spectra_step * n["outfreq"]);

      focus_w[(left_edges_rw[i] + left_edges_rw[i+1]) / (setup::real_t(2. * 1e-6) * si::metres)] =
        sum((tmp_w*rhod)(focusBox))
        * 1e-6 // per cm^{-3}
        / pow(2*box_size+1,3) 
        / ((left_edges_rw[i+1] - left_edges_rw[i]) / 1e-6 / si::metres); // per micrometre
    }
    const string name = "rw_rng" + zeropad(nsw + off) + "_mom0";
    auto tmp_w = plotter.h5load_timestep(name, spectra_step * n["outfreq"]);

    notice_macro("setting-up plot parameters");
    std::cout << "larger drops conc: " << (tmp_w * rhod)(x,y,z) * 1e-6 << endl;
//    gp << "set title 'larger drops conc: " << (tmp_w * rhod)(x,y,z) * 1e-6 <<"'" << endl;
    gp << "plot"
       << "'-' with histeps title 'wet radius' lw 3 lc rgb 'blue'," << endl;
//     << "'-' with histeps title 'dry radius' lw 1 lc rgb 'red' " << endl;
    gp.send(focus_w);
//    gp.send(focus_d);

  }
}

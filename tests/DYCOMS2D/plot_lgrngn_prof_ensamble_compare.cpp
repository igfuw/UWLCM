#include <blitz/array.h>
#include <fstream>
#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <boost/tuple/tuple.hpp>


BZ_USING_NAMESPACE(blitz)

int main(int argc, char* argv[]) 
{
  Array<double, 1> snap;
  Array<double, 1> res[2];
  Array<double, 1> res_up;
  Array<double, 1> res_down;
  blitz::firstIndex fi;
  double xmax;

  std::string prof_file_name="/out_lgrngn_mean_profiles.dat";
  std::set<std::string> profs({"00rtot", "rliq", "thl", "wvar", "w3rd", "prflux", "act_conc", "clfrac", "N_c", "non_gccn_rw_up", "gccn_rw_up", "gccn_rw_down", "non_gccn_rw_down", "sat_RH"});

  blitz::Array<float, 1> res_pos;

  int prof_ctr = 0;
  for (auto &plt : profs) 
  {
    Gnuplot gp;
    std::string file = argv[1] + std::string("/profiles_compare") + plt + std::string(".svg");
    std::cout << "out file: " << file << std::endl;
    init_prof(gp, file, 1, 1);
    gp << "set ylabel 'height / inversion height'\n";
    gp << "set key right bottom\n";

    if (plt == "rliq")
      gp << "set xlabel 'liquid water mixing ratio [g/kg]'\n";
    else if (plt == "00rtot")
      gp << "set xlabel 'total water mixing ratio [g/kg]'\n";
    else if (plt == "N_c")
      gp << "set xlabel 'cloud droplets concentration [1/cm^3]'\n";
    else if (plt == "thl")
      gp << "set xlabel 'liquid potential temp [K]'\n";
    else if (plt == "clfrac")
      gp << "set xlabel 'cloud fraction'\n";
    else if (plt == "prflux")
    {
      gp << "set xlabel 'precipitation flux [W/m^2]'\n";
      xmax=40;
    }
    else if (plt == "wvar")
      gp << "set xlabel 'variance of w [m^2 / s^2]'\n";
    else if (plt == "w3rd")
      gp << "set xlabel '3rd mom of w [m^3 / s^3]'\n";
    else if (plt == "act_rd")
    {
      gp << "set xlabel 'activated droplets mean dry radius [um]'\n";
      gp << "set xrange [0:0.15]\n";
    }
    else if (plt == "gccn_rw")
    {
      gp << "set xlabel 'GCCN-based droplets mean wet radius [um]'\n";
    }
    else if (plt == "gccn_rw_up")
    {
      gp << "set xlabel 'GCCN mean wet radius [μm] in updrafts'\n";
      gp << "set nokey\n";
      xmax = 35;
    }
    else if (plt == "non_gccn_rw_up")
    {
      gp << "set xlabel 'non-GCCN mean wet radius [μm] in updrafts'\n";
      gp << "set nokey\n";
      xmax = 7;
    }
    else if (plt == "gccn_rw_down")
    {
      gp << "set xlabel 'GCCN mean wet radius [μm] in downdrafts'\n";
      gp << "set nokey\n";
      xmax = 30;
    }
    else if (plt == "non_gccn_rw_down")
    {
      gp << "set xlabel 'non-GCCN mean wet radius [μm] in downdrafts'\n";
      gp << "set nokey\n";
      xmax = 6;
    }
    else if (plt == "sat_RH")
    {
      gp << "set title 'supersaturation in updrafts'\n";
      gp << "set yrange [0.45:1.]\n";
      gp << "set xrange [0.000:0.005]\n";
    }
 
    gp << "set arrow from 0, 0.48 to graph 1, first 0.48 nohead\n";


    const int lines_per_plot = plt == "act_conc" ? 2 : 1;

    for(int j=2; j<argc; j+=2)
      for(int lpp=0; lpp<lines_per_plot;++lpp)
      {
        if(j==2 && lpp==0)
          gp << "plot '-' with line t '" << argv[j] << "'";
        else
          gp << ", '-' with line t '" << argv[j] << "'";
      }
    gp << "\n";


    for(int i=1; i<argc; i+=2)
    {
      std::string file = argv[i] + prof_file_name;
      std::cout << "reading in " << i << ": " << file << std::endl;
      ifstream iprof_file(file);
      if (iprof_file.bad())
      {
        cerr << "Unable to open file: " << file << endl;
        continue;
      }
      // pass preceeding profile
      for(int k=0; k<prof_ctr; ++k)
        iprof_file >> snap;

      // read in interesting profiles
      for(int lpp=0; lpp<lines_per_plot;++lpp)
      {
        iprof_file >> snap;
        res[lpp].resize(snap.shape());
        res_pos.resize(snap.shape());
        res[lpp] = snap;
      }

      // pass rest of profiles
      for(int k=prof_ctr+lines_per_plot; k<profs.size()+1; ++k)
        iprof_file >> snap;
      double z_i;
      // read in inversion height
      iprof_file >> z_i;
      std::cout << z_i << std::endl;
      res_pos = fi * 5. / z_i; // hardcoded dz = 5m !!
      std::cout << i << " " << plt << " res_pos: " << res_pos;
      for(int lpp=0; lpp<lines_per_plot;++lpp)
      {
        std::cout << i << " " << plt << " res: " << res[lpp];
        gp.send1d(boost::make_tuple(res[lpp], res_pos));
      }
    }
    if(plt == "rv" || plt == "sat" || plt == "sat_RH" || plt == "gccn_rw_down"|| plt == "gccn_rw_up" || plt == "non_gccn_rw_down" || plt == "non_gccn_rw_up")
    {
      gp << "set yrange [0.:1.2]\n";
      gp << "set xrange [*:*]\n";
    }
    prof_ctr+=lines_per_plot;


    if(argc == 3) // if its only a GCCN plot, add the plot of up_down rw
    {
      res_up.resize(res[0].shape());
      res_down.resize(res[0].shape());
      if(plt == "gccn_rw_down" || plt == "non_gccn_rw_down")
        res_down = res[0];
      if(plt == "gccn_rw_up" || plt == "non_gccn_rw_up")
      {
        Gnuplot gp;
        std::string title;
        title = plt == "gccn_rw_up" ? "gccn_rw_up_down" : "non_gccn_rw_up_down";
        std::string file = argv[1] + std::string("/profiles_compare") + title + std::string(".svg");
        std::cout << "out file: " << file << std::endl;
        init_prof(gp, file, 1, 1);
        res_up = res[0];
        std::cout << res_pos << res_up << res_down;
        title = plt == "gccn_rw_up" ? "" : "non-";
        gp << "set xlabel '" << title << "GCCN mean wet radius [μm]'\n";
        gp << "set ylabel 'height / inversion height'\n";
        gp << "set arrow from 0, 0.48 to " << xmax << ", 0.48 nohead\n";
        gp << "plot '-' w l t 'updraft'";
        gp << ", '-' w l t 'downdraft'\n";
        gp.send1d(boost::make_tuple(res_up, res_pos));
        gp.send1d(boost::make_tuple(res_down, res_pos));
        gp << "set xrange [*:*]\n";
      }    
    }      
  }        
}          
           
           
           
           
           
           
           
           
           
           
           
           
           

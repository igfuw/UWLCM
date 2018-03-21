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
  Array<double, 1> res;
  blitz::firstIndex fi;

  std::string prof_file_name="/out_lgrngn_mean_series.dat";
  std::set<std::string> profs({"wvarmax", "clfrac", "lwp", "er", "surf_precip", "mass_dry", "acc_precip", "cl_nc"});

  blitz::Array<float, 1> res_pos;

  int prof_ctr = 0;
  for (auto &plt : profs) 
  {
    Gnuplot gp;
    std::string file = argv[1] + std::string("/series_compare_") + plt + std::string(".svg");
    std::cout << "out file: " << file << std::endl;
    init_prof(gp, file, 1, 1);

    gp << "set yrange[*:*]\n";
    gp << "set xrange[*:*]\n";
    gp << "set xlabel 'time [h]'\n";

    if (plt == "clfrac")
      gp << "set ylabel 'average cloud fraction'\n";
    else if (plt == "nc")
      gp << "set ylabel 'average cloud drop concentration [1/cm^3]'\n";
    else if (plt == "cl_nc")
      gp << "set ylabel 'average cloud drop concentration [1/cm^3] in cloudy cells'\n";
    else if (plt == "wvarmax")
      gp << "set ylabel 'max variance of w [m^2 / s^2]'\n";
    else if (plt == "surf_precip")
      gp << "set ylabel 'surface precipitation [mm/d]'\n";
    else if (plt == "acc_precip")
    {
      gp << "set key left top\n";
      gp << "set ylabel 'accumulated surface precipitation [mm]'\n";
    }
    else if (plt == "mass_dry")
      gp << "set ylabel 'total dry mass [g]'\n";
    else if (plt == "lwp")
      gp << "set ylabel 'liquid water path [g / m^2]'\n";
    else if (plt == "er")
      gp << "set ylabel 'entrainment rate [cm / s]'\n";

    gp << "plot '-' with line t '" << argv[2] << "'";
    for(int j=3; j<argc; j+=2)
      gp << ", '-' with line t '" << argv[j+1] << "'";
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
      for(int k=0; k<=prof_ctr; ++k)
        iprof_file >> snap;

      res.resize(snap.shape());
      res_pos.resize(snap.shape());
      res = snap;

      for(int k=prof_ctr+1; k<profs.size(); ++k)
        iprof_file >> snap;

      auto n = h5n(argv[i]);
      res_pos = fi * n["outfreq"] * n["dt"] / 3600.;
      std::cout << i << " " << plt << " res: " << res;
      std::cout << i << " " << plt << " res_pos: " << res_pos;
      gp.send1d(boost::make_tuple(res_pos, res));
    }
    if (plt == "acc_precip")
      gp << "set key right top'\n";
    prof_ctr++;
  }
}

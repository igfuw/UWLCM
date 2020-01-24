//#include <blitz/array.h>
//#include <fstream>
#include "common.hpp"
//#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include "plots.hpp"
#include "gnuplot_series_set_labels.hpp"
#include "gnuplot_profs_set_labels.hpp"

using namespace blitz;

enum whattoplot {series, profs};

// taken from https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

bool open_file(string filename, ifstream &stream)
{
  std::cout << "reading in " << filename << std::endl;
  stream.open(filename);
  if (stream.bad())
  {
    cerr << "Unable to open file: " << filename << endl;
    throw std::runtime_error("error opening file");
  }
}

void average(int argc, char* argv[], int wtp, std::vector<std::string> types, std::string suffix)
{
  Array<double, 1> snap, snap2;
  Array<double, 1> avg, pos, weight, std_dev;
  blitz::firstIndex fi;
  ifstream iprof_file;
  std::string line;

  std::string ofile(argv[1]);
  ofile.append(suffix);
  ofstream oprof_file(ofile + string(".dat"));

  // init the plotter
  Gnuplot gp;
  int hor = min<int>(types.size(), 4);
  int ver = double(types.size()) / 4. + 0.99999;
  init_prof(gp, ofile + string(".svg"), ver, hor); 

  types.insert(types.begin(), "position");

  int prof_ctr = 0;

  // init sizes of profiles/series assuming that they have same size in all files
  open_file(string(argv[4]).append(suffix).append(".dat"), iprof_file);
  std::getline(iprof_file, line); // discard array description
  iprof_file >> snap;
        cerr << "size init: " << line << endl;
        cerr << "size init: " << snap;
  iprof_file.close();
  avg.resize(snap.shape());
  std_dev.resize(snap.shape());
  pos.resize(snap.shape());
  weight.resize(snap.shape());
  snap2.resize(snap.shape());

  // actual averaging
  for (auto &plt : types)
  {
    // set gnuplot plot labels
    if(wtp == series)
      gnuplot_series_set_labels(gp, plt);
    else if (wtp == profs)
      gnuplot_profs_set_labels(gp, plt, false); // FIXME: we assume here that normalization is not done

    avg = 0;
    std_dev = 0;
    weight = 0;

    for(int i=4; i<argc; i+=1) // add value of the array from each file
    {
   //   try{
        open_file(string(argv[i]).append(suffix).append(".dat"), iprof_file);
     // } catch (...) {continue;}

      // discard arrays tored before the desired one
      for(int k=0; k<=prof_ctr; ++k)
      {
        std::getline(iprof_file, line); // discard array description
        cerr << "plt: " << plt << " " << line << endl;
        iprof_file >> snap; // read in the array
        cerr << "plt: " << plt  << " " << snap;
        std::getline(iprof_file, line); // blitz reading arrays doesnt move to the next line, need to do it manually here
      }

      if (plt == "base_prflux_vs_clhght") // we need weighted average for this plot - array of weights before array of values
      {
        weight += snap; // add to sum of weights
        snap2 = snap; // store weight of this one
        std::getline(iprof_file, line); // discard array description
        iprof_file >> snap; // read actual values
        avg += snap * snap2; // add weight*value
        std_dev += snap * snap * snap2; // weight * value^2
      }
      else // all values have same weight
      {
        avg += snap;
        std_dev += snap * snap;
        weight += 1;
      }
      iprof_file.close();
    }
    avg = where(weight > 0, avg / weight, 0);
    //std_dev = where(weight > 0, sqrt(std_dev / weight - avg * avg), 0); // without the Bessel correction
    std_dev = where(weight > 0, sqrt(weight / (weight - 1) * (std_dev / weight - avg * avg)), 0); // with the Bessel correction
    if(plt == "position")
      pos = avg;
    else
    {
      gp << "plot '-' with l \n";
      if(wtp == series)
        gp.send1d(boost::make_tuple(pos, avg));
      else if (wtp == profs)
        gp.send1d(boost::make_tuple(avg, pos));
    }

    std::cout << plt << " avg: " << avg;
    std::cout << plt << " std_dev: " << std_dev;
    if (plt == "base_prflux_vs_clhght") 
    {
      oprof_file << plt << " number of occurances" << endl;
      oprof_file << weight;
    }
    oprof_file << plt << endl;
    oprof_file << avg;
    oprof_file << plt << "_std_dev" << endl;
    oprof_file << std_dev;

    prof_ctr += 
      plt == "base_prflux_vs_clhght" ||  // this plot outputs two arrays: weights and averages
      plt == "ract_avg" ||               // following plots also output two arrays: averages and standard deviations
      plt == "cloud_avg_act_conc" ||
      plt == "cloud_avg_supersat" ||
      plt == "sd_conc_avg" ||
      plt == "sd_conc_act_avg"
      ? 2 : 1;
  }
}

int main(int argc, char* argv[])
{
 const string profiles_suffix = string("_profiles_")+argv[2]+string("_")+argv[3];
 string plot_type;
 // determine type of plots based on the name of the first file
 const string types[] = {"rico", "dycoms", "moist_thermal"};
 for(auto type : types)
 {
   if(hasEnding(string(argv[4]), type))
   {
     plot_type = type;
     break;
   }
 }

 cout << "Detected plot type: " << plot_type << endl;

 Plots plots(plot_type, true); // assume sgs is on

 try{ average(argc, argv, series, plots.series, "_series"); } catch(...){;}
 try{ average(argc, argv, profs, plots.profs, profiles_suffix);} catch(...){;}
// try{ average(argc, argv, {"base_prflux_vs_clhght"}, profiles_suffix); } catch(...){;} // TODO: add this profile to rico profiles, or give some other suffix here to distinguish outputted average from regular profiles
}

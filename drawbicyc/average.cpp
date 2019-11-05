//#include <blitz/array.h>
//#include <fstream>
#include "common.hpp"
//#include "PlotterMicro.hpp"
#include <boost/tuple/tuple.hpp>
#include "plots.hpp"

using namespace blitz;

void average(int argc, char* argv[], std::vector<std::string> types, std::string suffix)
{
  Array<double, 1> snap;
  Array<double, 1> avg;
  blitz::firstIndex fi;

  std::string ofile(argv[1]);
  ofile.append(suffix);//"_series.dat");
  ofstream oprof_file(ofile);

  types.insert(types.begin(), "position");

  int prof_ctr = 0;
  for (auto &plt : types)
  {
    avg = 0;
    int opened_files = 0;
    for(int i=2; i<argc; i+=1)
    {
      std::string file(argv[i]);
      file.append(suffix);//"_series.dat");
      std::cout << "reading in " << i << ": " << file << std::endl;
      ifstream iprof_file(file);
      if (iprof_file.bad())
      {
        cerr << "Unable to open file: " << file << endl;
        continue;
      }
      for(int k=0; k<=prof_ctr; ++k)
        iprof_file >> snap;

      if(i==2)
      {
        avg.resize(snap.shape());
        avg = 0.;
      }
      avg += snap;
      opened_files++;
    }
    avg /= opened_files;
    std::cout << plt << " avg: " << avg;
    oprof_file << avg;
    prof_ctr++;
  }
}

int main(int argc, char* argv[])
{
  average(argc, argv, series_dycoms, "_series.dat");
  average(argc, argv, profs_dycoms, "_profiles_7200_21600.dat");
}

#include"plot_series.hpp"
#include"plot_prof.hpp"
#include"plot_fields.hpp"

int main(int argc, char** argv)
{
  // make args global
  ac=argc;
  av=argv;

  // general opts
  opts_main.add_options()
    ("profs", po::value<bool>()->default_value(true), "plot profiles?")
    ("series", po::value<bool>()->default_value(true) , "plot series?")
    ("fields", po::value<bool>()->default_value(false) , "plot fields?")
    ("dir", po::value<std::string>()->required() , "directory containing out_lgrngn")
    ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
    ("type", po::value<std::string>()->required(), "one of: dycoms, moist_thermal")
  ;

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); //     ignores unknown

  // checking if all required options present
  po::notify(vm);

  // handling the "micro" option
  std::string micro = vm["micro"].as<std::string>();

  // handling the "type" option
  std::string type = vm["type"].as<std::string>();

  // parse dir name
  std::string
    dir = vm["dir"].as<std::string>(),
    h5  = dir + "out_" + micro;

  // reading required plot types
  bool flag_series = vm["series"].as<bool>(),
       flag_profiles = vm["profs"].as<bool>(),
       flag_fields = vm["fields"].as<bool>();

  // detecting input data dimensionality
  H5::H5File h5f(h5 + "/const.h5", H5F_ACC_RDONLY);
  H5::DataSet h5d = h5f.openDataSet("G");
  H5::DataSpace h5s = h5d.getSpace();
  int NDims = h5s.getSimpleExtentNdims();

  if(NDims == 2)
  {
    if(flag_series)   plot_series(PlotterMicro_t<2>(h5, micro), Plots(type));
    if(flag_profiles) plot_profiles(PlotterMicro_t<2>(h5, micro), Plots(type));
    if(flag_fields)   plot_fields(PlotterMicro_t<2>(h5, micro), Plots(type));
  }
  else if(NDims == 3)
  {
    if(flag_series)   plot_series(PlotterMicro_t<3>(h5, micro), Plots(type));
    if(flag_profiles) plot_profiles(PlotterMicro_t<3>(h5, micro), Plots(type));
    if(flag_fields)   plot_fields(PlotterMicro_t<3>(h5, micro), Plots(type));
  }
  else
    assert(false && "need 2d or 3d input data");

return 0;
} // main

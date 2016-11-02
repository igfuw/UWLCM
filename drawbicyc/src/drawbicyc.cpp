#include"plot_series.hpp"
#include"plot_prof.hpp"

int main(int argc, char** argv)
{
  // make args global
  ac=argc;
  av=argv;

  // general opts
  opts_main.add_options()
    ("prof", po::value<bool>()->default_value(true), "plot profiles?")
    ("series", po::value<bool>()->default_value(true) , "plot series?")
    ("dir", po::value<std::string>()->required() , "directory containing out_lgrngn")
    ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
  ;

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); //     ignores unknown

  // checking if all required options present
  po::notify(vm);

  // handling the "micro" option
  std::string micro = vm["micro"].as<std::string>();

  // parse dir name
  std::string
    dir = vm["dir"].as<std::string>(),
    h5  = dir + "out_" + micro;

  // reading required plot types
  bool flag_series = vm["series"].as<bool>(),
       flag_profiles = vm["prof"].as<bool>();

  // detecting input data dimensionality
  H5::H5File h5f(h5 + "/const.h5", H5F_ACC_RDONLY);
  H5::DataSet h5d = h5f.openDataSet("G");
  H5::DataSpace h5s = h5d.getSpace();
  int NDims = h5s.getSimpleExtentNdims();

  if(NDims == 2)
  {
    if(flag_series)   plot_series(PlotterMicro_t<2>(h5, micro));
    if(flag_profiles) plot_profiles(PlotterMicro_t<2>(h5, micro));
  }
  else if(NDims == 3)
  {
    if(flag_series)   plot_series(PlotterMicro_t<3>(h5, micro));
    if(flag_profiles) plot_profiles(PlotterMicro_t<3>(h5, micro));
  }
  else
    assert(false && "need 2d or 3d input data");

return 0;
} // main

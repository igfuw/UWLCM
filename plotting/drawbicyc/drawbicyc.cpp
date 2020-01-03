#include"plot_series.hpp"
#include"plot_prof.hpp"
#include"plot_fields.hpp"
#include"plot_qv_qc_2_6_10_min.hpp"
#include"plot_lgrngn_spec.hpp"

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
    ("spectra", po::value<bool>()->default_value(false) , "plot spectra?")
    ("qv_qc_2_6_10_min", po::value<bool>()->default_value(false) , "plot comparison of qv and qc fields at 2, 6 and 10 min?")
    ("dir", po::value<std::string>()->required() , "directory containing out_lgrngn")
    ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
    ("type", po::value<std::string>()->required(), "one of: dycoms, moist_thermal, rico")//, base_prflux_vs_clhght")
  ;

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); //     ignores unknown

  // checking if all required options present
  po::notify(vm);

  // handling the "micro" option
  std::string micro = vm["micro"].as<std::string>();

  // handling the "type" option
  std::string type = vm["type"].as<std::string>();
  if(type != "dycoms" && type != "moist_thermal" && type != "rico")// && type != "base_prflux_vs_clhght")
    throw std::runtime_error("Unrecognized 'type' option, only dycoms, rico, moist_thermal available now");//, base_prflux_vs_clhght available now");

  // should profiles be normalized by inversion height
  const bool normalize_prof = type == "dycoms";

  // parse dir name
  std::string
    dir = vm["dir"].as<std::string>(),
    h5  = dir + "out_" + micro;

  // reading required plot types
  bool flag_series = vm["series"].as<bool>(),
       flag_profiles = vm["profs"].as<bool>(),
       flag_fields = vm["fields"].as<bool>(),
       flag_lgrngn_spec = vm["spectra"].as<bool>(),
       flag_qv_qc_2_6_10_min = vm["qv_qc_2_6_10_min"].as<bool>();

  // detecting input data dimensionality
  H5::H5File h5f(h5 + "/const.h5", H5F_ACC_RDONLY);
  H5::DataSet h5d = h5f.openDataSet("G");
  H5::DataSpace h5s = h5d.getSpace();
  int NDims = h5s.getSimpleExtentNdims();
  
  // detecting if subgrid model was on
  bool sgs = true;
  try 
  {
    auto h5g = h5f.openGroup("sgs");
  }
  catch (...)
  {
    sgs = false;
  }

  Plots plots(type, sgs);

  if(NDims == 2)
  {

    if(flag_series)   plot_series(PlotterMicro_t<2>(h5, micro), plots, type);
    if(flag_profiles) plot_profiles(PlotterMicro_t<2>(h5, micro), plots, type, normalize_prof);
    if(flag_fields)   plot_fields(PlotterMicro_t<2>(h5, micro), plots, type);
    if(flag_qv_qc_2_6_10_min)   plot_qv_qc_2_6_10_min(PlotterMicro_t<2>(h5, micro));
  }
  else if(NDims == 3)
  {
    if(flag_series)   plot_series(PlotterMicro_t<3>(h5, micro), plots, type);
    if(flag_profiles) plot_profiles(PlotterMicro_t<3>(h5, micro), plots, type, normalize_prof);
    if(flag_fields)   plot_fields(PlotterMicro_t<3>(h5, micro), plots, type);
    if(flag_qv_qc_2_6_10_min)   plot_qv_qc_2_6_10_min(PlotterMicro_t<2>(h5, micro));
    if(flag_lgrngn_spec) {
      plot_lgrngn_spec_positions(PlotterMicro_t<3>(h5, "lgrngn"));
      plot_lgrngn_spec(PlotterMicro_t<3>(h5, "lgrngn"));
    }
  }
  else
    assert(false && "need 2d or 3d input data");

return 0;
} // main

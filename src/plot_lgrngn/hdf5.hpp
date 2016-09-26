#pragma once

#include <blitz/array.h>
#include <H5Cpp.h>
#include <map>


class PlotterCommon
{
  public:
  std::map<std::string, double> map;

  protected:
  H5::H5File h5f;
  H5::DataSet h5d;
  H5::DataSpace h5s;

  void h5load(
    const string &file, 
    const string &dataset
  )
  {
    notice_macro("about to open file: " << file)
    h5f.openFile(file, H5F_ACC_RDONLY);

    notice_macro("about to read dataset: " << dataset)
    h5d = h5f.openDataSet(dataset);
    h5s = h5d.getSpace();
  }

  //ctor
  PlotterCommon(const string &file)
  {
    // init dt and outfreq
    {
      h5load(file + "/const.h5", "T");

      float dt;
      {
        auto attr = h5d.openAttribute("dt");
        attr.read(attr.getDataType(), &dt);
      }
      map["dt"] = dt;

      const hsize_t two = 2, zero = 0;
      float tmp[2];
      h5s.selectHyperslab( H5S_SELECT_SET, &two, &zero);
      h5d.read(tmp, H5::PredType::NATIVE_FLOAT, H5::DataSpace(1, &two), h5s);
      map["outfreq"] = (tmp[1] - tmp[0]) / dt;
    }
  }
};

template<int NDims>
class Plotter_t : public PlotterCommon {};

template<>
class Plotter_t<2> : public PlotterCommon 
{
  protected:
  using parent_t = PlotterCommon;
  hsize_t n[2];
  enum {x, z};
  blitz::Array<float, 2> tmp;

  public:

  auto h5load(
    const string &file, 
    const string &dataset
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 2>() + 0))
  {
    parent_t::h5load(file, dataset);
    this->h5s.getSimpleExtentDims(n, NULL);
  
    hsize_t 
      cnt[2] = { n[x],  n[z] }, 
      off[2] = { 0,     0    };
    this->h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[2] = {
      hsize_t(tmp.extent(0)), 
      hsize_t(tmp.extent(1)) 
    };
    this->h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(2, ext), h5s);

    return blitz::safeToReturn(tmp + 0);
  }

  auto h5load_timestep(
    const string &file, 
    const string &dataset,
    int at
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 2>() + 0))
  {
    string timestep_file = file + "/timestep" + zeropad(at, 10) + ".h5";
    return h5load(timestep_file, dataset);
  }

  //ctor
  Plotter_t(const string &file):
    parent_t(file)
  {
    // read number of timesteps
    this->h5f.openDataSet("T").getSpace().getSimpleExtentDims(n, NULL);
    this->map["t"] = n[0];

    // read number of cells
    this->h5f.openDataSet("X").getSpace().getSimpleExtentDims(n, NULL); // X gives cell-border coordinates (+1)
    this->map["x"] = n[0]-1;
    this->map["z"] = n[1]-1;
    tmp.resize(n[0], n[1]);

    // read dx,dy,dz
    h5load(file + "/const.h5", "X");
    this->map["dx"] = tmp(1,0) - tmp(0,0);
    h5load(file + "/const.h5", "Y");
    this->map["dz"] = tmp(0,1) - tmp(0,0);
  }
};


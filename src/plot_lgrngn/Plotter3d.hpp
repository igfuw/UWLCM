#pragma once
#include "Plotter2d.hpp"

// 3d version
template<>
class Plotter_t<3> : public PlotterCommon 
{
  protected:
  using parent_t = PlotterCommon;
  hsize_t n[3];
  enum {x, y, z};
  blitz::Array<float, 3> tmp;

  public:

  auto h5load(
    const string &file, 
    const string &dataset
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 3>() + 0))
  {
    parent_t::h5load(file, dataset);
    this->h5s.getSimpleExtentDims(n, NULL);
  
    hsize_t 
      cnt[3] = { n[x],  n[y],  n[z] }, 
      off[3] = { 0,     0,     0    };
    this->h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[3] = {
      hsize_t(tmp.extent(0)), 
      hsize_t(tmp.extent(1)), 
      hsize_t(tmp.extent(2)) 
    };
    this->h5d.read(tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, ext), h5s);

    return blitz::safeToReturn(tmp + 0);
  }

  auto h5load_timestep(
    const string &file, 
    const string &dataset,
    int at
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 3>() + 0))
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
    this->map["y"] = n[1]-1;
    this->map["z"] = n[2]-1;
    tmp.resize(n[0], n[1], n[2]);

    // read dx,dy,dz
    h5load(file + "/const.h5", "X");
    this->map["dx"] = tmp(1,0,0) - tmp(0,0,0);
    h5load(file + "/const.h5", "Y");
    this->map["dy"] = tmp(0,1,0) - tmp(0,0,0);
    h5load(file + "/const.h5", "Z");
    this->map["dz"] = tmp(0,0,1) - tmp(0,0,0);
  }
};


#pragma once
#include "PlotterCommon.hpp"


template<int NDims>
class Plotter_t : public PlotterCommon {};

// 2d version
template<>
class Plotter_t<2> : public PlotterCommon 
{
  public:
  using arr_t = blitz::Array<float,2>;
  blitz::Array<int, 1> k_i;
  blitz::secondIndex LastIndex;

  protected:
  using parent_t = PlotterCommon;
  hsize_t n[2];
  enum {x, z};
  arr_t tmp;

  public:

  auto h5load(
    const string &file, 
    const string &dataset
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
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
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    string timestep_file = file + "/timestep" + zeropad(at, 10) + ".h5";
    return h5load(timestep_file, dataset);
  }

  auto horizontal_mean(
    const arr_t &data
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 1>() + 0))
  {
    auto tmp = blitz::mean(data(blitz::tensor::j, blitz::tensor::i), blitz::tensor::j);
    blitz::Array<float, 1> mean(tmp);
    return blitz::safeToReturn(mean + 0);
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
    k_i.resize(n[0]-1);
 
    // read dx,dy,dz
    h5load(file + "/const.h5", "X");
    this->map["dx"] = tmp(1,0) - tmp(0,0);
    h5load(file + "/const.h5", "Y");
    this->map["dz"] = tmp(0,1) - tmp(0,0);
    this->CellVol = this->map["dx"] * this->map["dz"];
    this->DomainSurf = this->map["dx"] * this->map["x"];

    // other dataset are of the size x*z, resize tmp
    tmp.resize(n[0]-1, n[1]-1);
  }
};


#pragma once
#include "common.hpp"
#include "PlotterCommon.hpp"

template<int NDims>
class Plotter_t : public PlotterCommon {};

// 2d version
template<>
class Plotter_t<2> : public PlotterCommon 
{
  public:
  static const int n_dims = 2;
  using arr_t = blitz::Array<float,2>;
  blitz::Array<int, 1> k_i;
  blitz::secondIndex LastIndex;

  protected:
  using parent_t = PlotterCommon;
  hsize_t n[2];
  enum {x, z};
  arr_t tmp, tmp_srfc;

  public:

  auto h5load(
    const string &file, 
    const string &dataset,
    bool srfc = false
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    parent_t::h5load(file, dataset);
    this->h5s.getSimpleExtentDims(n, NULL);
  
    hsize_t 
      cnt[2] = { n[x],  srfc ? 1 : n[z] }, 
      off[2] = { 0,     0    };
    this->h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[2] = {
      hsize_t(tmp.extent(0)), 
      hsize_t(srfc ? tmp_srfc.extent(1) : tmp.extent(1)) 
    };
    this->h5d.read(srfc ? tmp_srfc.data() : tmp.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(2, ext), h5s);

    return blitz::safeToReturn((srfc ? tmp_srfc : tmp) + 0);
  }

  auto h5load_timestep(
    const string &dataset,
    int at,
    bool srfc = false
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    string timestep_file = this->file + "/timestep" + zeropad(at, 10) + ".h5";
    return h5load(timestep_file, dataset, srfc);
  }

  auto horizontal_mean(
    const arr_t &data
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 1>() + 0))
  {
    auto tmp = blitz::mean(data(blitz::tensor::j, blitz::tensor::i), blitz::tensor::j);
    blitz::Array<float, 1> mean(tmp);
    return blitz::safeToReturn(mean + 0);
  }
  
  void subtract_horizontal_mean(
    arr_t &data
  )
  {
    blitz::Array<float, 1> mean(horizontal_mean(data));
    data = data(blitz::tensor::i, blitz::tensor::j) - mean(blitz::tensor::j);
  }

  auto horizontal_sum(
    const arr_t &data
  ) -> decltype(blitz::safeToReturn(blitz::Array<float, 1>() + 0))
  {
    auto tmp = blitz::sum(data(blitz::tensor::j, blitz::tensor::i), blitz::tensor::j);
    blitz::Array<float, 1> mean(tmp);
    return blitz::safeToReturn(mean + 0);
  }

  blitz::RectDomain<2> hrzntl_slice(const int &z)
  {
    return blitz::RectDomain<2>( blitz::TinyVector<blitz::Range, 2>(blitz::Range(0, this->map["x"]-1), blitz::Range(z,z)));
  }

  template <class gp_t, class data_t>
  void plot(gp_t &gp, const data_t &data)
  {
    parent_t::plot(gp);
    blitz::Array<float, 2> tmp(data);
  
//    gp << "set size 1.6\n";
    gp << "set xrange [0:" << tmp.extent(0)-1 << "]\n";
    gp << "set yrange [0:" << tmp.extent(1)-1 << "]\n";
//    gp << "set xrange [" << (tmp.extent(0)-1) * 1./3. <<":" << (tmp.extent(0)-1) * 2. / 3. << "]\n";
  //  gp << "set yrange [" << (tmp.extent(1)-1) * 1./5. <<":" << (tmp.extent(1)-1) * 3.5 / 5. << "]\n";

//  gp << "set xtics out scale .5 rotate by 60 ('1.2' 1200/dx, '1.5' 1500/dx, '1.8' 1800/dx, '2.1' 2100/dx, '2.4' 2400/dx)\n";//, '4.0' 4000/dx, '4.8' 4800/dx, '5.6' 5600/dx, '6.4' 6400/dx)\n"; 

  //gp << "set ytics out scale .5 rotate by 60 ('0.0' 0, '0.3' 300/dy, '0.6' 600/dy, '0.9' 900/dy, '1.2' 1200/dy, '1.5' 1500/dy)\n"; 

    gp << "splot '-' binary" << gp.binfmt(tmp.transpose(blitz::secondDim, blitz::firstDim)) << " scan=yx origin=(0,0,0) with image failsafe notitle\n";
    gp.sendBinary(tmp);
  }


  //ctor
  Plotter_t(const string &file):
    parent_t(file)
  {
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
    tmp_srfc.resize(n[0]-1, 1);
  }
};


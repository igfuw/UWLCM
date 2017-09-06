#pragma once
#include "Plotter2d.hpp"
#include "Plotter3d.hpp"

// 2d version
template<int NDims>
class PlotterMicro_t : public Plotter_t<NDims> 
{
  protected:
  using parent_t = Plotter_t<NDims>;

  public:
  using arr_t = typename parent_t::arr_t;

  protected:
  std::string micro;
  arr_t res;

  public:
  // cloud droplets mixing ratio
  auto h5load_rc_timestep(
    const string &file, 
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      auto snap = this->h5load_timestep(this->file, "cloud_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      res = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      auto snap = this->h5load_timestep(this->file, "rc", at);
      res = arr_t(snap);
    }
    return blitz::safeToReturn(res + 0);
  }

  // rain droplets mixing ratio
  auto h5load_rr_timestep(
    const string &file, 
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      auto snap = this->h5load_timestep(this->file, "rain_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      res = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      auto snap = this->h5load_timestep(this->file, "rc", at);
      res = arr_t(snap);
    }
    return blitz::safeToReturn(res + 0);
  }

  // activated drops mixing ratio
  auto h5load_ract_timestep(
    const string &file, 
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      auto snap = this->h5load_timestep(this->file, "actrw_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      res = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      {
        auto snap = this->h5load_timestep(this->file, "rc", at);
        res = arr_t(snap);
      }
      {
        auto snap = this->h5load_timestep(this->file, "rr", at);
        res += arr_t(snap);
      }
    }
    return blitz::safeToReturn(res + 0);
  }

  // height [m] of the center of mass of activated droplets
  double act_com_z_timestep(
    const string &file, 
    int at
  )
  {
    auto tmp = h5load_ract_timestep(this->file, at);
    arr_t ract(tmp);
    arr_t weighted(tmp);
    weighted = weighted * this->LastIndex * this->map["dz"];
    if(blitz::sum(ract) > 1e-3)
      return blitz::sum(weighted) / blitz::sum(ract);
    else
      return 0.; 
  }

  //ctor
  PlotterMicro_t(const string &file, const string &micro):
    parent_t(file),
    micro(micro),
    res(this->tmp.shape())
  {}
};


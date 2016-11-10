#pragma once
#include "Plotter2d.hpp"
#include "Plotter3d.hpp"

// 2d version
template<int NDims>
class PlotterMicro_t : public Plotter_t<NDims> 
{
  protected:
  using parent_t = Plotter_t<NDims>;
  std::string micro;

  public:
  using arr_t = typename parent_t::arr_t;

  // cloud droplets mixing ratio
  auto h5load_rc_timestep(
    const string &file, 
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      auto snap = this->h5load_timestep(this->file, "cloud_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      this->tmp = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      auto snap = this->h5load_timestep(this->file, "rc", at);
      this->tmp = arr_t(snap);
    }
    return blitz::safeToReturn(this->tmp + 0);
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
      this->tmp = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      auto snap = this->h5load_timestep(this->file, "rc", at);
      this->tmp = arr_t(snap);
    }
    return blitz::safeToReturn(this->tmp + 0);
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
      this->tmp = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      {
        auto snap = this->h5load_timestep(this->file, "rc", at);
        this->tmp = arr_t(snap);
      }
      {
        auto snap = this->h5load_timestep(this->file, "rr", at);
        this->tmp += arr_t(snap);
      }
    }
    return blitz::safeToReturn(this->tmp + 0);
  }

  //ctor
  PlotterMicro_t(const string &file, const string &micro):
    parent_t(file),
    micro(micro)
  {}
};


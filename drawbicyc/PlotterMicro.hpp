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
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      auto snap = this->h5load_timestep("cloud_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      res = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      auto snap = this->h5load_timestep("rc", at);
      res = arr_t(snap);
    }
    return blitz::safeToReturn(res + 0);
  }

  // rain droplets mixing ratio
  auto h5load_rr_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      auto snap = this->h5load_timestep("rain_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      res = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      auto snap = this->h5load_timestep("rc", at);
      res = arr_t(snap);
    }
    return blitz::safeToReturn(res + 0);
  }

  // activated drops mixing ratio
  auto h5load_ract_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      auto snap = this->h5load_timestep("actrw_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      res = arr_t(snap);
    }
    else if(this->micro == "blk_1m")
    {
      {
        auto snap = this->h5load_timestep("rc", at);
        res = arr_t(snap);
      }
      {
        auto snap = this->h5load_timestep("rr", at);
        res += arr_t(snap);
      }
    }
    return blitz::safeToReturn(res + 0);
  }

  // height [m] of the center of mass of activated droplets
  double act_com_z_timestep(int at)
  {
    auto tmp = h5load_ract_timestep(at);
    arr_t ract(tmp);
    arr_t weighted(tmp);
    weighted = weighted * this->LastIndex * this->map["dz"];
    if(blitz::sum(ract) > 1e-3)
      return blitz::sum(weighted) / blitz::sum(ract);
    else
      return 0.; 
  }

  // mean and std dev [g/kg] of the mixing ratio of activated dropelts in cloudy cells
  std::pair<double, double> cloud_ract_stats_timestep(int at)
  {
    std::pair<double, double> res;
    // read activated droplets mixing ratio 
    auto tmp = h5load_ract_timestep(at);
    arr_t ract(tmp);
    arr_t mask(tmp);
    
    mask = iscloudy_rc(mask);
    ract *= mask; // apply filter
    ract *= 1e3; // turn it into g/kg
    
    if(blitz::sum(mask) > 0.) 
      res.first = blitz::sum(ract) / blitz::sum(mask); 
    else
      res.first = 0.; 

    ract = pow(ract - res.first, 2); 
    ract *= mask; // apply filter
    if(res.first>0)
      res.second = sqrt(blitz::sum(ract) / blitz::sum(mask)); 
    else
      res.second = 0.;

    return res;
  }

  //ctor
  PlotterMicro_t(const string &file, const string &micro):
    parent_t(file),
    micro(micro),
    res(this->tmp.shape())
  {}
};


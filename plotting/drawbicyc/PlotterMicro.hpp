#pragma once
#include "Plotter2d.hpp"
#include "Plotter3d.hpp"

// TODO: make two: plotterlgrngn and plotter blk1m
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
  arr_t rhod;
  const double L_evap = 2264.76e3; // latent heat of evaporation [J/kg]

  public:
  // functions for diagnosing fields
  //
  // aerosol droplets mixing ratio
  auto h5load_ra_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      res = this->h5load_timestep("aerosol_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
    else if(this->micro == "blk_1m")
      res = 0;
    return blitz::safeToReturn(res + 0);
  }

  // cloud droplets mixing ratio
  auto h5load_rc_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      res = this->h5load_timestep("cloud_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
    else if(this->micro == "blk_1m")
      res = this->h5load_timestep("rc", at);
    return blitz::safeToReturn(res + 0);
  }

  // rain droplets mixing ratio
  auto h5load_rr_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      res = this->h5load_timestep("rain_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
    else if(this->micro == "blk_1m")
      res = this->h5load_timestep("rr", at);
    return blitz::safeToReturn(res + 0);
  }

  // activated drops mixing ratio
  auto h5load_ract_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      res = this->h5load_timestep("cloud_rw_mom3", at) * 4./3. * 3.1416 * 1e3;
      res += arr_t(this->h5load_timestep("rain_rw_mom3", at) * 4./3. * 3.1416 * 1e3);
    }
    else if(this->micro == "blk_1m")
    {
      res = this->h5load_timestep("rc", at);
      res += arr_t(this->h5load_timestep("rr", at));
    }
    return blitz::safeToReturn(res + 0);
  }

  // cloud droplets concentration [1/kg]
  auto h5load_nc_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      res = this->h5load_timestep("cloud_rw_mom0", at);
    else if(this->micro == "blk_1m")
      res = 0;
    return blitz::safeToReturn(res + 0);
  }

  // precipitation flux [W/m2]
  auto h5load_prflux_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
    {
      res = this->h5load_timestep("precip_rate", at)
              *  4./3 * 3.14 * 1e3 // to get mass
              / this->CellVol    // averaged over cell volume, TODO: make precip rate return specific moment? wouldnt need the dx and dy
              * L_evap;
    }
    else if(this->micro == "blk_1m")
      try
      {
        res = this->h5load_timestep("precip_rate", at); // precip_rate is the difference between influx and outflux
        for(int z = this->map["z"] - 2; z>=0; --z)
        {
          res(this->hrzntl_slice(z)) = res(this->hrzntl_slice(z+1)) - res(this->hrzntl_slice(z)); 
        }
        res *= rhod * this->map["dz"] * L_evap;
      }
      catch(...)
      {
        res = 0;
      }
    return blitz::safeToReturn(res + 0);
  }

  // RH
  auto h5load_RH_timestep(
    int at
  ) -> decltype(blitz::safeToReturn(arr_t() + 0))
  {
    if(this->micro == "lgrngn")
      res = this->h5load_timestep("RH", at);
    else if(this->micro == "blk_1m")
      res = 0;
    return blitz::safeToReturn(res + 0);
  }


  // functions for diagnosing statistics
  
  // helper function that calculates staistics (mean and std_dev) of a field only in cloudy cells
  std::pair<double, double> cloud_hlpr(arr_t arr, int at)
  {
    std::pair<double, double> res;
    // read activated droplets mixing ratio 
    arr_t mask(h5load_ract_timestep(at));
    mask = iscloudy_rc(mask);
    arr *= mask; // apply filter
    
    if(blitz::sum(mask) > 0.) 
      res.first = blitz::sum(arr) / blitz::sum(mask); 
    else
      res.first = 0.; 

    arr = pow(arr - res.first, 2); 
    arr *= mask; // apply filter
    if(res.first>0)
      res.second = sqrt(blitz::sum(arr) / blitz::sum(mask)); 
    else
      res.second = 0.;

    return res;
  }

  
  // height [m] of the center of mass of activated droplets
  double act_com_z_timestep(int at)
  {
    arr_t ract(h5load_ract_timestep(at));
    arr_t weighted(ract.copy());
    weighted = weighted * this->LastIndex * this->map["dz"];
    if(blitz::sum(ract) > 1e-3)
      return blitz::sum(weighted) / blitz::sum(ract);
    else
      return 0.; 
  }

  // mean and std dev [g/kg] of the mixing ratio of activated dropelts in cloudy cells (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_ract_stats_timestep(int at)
  {
    // read activated droplets mixing ratio 
    arr_t ract(h5load_ract_timestep(at));
    ract *= 1e3; // turn it into g/kg
    return cloud_hlpr(ract, at);
  }

  // mean and std_dev of concentration of activated droplets in cloudy cells [1/cm^3] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_actconc_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    // read concentration of activated droplets
    arr_t actconc(this->h5load_timestep("actrw_rw_mom0", at));
    actconc *= rhod; // b4 it was specific moment
    actconc /= 1e6; // per cm^3
    return cloud_hlpr(actconc, at);
  } 

  // mean and std_dev of supersaturation in cells with positive supersaturation [%] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> positive_supersat_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    std::pair<double, double> res;

    // read RH 
    arr_t RH(this->h5load_timestep("RH", at));
    RH -= 1.; 
    arr_t sat(RH.copy());
    sat = iscloudy_sat(RH);
    RH *= sat; //apply the cloudiness mask
    RH *= 100; // to get %
    if(blitz::sum(sat) > 0)
      res.first = blitz::sum(RH) / blitz::sum(sat); 
    else
      res.first = 0;
  
    RH = pow(RH - res.first, 2); 
    RH *= sat; // apply filter
    if(res.first>0)
      res.second = sqrt(blitz::sum(RH) / blitz::sum(sat)); 
    else
      res.second = 0.; 

    return res;
  }

  // mean and std_dev of number of SDs in cloudy cells (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_sdconc_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    arr_t sdconc(this->h5load_timestep("sd_conc", at));
    return cloud_hlpr(sdconc, at);
  }

  // mean and std_dev of number of activated SDs in cloudy cells (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_sdconc_act_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    arr_t sdconc_act(this->h5load_timestep("sd_conc_act", at));
    return cloud_hlpr(sdconc_act, at);
  }

  // mean and std_dev of mean radius of activated droplets in cloudy cells [um] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_meanr_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    // read act drop 0th raw moment / mass [1/kg]
    arr_t act0th(this->h5load_timestep("actrw_rw_mom0", at)); 
    // read act drop 1st raw moment / mass [um/kg]
    arr_t act1st(this->h5load_timestep("actrw_rw_mom1", at) * 1e6);
    // calculate mean radius, store in act1st
    act1st = where(act0th > 0, act1st / act0th, 0.);
    return cloud_hlpr(act1st, at);
  }

  // mean and std_dev of std_dev of radius of activated droplets in cloudy cells [um] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_stddevr_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    // read act drop 0th raw moment / mass [1/kg]
    arr_t act0th(this->h5load_timestep("actrw_rw_mom0", at)); 
    // read act drop 1st raw moment / mass [um/kg]
    arr_t act1st(this->h5load_timestep("actrw_rw_mom1", at) * 1e6);
    // read act drop 2nd raw moment / mass [um^2/kg]
    arr_t act2nd(this->h5load_timestep("actrw_rw_mom2", at) * 1e12);
    // calculate stddev of radius, store in act1st
    act1st = where(act0th > 0, 
      act2nd / act0th - act1st / act0th * act1st / act0th, 0.);
    // might be slightly negative due to numerical errors
    act1st = where(act1st < 0, 0, act1st);
    act1st = sqrt(act1st);
    return cloud_hlpr(act1st, at);
  }

  // surface precipitation since last output [mm/day]
  double calc_surf_precip(double prec_vol_diff)
  {
    if(this->micro == "lgrngn")
      return prec_vol_diff / this->DomainSurf / (double(this->map["outfreq"]) * this->map["dt"] / 3600. / 24.) * 1e3; // SDM
    if(this->micro == "blk_1m")
      return prec_vol_diff / double(this->map["outfreq"]) // flux in [kg / m^3 / s] averaged over time since last output and over cells on the bottom
                     / (this->map["x"] * this->map["y"])
                     * 3600. * 24. // per day
                     * this->map["dz"]     // per m^2
                     / 1e3         // to m^3 of water
                     * 1e3;        // to mm
  }

  // accumulated surface precipitation [mm]
  double calc_acc_surf_precip(double prec_vol)
  {
    if(this->micro == "lgrngn")
      return prec_vol / this->DomainSurf * 1e3; 
    if(this->micro == "blk_1m")
      return prec_vol * this->map["dt"] 
                     / (this->map["x"] * this->map["y"])
                     * this->map["dz"]     // per m^2
                     / 1e3         // to m^3 of water
                     * 1e3;        // to mm
  }

  //ctor
  PlotterMicro_t(const string &file, const string &micro):
    parent_t(file),
    micro(micro),
    res(this->tmp.shape()),
    rhod(this->h5load(file + "/const.h5", "G"))
  {
  }
};


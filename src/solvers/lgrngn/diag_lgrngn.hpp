#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::diag()
{
  parent_t::diag();

  // recording super-droplet concentration per grid cell 
  prtcls->diag_all();
  prtcls->diag_sd_conc();
  this->record_aux("sd_conc", prtcls->outbuf());

  // recording concentration of SDs that represent activated droplets 
  /*
  prtcls->diag_rw_ge_rc();
  prtcls->diag_sd_conc();
  this->record_aux("sd_conc_act", prtcls->outbuf());
*/

  // recording relative humidity
  prtcls->diag_RH();
  this->record_aux("RH", prtcls->outbuf());

  // recording pressure
  /*
  prtcls->diag_pressure();
  this->record_aux("libcloud_pressure", prtcls->outbuf());

  // recording temperature
  prtcls->diag_temperature();
  this->record_aux("libcloud_temperature", prtcls->outbuf());
  */

  // recording precipitation rate per grid cel
  prtcls->diag_water();
  prtcls->diag_precip_rate();
  this->record_aux("precip_rate", prtcls->outbuf());
  prtcls->diag_ice();
  prtcls->diag_precip_rate_ice();
  this->record_aux("precip_rate_ice_mass", prtcls->outbuf());

  // recording 0th mom of rw of rd>=0.8um
//  prtcls->diag_dry_rng(0.7999e-6, 1);
//  prtcls->diag_wet_mom(0);
//  this->record_aux("rd_geq_0.8um_rw_mom0", prtcls->outbuf());
//
//  // recording 0th mom of rw of rd>=0.8um
//  prtcls->diag_dry_rng(0, 0.8e-6);
//  prtcls->diag_wet_mom(0);
//  this->record_aux("rd_lt_0.8um_rw_mom0", prtcls->outbuf());

//    // recording 1st mom of rw of gccns
//    prtcls->diag_dry_rng(2e-6, 1);
//    prtcls->diag_wet_mom(1);
//    this->record_aux("gccn_rw_mom1", prtcls->outbuf());
//
//    // recording 0th mom of rw of gccns
//    prtcls->diag_dry_rng(2e-6, 1);
//    prtcls->diag_wet_mom(0);
//    this->record_aux("gccn_rw_mom0", prtcls->outbuf());
//
//    // recording 1st mom of rw of non-gccns
//    prtcls->diag_dry_rng(0., 2e-6);
//    prtcls->diag_wet_mom(1);
//    this->record_aux("non_gccn_rw_mom1", prtcls->outbuf());
//
//    // recording 0th mom of rw of gccns
//    prtcls->diag_dry_rng(0., 2e-6);
//    prtcls->diag_wet_mom(0);
//    this->record_aux("non_gccn_rw_mom0", prtcls->outbuf());

  prtcls->diag_water();

  // recording 0th mom of rw of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_wet_mom(0);
  this->record_aux("actrw_rw_mom0", prtcls->outbuf());

  // recording 1st mom of rw of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_wet_mom(1);
  this->record_aux("actrw_rw_mom1", prtcls->outbuf());

  // recording 2nd mom of rw of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_wet_mom(2);
  this->record_aux("actrw_rw_mom2", prtcls->outbuf());

  // recording 3rd mom of rw of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_wet_mom(3);
  this->record_aux("actrw_rw_mom3", prtcls->outbuf());
/*
  // recording 1st mom of rd of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_dry_mom(1);
  this->record_aux("actrw_rd_mom1", prtcls->outbuf());

  // recording 0th mom of rd of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_dry_mom(0);
  this->record_aux("actrw_rd_mom0", prtcls->outbuf());
  
  // recording 1st mom of rd of activated drops
  prtcls->diag_RH_ge_Sc();
  prtcls->diag_dry_mom(1);
  this->record_aux("actRH_rd_mom1", prtcls->outbuf());
 
  // recording 3rd mom of rw of activated drops
  prtcls->diag_RH_ge_Sc();
  prtcls->diag_wet_mom(3);
  this->record_aux("actRH_rw_mom3", prtcls->outbuf());

  // recording 0th mom of rd of activated drops
  prtcls->diag_RH_ge_Sc();
  prtcls->diag_dry_mom(0);
  this->record_aux("actRH_rd_mom0", prtcls->outbuf());
  */

  // recording 0th wet mom of radius of rain drops (r>25um)
  prtcls->diag_wet_rng(25.e-6, 1);
  prtcls->diag_wet_mom(0);
  this->record_aux("rain_rw_mom0", prtcls->outbuf());
  
  // recording 3rd wet mom of radius of rain drops (r>25um)
  prtcls->diag_wet_rng(25.e-6, 1);
  prtcls->diag_wet_mom(3);
  this->record_aux("rain_rw_mom3", prtcls->outbuf());

  // recording 0th wet mom of radius of cloud drops (.5um< r < 25um)
  prtcls->diag_wet_rng(.5e-6, 25.e-6);
  prtcls->diag_wet_mom(0);
  this->record_aux("cloud_rw_mom0", prtcls->outbuf());

  // recording 3rd wet mom of radius of cloud drops (.5um< r < 25um)
  prtcls->diag_wet_rng(.5e-6, 25.e-6);
  prtcls->diag_wet_mom(3);
  this->record_aux("cloud_rw_mom3", prtcls->outbuf());

  // recording 0th wet mom of radius of aerosols (r < .5um)
//  prtcls->diag_wet_rng(0., .5e-6);
//  prtcls->diag_wet_mom(0);
//  this->record_aux("aerosol_rw_mom0", prtcls->outbuf());
//
//  // recording 3rd wet mom of radius of aerosols (r < .5um)
//  prtcls->diag_wet_rng(0., .5e-6);
//  prtcls->diag_wet_mom(3);
//  this->record_aux("aerosol_rw_mom3", prtcls->outbuf());
//
//  // recording 1st wet mom of radius of all particles
//  prtcls->diag_all();
//  prtcls->diag_wet_mom(1);
//  this->record_aux("all_rw_mom1", prtcls->outbuf());
//
//  // recording 2nd wet mom of radius of all particles
//  prtcls->diag_all();
//  prtcls->diag_wet_mom(2);
//  this->record_aux("all_rw_mom2", prtcls->outbuf());
//
//  // recording 6th wet mom of radius of all particles
//  prtcls->diag_all();
//  prtcls->diag_wet_mom(6);
//  this->record_aux("all_rw_mom6", prtcls->outbuf());

/*    
    // recording divergence of the velocity field
    prtcls->diag_vel_div();
    this->record_aux("vel_div", prtcls->outbuf());

    // record 0th moment of cloud droplets with kappa < 0.5
    prtcls->diag_kappa_rng(0.0,0.5);
    prtcls->diag_dry_mom(0);
    this->record_aux("smallkappa_rd_mom0", prtcls->outbuf());    

    // record 0th moment of cloud droplets with kappa >= 0.5
    prtcls->diag_kappa_rng(0.5,2.0);
    prtcls->diag_dry_mom(0);
    this->record_aux("bigkappa_rd_mom0", prtcls->outbuf());    
*/

//  // recording 0th wet mom of radius of big rain drops (r>40um)
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_wet_mom(0);
//  this->record_aux("bigrain_rw_mom0", prtcls->outbuf());
//
//  // recording 1st mom of incloud_time of big rain drops (r>40um)
//  if(params.cloudph_opts_init.diag_incloud_time)
//  {
//    prtcls->diag_wet_rng(40.e-6, 1);
//    prtcls->diag_incloud_time_mom(1);
//    this->record_aux("bigrain_incloud_time_mom1", prtcls->outbuf());
//  }
//
//  // recording 1st mom of kappa of big rain drops (r>40um)
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_kappa_mom(1);
//  this->record_aux("bigrain_kappa_mom1", prtcls->outbuf());
//
//  // recording 1st mom of rd of big rain drops (r>40um)
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_dry_mom(1);
//  this->record_aux("bigrain_rd_mom1", prtcls->outbuf());
//
//  // recording 0th mom of rw of big rain drops (r>40um) with kappa > 0.61
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_kappa_rng_cons(0.61000001, 10);
//  prtcls->diag_wet_mom(0);
//  this->record_aux("bigrain_gccn_rw_mom0", prtcls->outbuf());

  if (params.cloudph_opts_init.ice_switch)
  {
    prtcls->diag_ice();

    // recording ice mixing ratio
    prtcls->diag_ice_mass();
    this->record_aux("r_i", prtcls->outbuf());

    // recording 0th moment of ice
    prtcls->diag_ice_a_mom(0);
    this->record_aux("ice_mom0", prtcls->outbuf());

    // recording 1st moment of ice_a
    prtcls->diag_ice_a_mom(1);
    this->record_aux("ice_a_mom1", prtcls->outbuf());

    // recording 1st moment of ice_c
    prtcls->diag_ice_c_mom(1);
    this->record_aux("ice_c_mom1", prtcls->outbuf());
  }
  // recording requested statistical moments
  if ((this->timestep ) % static_cast<int>(params.outfreq_spec) == 0)
  {
    // dry
    int rng_num = 0;
    for (auto &rng_moms : params.out_dry)
    {
      auto &rng(rng_moms.first);
      prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
      for (auto &mom : rng_moms.second)
      {
        prtcls->diag_dry_mom(mom);
        this->record_aux(aux_name("rd", rng_num, mom), prtcls->outbuf());
      }
      rng_num++;
    }

    // wet
    prtcls->diag_water();
    rng_num = 0;
    for (auto &rng_moms : params.out_wet)
    {
      auto &rng(rng_moms.first);
      prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
      for (auto &mom : rng_moms.second)
      {
        prtcls->diag_wet_mom(mom);
        this->record_aux(aux_name("rw", rng_num, mom), prtcls->outbuf());
      }
      rng_num++;
    }

    if (params.cloudph_opts_init.ice_switch)
    {
      prtcls->diag_ice();
      // ice_a
      rng_num = 0;
      for (auto &rng_moms : params.out_ice)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_ice_a_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_ice_a_mom(mom);
          this->record_aux(aux_name("ice_a", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
      // ice_c
      rng_num = 0;
      for (auto &rng_moms : params.out_ice)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_ice_c_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_ice_c_mom(mom);
          this->record_aux(aux_name("ice_c", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
  }

  diag_prev_step = 1;
} 

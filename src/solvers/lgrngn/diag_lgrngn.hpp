#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::diag()
{
  parent_t::diag();

  // ---- DEBUGGING: record interpolated refined courants ----
  // we dont have a function for recording variables are at the edges, so we crop last courant and shift values to the center
  // making a copy is not efficient, but this is just for debugging
  // also wont work with MPI

  {
    rng_t rng_m1(this->mem->grid_size_ref[0].first(), this->mem->grid_size_ref[0].last());
    typename parent_t::arr_t contiguous_arr(rng_m1, this->mem->grid_size_ref[1], this->mem->grid_size_ref[2]);
    contiguous_arr = this->courants[0](rng_m1, blitz::Range::all(), blitz::Range::all());
    this->record_aux_refined("courants[0]", contiguous_arr.data()); 
  }
  {
    rng_t rng_m1(this->mem->grid_size_ref[1].first(), this->mem->grid_size_ref[1].last());
    typename parent_t::arr_t contiguous_arr(this->mem->grid_size_ref[0], rng_m1, this->mem->grid_size_ref[2]);
    contiguous_arr = this->courants[1](blitz::Range::all(), rng_m1, blitz::Range::all());
    this->record_aux_refined("courants[1]", contiguous_arr.data()); 
  }
  {
    rng_t rng_m1(this->mem->grid_size_ref[2].first(), this->mem->grid_size_ref[2].last());
    typename parent_t::arr_t contiguous_arr(this->mem->grid_size_ref[0], this->mem->grid_size_ref[1], rng_m1);
    contiguous_arr = this->courants[2](blitz::Range::all(), blitz::Range::all(), rng_m1);
    this->record_aux_refined("courants[2]", contiguous_arr.data()); 

//    std::cerr << "courants[2](blitz::Range::all(), blitz::Range::all(), rng_m1): " << this->courants[2](blitz::Range::all(), blitz::Range::all(), rng_m1) << std::endl;
  }

//  std::cerr << "courants[2]: " << this->courants[2] << std::endl;

  this->record_aux_dsc_refined("refined u", this->mem->refinee(this->ix_r2r.at(ix::u)));
  this->record_aux_dsc_refined("refined v", this->mem->refinee(this->ix_r2r.at(ix::v)));
  this->record_aux_dsc_refined("refined w", this->mem->refinee(this->ix_r2r.at(ix::w)));

  // ---- END OF COURANT DEBUGGING ----


  // recording super-droplet concentration per grid cell 
  prtcls->diag_all();
  prtcls->diag_sd_conc();
  this->record_aux_refined("sd_conc", prtcls->outbuf());

  // recording concentration of SDs that represent activated droplets 
  /*
  prtcls->diag_rw_ge_rc();
  prtcls->diag_sd_conc();
  this->record_aux_refined("sd_conc_act", prtcls->outbuf());
*/

  // recording relative humidity
  prtcls->diag_RH();
  this->record_aux_refined("RH", prtcls->outbuf());

  // recording pressure
  /*
  prtcls->diag_pressure();
  this->record_aux_refined("libcloud_pressure", prtcls->outbuf());

  // recording temperature
  prtcls->diag_temperature();
  this->record_aux_refined("libcloud_temperature", prtcls->outbuf());
  */

  // recording precipitation rate per grid cel
  prtcls->diag_all();
  prtcls->diag_precip_rate();
  this->record_aux_refined("precip_rate", prtcls->outbuf());

//  // recording 0th mom of rw of rd>=0.8um
//  prtcls->diag_dry_rng(0.7999e-6, 1);
//  prtcls->diag_wet_mom(0);
//  this->record_aux_refined("rd_geq_0.8um_rw_mom0", prtcls->outbuf());
//
//  // recording 0th mom of rw of rd>=0.8um
//  prtcls->diag_dry_rng(0, 0.8e-6);
//  prtcls->diag_wet_mom(0);
//  this->record_aux_refined("rd_lt_0.8um_rw_mom0", prtcls->outbuf());

//    // recording 1st mom of rw of gccns
//    prtcls->diag_dry_rng(2e-6, 1);
//    prtcls->diag_wet_mom(1);
//    this->record_aux_refined("gccn_rw_mom1", prtcls->outbuf());
//
//    // recording 0th mom of rw of gccns
//    prtcls->diag_dry_rng(2e-6, 1);
//    prtcls->diag_wet_mom(0);
//    this->record_aux_refined("gccn_rw_mom0", prtcls->outbuf());
//
//    // recording 1st mom of rw of non-gccns
//    prtcls->diag_dry_rng(0., 2e-6);
//    prtcls->diag_wet_mom(1);
//    this->record_aux_refined("non_gccn_rw_mom1", prtcls->outbuf());
//
//    // recording 0th mom of rw of gccns
//    prtcls->diag_dry_rng(0., 2e-6);
//    prtcls->diag_wet_mom(0);
//    this->record_aux_refined("non_gccn_rw_mom0", prtcls->outbuf());

  // recording 0th mom of rw of activated drops
//  prtcls->diag_rw_ge_rc();
//  prtcls->diag_wet_mom(0);
//  this->record_aux_refined("actrw_rw_mom0", prtcls->outbuf());
//
//  // recording 1st mom of rw of activated drops
//  prtcls->diag_rw_ge_rc();
//  prtcls->diag_wet_mom(1);
//  this->record_aux_refined("actrw_rw_mom1", prtcls->outbuf());
//
//  // recording 2nd mom of rw of activated drops
//  prtcls->diag_rw_ge_rc();
//  prtcls->diag_wet_mom(2);
//  this->record_aux_refined("actrw_rw_mom2", prtcls->outbuf());
//
//  // recording 3rd mom of rw of activated drops
//  prtcls->diag_rw_ge_rc();
//  prtcls->diag_wet_mom(3);
//  this->record_aux_refined("actrw_rw_mom3", prtcls->outbuf());
/*
  // recording 1st mom of rd of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_dry_mom(1);
  this->record_aux_refined("actrw_rd_mom1", prtcls->outbuf());

  // recording 0th mom of rd of activated drops
  prtcls->diag_rw_ge_rc();
  prtcls->diag_dry_mom(0);
  this->record_aux_refined("actrw_rd_mom0", prtcls->outbuf());
  
  // recording 1st mom of rd of activated drops
  prtcls->diag_RH_ge_Sc();
  prtcls->diag_dry_mom(1);
  this->record_aux_refined("actRH_rd_mom1", prtcls->outbuf());
 
  // recording 3rd mom of rw of activated drops
  prtcls->diag_RH_ge_Sc();
  prtcls->diag_wet_mom(3);
  this->record_aux_refined("actRH_rw_mom3", prtcls->outbuf());

  // recording 0th mom of rd of activated drops
  prtcls->diag_RH_ge_Sc();
  prtcls->diag_dry_mom(0);
  this->record_aux_refined("actRH_rd_mom0", prtcls->outbuf());
  */

  // recording 0th wet mom of radius of rain drops (r>25um)
  prtcls->diag_wet_rng(25.e-6, 1);
  prtcls->diag_wet_mom(0);
  this->record_aux_refined("rain_rw_mom0", prtcls->outbuf());
  
  // recording 3rd wet mom of radius of rain drops (r>25um)
  prtcls->diag_wet_rng(25.e-6, 1);
  prtcls->diag_wet_mom(3);
  this->record_aux_refined("rain_rw_mom3", prtcls->outbuf());

  // recording 0th wet mom of radius of cloud drops (.5um< r < 25um)
  prtcls->diag_wet_rng(.5e-6, 25.e-6);
  prtcls->diag_wet_mom(0);
  this->record_aux_refined("cloud_rw_mom0", prtcls->outbuf());

  // recording 1th wet mom of radius of cloud drops (.5um< r < 25um)
  prtcls->diag_wet_rng(.5e-6, 25.e-6);
  prtcls->diag_wet_mom(1);
  this->record_aux_refined("cloud_rw_mom1", prtcls->outbuf());

  // recording 2th wet mom of radius of cloud drops (.5um< r < 25um)
  prtcls->diag_wet_rng(.5e-6, 25.e-6);
  prtcls->diag_wet_mom(2);
  this->record_aux_refined("cloud_rw_mom2", prtcls->outbuf());

  // recording 3rd wet mom of radius of cloud drops (.5um< r < 25um)
  prtcls->diag_wet_rng(.5e-6, 25.e-6);
  prtcls->diag_wet_mom(3);
  this->record_aux_refined("cloud_rw_mom3", prtcls->outbuf());

//    // recording 0th wet mom of radius of aerosols (r < .5um)
//    prtcls->diag_wet_rng(0., .5e-6);
//    prtcls->diag_wet_mom(0);
//    this->record_aux_refined("aerosol_rw_mom0", prtcls->outbuf());

//    // recording 3rd wet mom of radius of aerosols (r < .5um)
//    prtcls->diag_wet_rng(0., .5e-6);
//    prtcls->diag_wet_mom(3);
//    this->record_aux_refined("aerosol_rw_mom3", prtcls->outbuf());

/*    
    // recording divergence of the velocity field
    prtcls->diag_vel_div();
    this->record_aux_refined("vel_div", prtcls->outbuf());

    // record 0th moment of cloud droplets with kappa < 0.5
    prtcls->diag_kappa_rng(0.0,0.5);
    prtcls->diag_dry_mom(0);
    this->record_aux_refined("smallkappa_rd_mom0", prtcls->outbuf());    

    // record 0th moment of cloud droplets with kappa >= 0.5
    prtcls->diag_kappa_rng(0.5,2.0);
    prtcls->diag_dry_mom(0);
    this->record_aux_refined("bigkappa_rd_mom0", prtcls->outbuf());    
*/

//  // recording 0th wet mom of radius of big rain drops (r>40um)
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_wet_mom(0);
//  this->record_aux_refined("bigrain_rw_mom0", prtcls->outbuf());
//
//  // recording 1st mom of incloud_time of big rain drops (r>40um)
//  if(params.cloudph_opts_init.diag_incloud_time)
//  {
//    prtcls->diag_wet_rng(40.e-6, 1);
//    prtcls->diag_incloud_time_mom(1);
//    this->record_aux_refined("bigrain_incloud_time_mom1", prtcls->outbuf());
//  }
//
//  // recording 1st mom of kappa of big rain drops (r>40um)
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_kappa_mom(1);
//  this->record_aux_refined("bigrain_kappa_mom1", prtcls->outbuf());
//
//  // recording 1st mom of rd of big rain drops (r>40um)
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_dry_mom(1);
//  this->record_aux_refined("bigrain_rd_mom1", prtcls->outbuf());
//
//  // recording 0th mom of rw of big rain drops (r>40um) with kappa > 0.61
//  prtcls->diag_wet_rng(40.e-6, 1);
//  prtcls->diag_kappa_rng_cons(0.61000001, 10);
//  prtcls->diag_wet_mom(0);
//  this->record_aux_refined("bigrain_gccn_rw_mom0", prtcls->outbuf());

  // recording requested statistical moments
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
        this->record_aux_refined(aux_name("rd", rng_num, mom), prtcls->outbuf());
      }
      rng_num++;
    }
  }
  {
    // wet
    int rng_num = 0;
    for (auto &rng_moms : params.out_wet)
    {
      auto &rng(rng_moms.first);
      prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
      for (auto &mom : rng_moms.second)
      {
        prtcls->diag_wet_mom(mom);
        this->record_aux_refined(aux_name("rw", rng_num, mom), prtcls->outbuf());
      }
      rng_num++;
    }
  }
} 

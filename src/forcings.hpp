#pragma once

//TODO: make these functions return arrays
//

// Grabowski & Smolarkiewicz 1995 "two-time semi-lagrangian modeling of precipitating clouds" eq. (2)
template <class ct_params_t>
void slvr_common<ct_params_t>::buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv)
{
  const auto &ijk = this->ijk;

  namespace moist_air = libcloudphxx::common::moist_air;
  const real_t eps = moist_air::R_v<real_t>() / moist_air::R_d<real_t>() - 1.;
  if(params.buoyancy_wet)
    F(ijk).reindex(this->zero) = 
      (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared) * (
        (th(ijk).reindex(this->zero) - (*params.th_e)(this->vert_idx)) / (*params.th_ref)(this->vert_idx)
        + eps * (rv(ijk).reindex(this->zero) - (*params.rv_e)(this->vert_idx)) 
        - (r_l(ijk).reindex(this->zero) - (*params.rl_e)(this->vert_idx))
      );
  else
    F(ijk).reindex(this->zero) = 
      (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared) * (
        (th(ijk).reindex(this->zero) - (*params.th_e)(this->vert_idx)) / (*params.th_ref)(this->vert_idx)
      );
}

template <class ct_params_t>
void slvr_common<ct_params_t>::radiation(typename parent_t::arr_t &rv)
// calc upward radiative flux through the bottom of the cells
{
  const auto &ijk = this->ijk;
  if(params.radiation)
  {
    namespace idxperm = libmpdataxx::idxperm;
    using ix = typename ct_params_t::ix;
    constexpr int perm_no=ix::w; // 1 for 2D, 2 for 3D
  
    int nz = this->mem->grid_size[perm_no].length(); 
  
    // index of first cell above inversion
    tmp1(ijk)  = rv(ijk) + r_l(ijk);
  
    k_i.reindex(this->zero_plane) = blitz::first( tmp1(ijk).reindex(this->zero) < params.ForceParameters.q_i, this->vert_idx); // vertical index of first cell above inversion (inversion is at the lower edge of this cell)
  
    // calc Eqs. 5 and 6 from Ackerman et al 2009
    // calc sum of r_l above certain level and store it in tmp1
    tmp1(ijk) = r_l(ijk);
  
    tmp1(ijk).reindex(this->zero) *= -params.dz * params.ForceParameters.heating_kappa * (*params.rhod)(this->vert_idx);
  
    for(int z = nz-2 ; z >= 0; --z)
      tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) += tmp1(idxperm::pi<perm_no>(z+1, this->hrzntl_subdomain));
  
    auto ground = idxperm::pi<perm_no>(0, this->hrzntl_subdomain);
    auto noground = idxperm::pi<perm_no>(rng_t(1, nz-1), this->hrzntl_subdomain);
    auto notop = idxperm::pi<perm_no>(rng_t(0, nz-2), this->hrzntl_subdomain);
  
    radiative_flux(ijk) = params.ForceParameters.F_0 * exp(tmp1(ijk)); 
  
    // calc sum of r_l below certain level and store it in tmp1
    tmp1(ijk) = r_l(ijk);
  
    tmp1(ijk).reindex(this->zero) *= - params.dz * params.ForceParameters.heating_kappa * (*params.rhod)(this->vert_idx);
  
    // copy one cell upwards
    for(int z = nz-1 ; z >= 1; --z)
      tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) = tmp1(idxperm::pi<perm_no>(z-1, this->hrzntl_subdomain));
    tmp1(ground) = 0.;
  
    for(int z = 1 ; z <= nz-1; ++z)
      tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) += tmp1(idxperm::pi<perm_no>(z-1, this->hrzntl_subdomain));

    radiative_flux(ijk) += params.ForceParameters.F_1 * exp(tmp1(ijk));
  
    // free atmosphere part
    radiative_flux(ijk).reindex(this->zero) += where(this->vert_idx > k_i.reindex(this->zero_plane)(blitz::tensor::i, blitz::tensor::j),  // works even in 2D ?!?!
        (libcloudphxx::common::moist_air::c_pd<setup::real_t>() / si::joules * si::kilograms * si::kelvins) * params.ForceParameters.rho_i * params.ForceParameters.D *
        (0.25 * pow((this->vert_idx - 0.5) * params.dz - (k_i.reindex(this->zero_plane)(blitz::tensor::i, blitz::tensor::j) - .5) * params.dz, 4./3) +
        (k_i.reindex(this->zero_plane)(blitz::tensor::i, blitz::tensor::j) - .5) * params.dz * pow((this->vert_idx - 0.5) * params.dz - (k_i.reindex(this->zero_plane)(blitz::tensor::i, blitz::tensor::j) - .5) * params.dz, 1./3))
      , 0);
  }
  else
    radiative_flux(ijk)=0.;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens_impl(iles_tag)
{
  params.update_surf_flux_sens(surf_flux_sens(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  //F(ijk).reindex(this->zero) = surf_flux_sens(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j) 
  //                             * (*params.hgt_fctr_sclr)(this->vert_idx);
  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_sens(this->hrzntl_slice(0)) * (*params.hgt_fctr_sclr)(k);
  }

//  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  //this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens_impl(smg_tag)
{
  params.update_surf_flux_sens(surf_flux_sens(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  F(this->ijk) = 0;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens()
{
  surf_sens_impl(sgs_tag{});
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent_impl(iles_tag)
{
  params.update_surf_flux_lat(surf_flux_lat(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  //F(ijk).reindex(this->zero) = surf_flux_lat(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j)  
  //                             * (*params.hgt_fctr_sclr)(this->vert_idx);
  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_lat(this->hrzntl_slice(0)) * (*params.hgt_fctr_sclr)(k);
  }

//  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  //this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent_impl(smg_tag)
{
  params.update_surf_flux_lat(surf_flux_lat(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  F(this->ijk) = 0;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent()
{
  surf_latent_impl(sgs_tag{});
}

template <class ct_params_t>
void slvr_common<ct_params_t>::subsidence(const int &type) // large-scale vertical wind
{
  const auto &ijk = this->ijk;
  if(params.subsidence)
  {
    tmp1(ijk) = this->state(type)(ijk);
    this->vert_grad_cnt(tmp1, F, params.dz);
    F(ijk).reindex(this->zero) *= - (*params.w_LS)(this->vert_idx);

//    tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  //  this->smooth(tmp1, F);
  }
  else
    F(ijk)=0.;
}

// Coriolis force including large-scale geostrophic wind
// F_u = + coriolis_parameter * (v - v_geostrophic)
// F_v = - coriolis_parameter * (u - u_geostrophic)
// the +- sign is handled by calc_forces
template <class ct_params_t>
void slvr_common<ct_params_t>::coriolis(
  const int &vel_idx
) 
{
  const auto &ijk = this->ijk;
  if(params.coriolis && ct_params_t::n_dims==3) // TODO: n_dims is known at compile time
  {
    F(ijk).reindex(this->zero) = params.ForceParameters.coriolis_parameter *
      (this->state(vel_idx)(ijk).reindex(this->zero) - (*params.geostr[vel_idx])(this->vert_idx));
  }
  else
    F(ijk)=0.;
}

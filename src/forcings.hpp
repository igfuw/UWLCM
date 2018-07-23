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
    tmp1(ijk).reindex(this->zero) = 
      (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared) * (
        (th(ijk).reindex(this->zero) - (*params.th_e)(this->vert_idx)) / (*params.th_ref)(this->vert_idx)
        + eps * (rv(ijk).reindex(this->zero) - (*params.rv_e)(this->vert_idx)) 
        - r_l(ijk).reindex(this->zero)
      );
  else
    tmp1(ijk).reindex(this->zero) = 
      (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared) * (
        (th(ijk).reindex(this->zero) - (*params.th_e)(this->vert_idx)) / (*params.th_ref)(this->vert_idx)
      );

//  this->smooth(tmp1, F);
  F(ijk) = tmp1(ijk);
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
  
    k_i(this->hrzntl_subdomain) = blitz::first( tmp1(ijk).reindex(this->zero) < params.ForceParameters.q_i, this->vert_idx); // vertical index of first cell above inversion (inversion is at the lower edge of this cell)
  
    // calc Eqs. 5 and 6 from Ackerman et al 2009
    // calc sum of r_l above certain level and store it in tmp1
    tmp1(ijk) = r_l(ijk);
  
    tmp1(ijk).reindex(this->zero) *= params.ForceParameters.heating_kappa * (*params.rhod)(this->vert_idx);
  
    for(int z = nz-2 ; z >= 0; --z)
      tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) += tmp1(idxperm::pi<perm_no>(z+1, this->hrzntl_subdomain));
  
    auto ground = idxperm::pi<perm_no>(0, this->hrzntl_subdomain);
    auto noground = idxperm::pi<perm_no>(rng_t(1, nz-1), this->hrzntl_subdomain);
    auto notop = idxperm::pi<perm_no>(rng_t(0, nz-2), this->hrzntl_subdomain);
  
    // multiply by distance from the bottom of the cell to the top of the domain
    tmp1(ground) *= - (nz - 1) * params.dz * tmp1(ground);
    tmp1(noground).reindex(this->zero) *= - (nz - this->vert_idx - 1.5) * params.dz;  // vert_idx starts from 0, but its 2nd cell from ground; do not merge this line with F_0 * exp(...) since it gives some strange values!
  
    F(ijk) = params.ForceParameters.F_0 * exp(tmp1(ijk)); 
  
    // calc sum of r_l below certain level and store it in tmp1
    tmp1(ijk) = r_l(ijk);
  
    tmp1(ijk).reindex(this->zero) *= params.ForceParameters.heating_kappa * (*params.rhod)(this->vert_idx);
  
    // copy one cell upwards
    for(int z = nz-1 ; z >= 1; --z)
      tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) = tmp1(idxperm::pi<perm_no>(z-1, this->hrzntl_subdomain));
    tmp1(ground) = 0.;
  
    for(int z = 1 ; z <= nz-1; ++z)
      tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) += tmp1(idxperm::pi<perm_no>(z-1, this->hrzntl_subdomain));
  
    // multiply by distance from the bottom of the cell to the bottom of the domain
    tmp1(noground).reindex(this->zero) *= - (this->vert_idx + 0.5) * params.dz;
    tmp1(ground) = 0.;
    F(ijk) += params.ForceParameters.F_1 * exp(tmp1(ijk));
  
    // free atmosphere part
    F(ijk).reindex(this->zero) += where(this->vert_idx > k_i(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j),  // works even in 2D ?!?!
        (libcloudphxx::common::moist_air::c_pd<setup::real_t>() / si::joules * si::kilograms * si::kelvins) * params.ForceParameters.rho_i * params.ForceParameters.D *
        (0.25 * pow((this->vert_idx - 0.5) * params.dz - (k_i(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j) - .5) * params.dz, 4./3) +
        (k_i(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j) - .5) * params.dz * pow((this->vert_idx - 0.5) * params.dz - (k_i(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j) - .5) * params.dz, 1./3))
        , 0);
  //  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  //  this->smooth(tmp1, F);
  }
  else
    F(ijk)=0.;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens()
{
  const auto &ijk = this->ijk;
  //TODO: each thread has surf_flux_sens of the size of the domain of all threads and each updates all of it
  //      either make it shared among threads and updated by one all make it of the size of hrzntl_subdomain
  params.update_surf_flux_sens(surf_flux_sens, this->timestep, this->dt);
  F(ijk).reindex(this->zero) = surf_flux_sens(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j) 
                               * (*params.hgt_fctr_sclr)(this->vert_idx);

//  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  //this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent()
{
  const auto &ijk = this->ijk;
  //TODO: each thread has surf_flux_sens of the size of the domain of all threads and each updates all of it
  //      either make it shared among threads and updated by one all make it of the size of hrzntl_subdomain
  params.update_surf_flux_lat(surf_flux_lat, this->timestep, this->dt);
  F(ijk).reindex(this->zero) = surf_flux_lat(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j)  
                               * (*params.hgt_fctr_sclr)(this->vert_idx);

//  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  //this->smooth(tmp1, F);
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

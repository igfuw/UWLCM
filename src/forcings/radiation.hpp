#pragma once

//TODO: make these functions return arrays

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

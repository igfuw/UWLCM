#pragma once

//TODO: make these functions return arrays

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv)
{
  const auto &ijk = this->ijk;
  const real_t g = 9.81; 

  namespace moist_air = libcloudphxx::common::moist_air;
  const real_t eps = moist_air::R_v<real_t>() / moist_air::R_d<real_t>() - 1.;
  tmp1(ijk).reindex(this->zero) = g * ((th(ijk).reindex(this->zero) - (*params.th_e)(this->vert_idx)) / (*params.th_ref)(this->vert_idx)) + eps * (rv(ijk).reindex(this->zero) - (*params.rv_e)(this->vert_idx)) - r_l(ijk).reindex(this->zero);

  this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::radiation(typename parent_t::arr_t &rv)
// calc upward radiative flux through the bottom of the cells
{
  const auto &ijk = this->ijk;
  int nz = this->mem->grid_size[this->vert_dim].length(); 
  namespace idxperm = libmpdataxx::idxperm;

  // index of first cell above inversion
  tmp1(ijk)  = rv(ijk) + r_l(ijk);
  k_i(this->hrzntl_subdomain) = blitz::first( tmp1 < setup::q_i, this->vert_idx); 

  // calc Eqs. 5 and 6 from Ackerman et al 2009
  const int perm_no = this->vert_dim;
  // calc sum of r_l above certain level and store it in tmp1
  tmp1(ijk) = r_l(ijk);
  for(int z = nz-2 ; z >= 0; --z)
    tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) += tmp1(idxperm::pi<perm_no>(z+1, this->hrzntl_subdomain));
  auto ground = idxperm::pi<perm_no>(0, this->hrzntl_subdomain);
  auto noground = idxperm::pi<perm_no>(rng_t(1, nz-1), this->hrzntl_subdomain);
  auto notop = idxperm::pi<perm_no>(rng_t(0, nz-2), this->hrzntl_subdomain);
  F(noground) = setup::F_0 * exp(- (nz - z - 0.5) * params.dz * tmp1(noground)); 
  F(ground) = setup::F_0 * exp(- (nz - z - 1) * params.dz * tmp1(ground));

  // calc sum of r_l below certain level and store it in tmp1
  tmp1 = r_l;
  for(int z = 1 ; z < nz-1; ++z)
    tmp1(idxperm::pi<perm_no>(z, this->hrzntl_subdomain)) += tmp1(idxperm::pi<perm_no>(z-1, this->hrzntl_subdomain));
  F(noground) += setup::F_1 * exp(- (z - 0.5) * params.dz * tmp1(notop));

  F(ijk).reindex(this->zero) += where(this->vert_idx > k_i(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j),  // works even in 2D ?!?!
      setup::c_p * setup::rho_i * setup::D *
      (0.25 * pow((z - 0.5) * params.dz - (k_i(this->hrzntl_subdomain) - .5) * params.dz, 4./3) +
      (k_i(this->hrzntl_subdomain) - .5) * params.dz * pow((z - 0.5) * params.dz - (k_i(this->hrzntl_subdomain) - .5) * params.dz, 1./3))
      , 0);
  }
  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::surf_sens()
{
  const auto &ijk = this->ijk;
  F(ijk).reindex(this->zero) = setup::F_sens * (*params.hgt_fctr_sclr)(this->vert_idx);

  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::surf_latent()
{
  const auto &ijk = this->ijk;
  F(ijk).reindex(this->zero) =  setup::F_lat * (*params.hgt_fctr_sclr)(this->vert_idx); // we need to use a reindexed view, because the profile's base is 0

  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::subsidence(const int &type) // large-scale vertical wind
{
  const auto &ijk = this->ijk;
  tmp1(ijk) = this->state(type)(ijk);
  this->vert_grad_cnt(tmp1, F, params.dz);
  F(ijk).reindex(this->zero) *= - (*params.w_LS)(this->vert_idx);

  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  this->smooth(tmp1, F);
}

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
  const auto &i = this->i;
  const auto &j = this->j;

  int nz = this->mem->grid_size[1].length(); //76

  // index of first cell above inversion
  tmp1(ijk)  = rv(ijk) + r_l(ijk);
  k_i(this->horizontal_subdomain) = blitz::first( tmp1 < setup::q_i, this->vert_idx); 
  F(ijk) = 0;

  // calc Eqs. 5 and 6 from Ackerman et al 2009
  auto sum = tmp1(this->horizontal_subdomain, 0); // alias for a 1D/2D temp array
  for(int z = 0 ; z < nz; ++z)
  {
    sum = blitz::sum(r_l(this->horizontal_subdomain, blitz::Range(z, nz-1)), this->vert_idx);
    if(z==0)
      F(this->horizontal_subdomain, z) += setup::F_0 * exp(- (nz - z - 1) * params.dz * sum); 
    else
      F(this->horizontal_subdomain, z) += setup::F_0 * exp(- (nz - z - 0.5) * params.dz * sum); 

    if(z > 0)
    {
      sum = blitz::sum(r_l(this->horizontal_subdomain, blitz::Range(0, z-1)), this->vert_idx);
      F(this->horizontal_subdomain, z) += setup::F_1 * exp(- (z - 0.5) * params.dz * sum);
    }

    F(this->horizontal_subdomain, z) += where(z > k_i(this->horizontal_subdomain), 
      setup::c_p * setup::rho_i * setup::D *
      (0.25 * pow((z - 0.5) * params.dz - (k_i(this->horizontal_subdomain) - .5) * params.dz, 4./3) +
      (k_i(this->horizontal_subdomain) - .5) * params.dz * pow((z - 0.5) * params.dz - (k_i(this->horizontal_subdomain) - .5) * params.dz, 1./3))
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

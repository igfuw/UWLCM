#pragma once

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv)
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;

  const real_t g = 9.81; 

  namespace moist_air = libcloudphxx::common::moist_air;
  const real_t eps = moist_air::R_v<real_t>() / moist_air::R_d<real_t>() - 1.;
  tmp1(ijk).reindex({0,0}) = g * ((th(ijk).reindex({0,0}) - (*params.th_e)(blitz::tensor::j)) / (*params.th_ref)(blitz::tensor::j)) + eps * (rv(ijk).reindex({0,0}) - (*params.rv_e)(blitz::tensor::j)) - r_l(ijk).reindex({0,0});

// smoothing
  this->xchng_sclr(tmp1, i, j); 
  F(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
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
  blitz::secondIndex ki;
  tmp1(ijk)  = rv(ijk) + r_l(ijk);
  k_i(i) = blitz::first( tmp1 < setup::q_i, ki); 
  F(ijk) = 0;

  // calc Eqs. 5 and 6 from Ackerman et al 2009
  for(int x = i.first() ; x <= i.last(); ++x) // 0..75 || 0..37 38..75 || ...
  {
    for(int z = 0 ; z < nz; ++z)
    {
      setup::real_t sum = blitz::sum(r_l(x, blitz::Range(z, nz-1)));
      if(z==0)
        F(x, z) += setup::F_0 * exp(- (nz - z - 1) * this->dj * sum); 
      else
        F(x, z) += setup::F_0 * exp(- (nz - z - 0.5) * this->dj * sum); 

      if(z > 0)
      {
        sum = blitz::sum(r_l(x, blitz::Range(0, z-1)));
        F(x, z) += setup::F_1 * exp(- (z - 0.5) * this->dj * sum);
      }

      if(z > k_i(x) )
      {
        real_t z_i = (k_i(x) - .5) * this->dj; // bottom of the first cell above inversion, z=0 at k=0.5
        real_t z_d = (z - 0.5) * this->dj - z_i;
        F(x, z) += setup::c_p * setup::rho_i * setup::D * (0.25 * pow(z_d, 4./3) + z_i * pow(z_d, 1./3)); 
      }
    }
  }
// smoothing
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j); 
  F(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::surf_sens()
{
  const auto &ijk = this->ijk;
  F(ijk).reindex({0,0}) = setup::F_sens * (*params.hgt_fctr_sclr)(blitz::tensor::j);
// smoothing
  const auto &i = this->i;
  const auto &j = this->j;
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j); 
  F(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::surf_latent()
{
  const auto &ijk = this->ijk;
  std::cout << ijk[0] << std::endl;
  std::cout << ijk[1] << std::endl;
  std::cout << F << std::endl;
  std::cout << (*params.hgt_fctr_sclr) << std::endl;
  F(ijk).reindex({0,0}) =  setup::F_lat * (*params.hgt_fctr_sclr)(blitz::tensor::j); // we need to use a reindexed view, because the profile's base is 0
// smoothing
  const auto &i = this->i;
  const auto &j = this->j;
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j); 
  F(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
}

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::subsidence(const int &type) // large-scale vertical wind
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  tmp1(ijk) = this->state(type)(ijk);
  this->xchng_sclr(tmp1, i, j);
  F(i, j).reindex({0,0}) = - (*params.w_LS)(blitz::tensor::j) * (tmp1(i, j + 1).reindex({0,0})- tmp1(i, j - 1).reindex({0,0})) / (2. * this->dj); 
// smoothing
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j); 
  F(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
}

#pragma once
#include "../formulae/nabla_formulae.hpp"
#include "../detail/subs_t.hpp"

//TODO: make these functions return arrays

template <class ct_params_t>
void slvr_common<ct_params_t>::subsidence(const int &type) // large-scale vertical wind
{
  const auto &ijk = this->ijk;
  if(params.subsidence == subs_t::local)
  {
    tmp1(ijk) = this->state(type)(ijk);
    this->vert_grad_cnt(tmp1, F, params.dz); 
    F(ijk).reindex(this->zero) *= - (*params.w_LS)(this->vert_idx);
    //this->smooth(tmp1, F);
  }
  else if(params.subsidence == subs_t::mean)
  {
    setup::arr_1D_t mean(this->vert_rng.length()+1),
                    grad(this->vert_rng.length());
    this->hrzntl_mean(this->state(type), mean);
    grad_fwd(mean, grad, params.dz, this->vert_rng); 
    F(ijk).reindex(this->zero) = - grad(this->vert_idx) * (*params.w_LS)(this->vert_idx);
    //this->smooth(tmp1, F);
  }
  else
    F(ijk)=0.;
}

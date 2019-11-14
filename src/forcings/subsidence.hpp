#pragma once

//TODO: make these functions return arrays

template <class ct_params_t>
void slvr_common<ct_params_t>::subsidence(const int &type) // large-scale vertical wind
{
  const auto &ijk = this->ijk;
  if(params.subsidence)
  {
    F(ijk) = this->state(type)(ijk);
    this->vert_grad_cnt(F, tmp1, params.dz); 
    tmp1(ijk).reindex(this->zero) *= - (*params.w_LS)(this->vert_idx);
    this->smooth(tmp1, F);
  }
  else
    F(ijk)=0.;
}

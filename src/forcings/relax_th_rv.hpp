#pragma once

template <class ct_params_t>
void slvr_common<ct_params_t>::relax_th_rv(const int &type)
{
  const auto &ijk = this->ijk;
  if(params.relax_th_rv)
  {
    if(type == ix::th)
      F(ijk).reindex(this->zero) = (*params.relax_th_rv_coeff)(this->vert_idx) * ((*params.th_e)(this->vert_idx) - th_mean_prof(this->vert_idx));// / this->n_cell_per_level;
    else if(type == ix::rv)
      F(ijk).reindex(this->zero) = (*params.relax_th_rv_coeff)(this->vert_idx) * ((*params.rv_e)(this->vert_idx) - rv_mean_prof(this->vert_idx));// / this->n_cell_per_level;
  }
  else
    F(ijk)=0.;
}

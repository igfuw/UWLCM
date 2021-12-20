#pragma once

template <class ct_params_t>
void slvr_common<ct_params_t>::nudging(const int &type)
{
  const auto &ijk = this->ijk;
  if(params.nudging)
  {
    if(type == ix::th)
      F(ijk).reindex(this->zero) = (*params.nudging_coeff)(this->vert_idx) * ((*params.th_e)(this->vert_idx) - th_mean_prof(this->vert_idx)) / this->n_cell_per_level;
    else if(type == ix::rv)
      F(ijk).reindex(this->zero) = (*params.nudging_coeff)(this->vert_idx) * ((*params.rv_e)(this->vert_idx) - rv_mean_prof(this->vert_idx)) / this->n_cell_per_level;
  }
  else
    F(ijk)=0.;
}

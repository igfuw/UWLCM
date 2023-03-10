#pragma once

// TODO: rv_LS and th_LS are very similar, make one function

template <class ct_params_t>
void slvr_common<ct_params_t>::rv_LS()
{
  params.update_rv_LS(
    params.profs.rv_LS, this->timestep, params.dt, params.dz
  );

  F(this->ijk).reindex(this->zero) = (params.profs.rv_LS)(this->vert_idx);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::th_LS()
{
  params.update_th_LS(
    params.profs.th_LS, this->timestep, params.dt, params.dz
  );

  F(this->ijk).reindex(this->zero) = (params.profs.th_LS)(this->vert_idx);
}

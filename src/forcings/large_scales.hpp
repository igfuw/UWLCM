#pragma once

// TODO: rv_LS and th_LS are very similar, make one function

template <class ct_params_t>
void slvr_common<ct_params_t>::calc_rv_LS()
{
  params.update_rv_LS(
    rv_LS(this->ijk).reindex(this->origin), this->timestep * params.dt, this->di, this->dj, this->dk
  );
  F(this->ijk) = rv_LS(this->ijk);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::calc_th_LS()
{
  params.update_th_LS(
    th_LS(this->ijk).reindex(this->origin), this->timestep * params.dt, this->di, this->dj, this->dk
  );

  F(this->ijk) = th_LS(this->ijk);
}

#pragma once

// TODO: rv_LS and th_LS are very similar, make one function

template <class ct_params_t>
void slvr_common<ct_params_t>::rv_LS()
{
  this->mem->barrier();
  if(this->rank == 0)
  {
    params.update_rv_LS(this->timestep, this->dt);
  }
  this->mem->barrier();

  F(this->ijk).reindex(this->zero) = (params.profs.rv_LS)(this->vert_idx);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::th_LS()
{
  this->mem->barrier();
  if(this->rank == 0)
  {
    params.update_th_LS(this->timestep, this->dt);
  }
  this->mem->barrier();

  F(this->ijk).reindex(this->zero) = (params.profs.th_LS)(this->vert_idx);
}

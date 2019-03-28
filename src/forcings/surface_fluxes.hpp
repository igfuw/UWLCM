#pragma once

//TODO: make these functions return arrays

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens_impl(iles_tag)
{
  params.update_surf_flux_sens(surf_flux_sens(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  //F(ijk).reindex(this->zero) = surf_flux_sens(this->hrzntl_subdomain)(blitz::tensor::i, blitz::tensor::j) 
  //                             * (*params.hgt_fctr_sclr)(this->vert_idx);
  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_sens(this->hrzntl_slice(0)) * (*params.hgt_fctr_sclr)(k);
  }

//  tmp1(ijk)=F(ijk); //TODO: unnecessary copy
  //this->smooth(tmp1, F);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens_impl(smg_tag)
{
  params.update_surf_flux_sens(surf_flux_sens(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  F(this->ijk) = 0;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens()
{
  surf_sens_impl(sgs_tag{});
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent_impl(iles_tag)
{
  params.update_surf_flux_lat(surf_flux_lat(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_lat(this->hrzntl_slice(0)) * (*params.hgt_fctr_sclr)(k);
  }
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent_impl(smg_tag)
{
  params.update_surf_flux_lat(surf_flux_lat(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  F(this->ijk) = 0;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent()
{
  surf_latent_impl(sgs_tag{});
}


template <class ct_params_t>
void slvr_common<ct_params_t>::surf_u_impl(iles_tag)
{
  params.update_surf_flux_u(surf_flux_u(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_u(this->hrzntl_slice(0)) * (*params.hgt_fctr_vctr)(k);
  }
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_u_impl(smg_tag)
{
  params.update_surf_flux_u(surf_flux_u(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  F(this->ijk) = 0;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_u()
{
  surf_u_impl(sgs_tag{});
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_v_impl(iles_tag)
{
  params.update_surf_flux_v(surf_flux_v(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_v(this->hrzntl_slice(0)) * (*params.hgt_fctr_vctr)(k);
  }
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_v_impl(smg_tag)
{
  params.update_surf_flux_u(surf_flux_v(this->hrzntl_slice(0)).reindex(this->origin), this->timestep, this->dt, this->di, this->dj);
  F(this->ijk) = 0;
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_v()
{
  surf_v_impl(sgs_tag{});
}

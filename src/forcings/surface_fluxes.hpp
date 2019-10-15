#pragma once

//TODO: make these functions return arrays

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens_impl(iles_tag)
{
  params.update_surf_flux_sens(
    surf_flux_sens(this->hrzntl_slice(0)).reindex(this->origin),
    this->state(ix::th)(this->hrzntl_slice(1)).reindex(this->origin),
    U_ground(this->hrzntl_slice(0)).reindex(this->origin), params.dz,
    this->timestep, this->dt, this->di, this->dj
  );

  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_sens(this->hrzntl_slice(0)) * (*params.hgt_fctr)(k);
  }
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_sens_impl(smg_tag)
{
  params.update_surf_flux_sens(
    surf_flux_sens(this->hrzntl_slice(0)).reindex(this->origin),
    this->state(ix::th)(this->hrzntl_slice(1)).reindex(this->origin),
    U_ground(this->hrzntl_slice(0)).reindex(this->origin), params.dz,
    this->timestep, this->dt, this->di, this->dj
  );

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
  params.update_surf_flux_lat(
    surf_flux_lat(this->hrzntl_slice(0)).reindex(this->origin),
    this->state(ix::rv)(this->hrzntl_slice(1)).reindex(this->origin) + this->r_l(this->hrzntl_slice(1)).reindex(this->origin), // when was r_l diagnosed? this or previous step? probably makes no difference (we are on the ground)
    U_ground(this->hrzntl_slice(0)).reindex(this->origin), params.dz,
    this->timestep, this->dt, this->di, this->dj
  );

  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_lat(this->hrzntl_slice(0)) * (*params.hgt_fctr)(k);
  }
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_latent_impl(smg_tag)
{
  params.update_surf_flux_lat(
    surf_flux_lat(this->hrzntl_slice(0)).reindex(this->origin),
    this->state(ix::rv)(this->hrzntl_slice(1)).reindex(this->origin) + this->r_l(this->hrzntl_slice(1)).reindex(this->origin), // when was r_l diagnosed? this or previous step? probably makes no difference (we are on the ground)
    U_ground(this->hrzntl_slice(0)).reindex(this->origin), params.dz,
    this->timestep, this->dt, this->di, this->dj
  );

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
  params.update_surf_flux_uv(
    surf_flux_u(this->hrzntl_slice(0)).reindex(this->origin),
    this->state(ix::vip_i)(this->hrzntl_slice(1)).reindex(this->origin),
    U_ground(this->hrzntl_slice(0)).reindex(this->origin), params.dz,
    this->timestep, this->dt, this->di, this->dj
  );

  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_u(this->hrzntl_slice(0)) * (*params.hgt_fctr)(k);
  }
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_u_impl(smg_tag)
{
  throw std::runtime_error("momentum surface (u) flux called in a Smagorinsky simulation.");
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_u()
{
  surf_u_impl(sgs_tag{});
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_v_impl(iles_tag)
{
  params.update_surf_flux_uv(
    surf_flux_v(this->hrzntl_slice(0)).reindex(this->origin),
    this->state(ix::vip_j)(this->hrzntl_slice(1)).reindex(this->origin),
    U_ground(this->hrzntl_slice(0)).reindex(this->origin), params.dz,
    this->timestep, this->dt, this->di, this->dj
  );

  for (auto k = this->vert_rng.first(); k <= this->vert_rng.last(); ++k)
  {
    F(this->hrzntl_slice(k)) = surf_flux_v(this->hrzntl_slice(0)) * (*params.hgt_fctr)(k);
  }
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_v_impl(smg_tag)
{
  throw std::runtime_error("momentum surface flux (v) called in a Smagorinsky simulation.");
}

template <class ct_params_t>
void slvr_common<ct_params_t>::surf_v()
{
  surf_v_impl(sgs_tag{});
}

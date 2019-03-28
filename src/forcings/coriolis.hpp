#pragma once

//TODO: make these functions return arrays
//

// Coriolis force including large-scale geostrophic wind
// F_u = + coriolis_parameter * (v - v_geostrophic)
// F_v = - coriolis_parameter * (u - u_geostrophic)
// the +- sign is handled by calc_forces
template <class ct_params_t>
void slvr_common<ct_params_t>::coriolis(
  const int &vel_idx
) 
{
  const auto &ijk = this->ijk;
  if(params.coriolis && ct_params_t::n_dims==3) // TODO: n_dims is known at compile time
  {
    F(ijk).reindex(this->zero) = params.ForceParameters.coriolis_parameter *
      (this->state(vel_idx)(ijk).reindex(this->zero) - (*params.geostr[vel_idx])(this->vert_idx));
  }
  else
    F(ijk)=0.;
}

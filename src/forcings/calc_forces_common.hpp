//TODO: move calc_forces and forcings to case class?
#pragma once
#include "buoyancy.hpp"
#include "coriolis.hpp"
#include "radiation.hpp"
#include "subsidence.hpp"
#include "surface_fluxes.hpp"
#include "../solvers/slvr_common.hpp"
#include "../detail/blitz_hlpr_fctrs.hpp"

// common forcing functions
// TODO: make functions return blitz arrays to avoid unnecessary copies
template <class ct_params_t>
void slvr_common<ct_params_t>::rv_src()
{
  const auto &ijk = this->ijk;
  if(params.rv_src)
  {
    // surface flux
    surf_latent();
    alpha(ijk) = F(ijk);

    // large-scale vertical wind
    subsidence(ix::rv);
    alpha(ijk) += F(ijk);

    // large-scale horizontal advection
    alpha(ijk).reindex(this->zero) += (*params.rv_LS)(this->vert_idx);
  }
  else
    alpha(ijk) = 0.;

  beta(ijk) = 0.;
  // nudging, todo: use some other coeff than vab_coeff
//  alpha(ijk).reindex(this->zero) += (*this->mem->vab_coeff)(ijk).reindex(this->zero) * (*params.rv_e)(this->vert_idx); // TODO: its a constant, cache it
//  beta(ijk) = - (*this->mem->vab_coeff)(ijk);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::th_src(typename parent_t::arr_t &rv)
{
  const auto &ijk = this->ijk;
  if(params.th_src)
  {
    // -- heating --

    // radiation
    radiation(rv);
    nancheck(radiative_flux(ijk), "radiation");
    
    // sum of th flux, F(j) is upward flux through the bottom of the j-th cell
    this->vert_grad_fwd(radiative_flux, alpha, params.dz);
    alpha(ijk) *= -1; // negative gradient means inflow
    nancheck(alpha(ijk), "sum of th flux");
    
    // change of theta[K/s] = heating[W/m^3] / exner / c_p[J/K/kg] / this->rhod[kg/m^3]
    alpha(ijk).reindex(this->zero) /=  calc_exner()((*params.p_e)(this->vert_idx)) *
      calc_c_p()(rv(ijk).reindex(this->zero)) *
      (*params.rhod)(this->vert_idx);
    nancheck2(alpha(ijk), this->state(ix::th)(ijk), "change of theta");

    // surf flux = d/dz mean(theta*w) [K/s]
    surf_sens();
    nancheck(F(ijk), "sensible surf forcing");
    alpha(ijk) += F(ijk);
  
    // large-scale vertical wind
    subsidence(ix::th);
    nancheck(F(ijk), "subsidence");
    alpha(ijk) += F(ijk);
    nancheck(alpha(ijk), "alpha in th_src");

    // large-scale horizontal advection
    alpha(ijk).reindex(this->zero) += (*params.th_LS)(this->vert_idx);
  }
  else
    alpha(ijk) = 0.;

  beta(ijk) = 0.;
  // nudging, todo: use some other coeff than vab_coeff
  //alpha(ijk).reindex(this->zero) += (*this->mem->vab_coeff)(ijk).reindex(this->zero) * (*params.th_e)(this->vert_idx);
  //beta(ijk) = - (*this->mem->vab_coeff)(ijk);
}

template <class ct_params_t>
void slvr_common<ct_params_t>::w_src(typename parent_t::arr_t &th, typename parent_t::arr_t &rv, const int at)
{
  const auto &ijk = this->ijk;
  // buoyancy
  // TODO: buoyancy is now calculated twice, at n and at n+1, make it so that it is calculated once (will need to remove zeroing-out of w rhs in parent:update-rhs)
  buoyancy(th, rv);
  alpha(ijk) = 0.5 * F(ijk); // halved, because it is applied trapezoidaly
  if(at == 0 && params.vel_subsidence) // subsidence added explicitly, so updated only at n
  {
    // large-scale vertical wind
    subsidence(ix::w);
    alpha(ijk) += F(ijk);
  }
}

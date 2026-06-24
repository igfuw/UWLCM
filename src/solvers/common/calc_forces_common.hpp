//TODO: move calc_forces and forcings to case class?
#pragma once
#include "../slvr_common.hpp"
#include "../../detail/blitz_hlpr_fctrs.hpp"
#include "../../forcings/buoyancy.hpp"
#include "../../forcings/coriolis.hpp"
#include "../../forcings/radiation.hpp"
#include "../../forcings/subsidence.hpp"
#include "../../forcings/relax_th_rv.hpp"
#include "../../forcings/surface_fluxes.hpp"
#include "../../forcings/large_scales.hpp"

/**
 * @brief Apply source term for water vapor (rv) due to surface fluxes, large-scale vertical motion,
 * horizontal advection, and per-level nudging.
 *
 * If `params.rv_src` is true, the function sequentially applies:
 *  - surface latent heat flux,
 *  - subsidence (large-scale vertical wind),
 *  - large-scale horizontal advection,
 *  - nudging of the mean water vapor.
 * The forcing coefficient `alpha` accumulates contributions from all sources, and `beta` is set to zero.
 */
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
    rv_LS();
    alpha(ijk) += F(ijk);

    // per-level nudging of the mean
    relax_th_rv(ix::rv);
    alpha(ijk) += F(ijk);
  }
  else
    alpha(ijk) = 0.;

  beta(ijk) = 0.;
}


/**
 * @brief Apply source term for potential temperature (th) due to radiation, surface fluxes,
 * large-scale motions, and per-level nudging.
 *
 * If `params.th_src` is true, the function sequentially applies:
 *  - radiative heating,
 *  - vertical flux divergence to compute local heating rate,
 *  - surface sensible heat flux,
 *  - subsidence (large-scale vertical wind),
 *  - large-scale horizontal advection,
 *  - nudging of the mean potential temperature.
 *
 * The forcing coefficient `alpha` accumulates contributions from all sources, corrected for
 * specific heat and density. The `beta` coefficient is set to zero.
 *
 * @param rv Array of water vapor used for computing heat capacities and radiative effects.
 */
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
    th_LS();
    alpha(ijk) += F(ijk);

    // per-level nudging of the mean
    relax_th_rv(ix::th);
    alpha(ijk) += F(ijk);
  }
  else
    alpha(ijk) = 0.;

  beta(ijk) = 0.;
}


/**
 * @brief Apply source term for vertical velocity (w) due to buoyancy and optionally subsidence.
 *
 * The function computes buoyancy forcing based on `th` and `rv`, applies trapezoidal scaling
 * (halving `alpha`), and optionally adds large-scale vertical motion if `at == 0` and
 * `params.vel_subsidence` is true.
 *
 * @param th Array of potential temperature.
 * @param rv Array of water vapor.
 * @param at Integration stage (0 for first stage, 1 for trapezoidal second stage).
 */
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

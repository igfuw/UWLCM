#pragma once

namespace formulas
{
  // scaling of bulk surface flux coefficients following Stevens et al. JAS 2001
  template<class real_t>
  real_t surf_flux_coeff_scaling(const real_t &z, const real_t &z_ref) // [m]
  {
    static const real_t z_0 = 1.5e-4; // [m]
    return pow(log(z_ref / z_0) / log(z / z_0), real_t(2));
  }

  // function used for sensible and latent heat fluxes from Grabowski et al. (2006)
  template<class real_t>
  real_t surf_flux_function(real_t t) //[s]
  {
    real_t t_hours = t / real_t(3600);
    real_t hrl = real_t(7.5) + t_hours;
    real_t thea = boost::math::constants::pi<real_t>() / real_t(2) * (real_t(12.75) - hrl) / real_t(5.25);
    real_t xfact = cos(thea);
    if (xfact < real_t(0)) xfact = real_t(0);
    return xfact;
  }

};

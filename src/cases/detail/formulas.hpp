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
  real_t surf_flux_function(const real_t &t) //[s]
  {
    static const real_t t_hours = t / real_t(3600);
    return std::max(real_t(0), cos(boost::math::constants::pi<real_t>() * (real_t(5.25) - t_hours) / real_t(10.5)));
  }

};

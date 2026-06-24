// Cases in which envionmental profile is set like in babyEulag

#pragma once
#include "CasesCommon.hpp"

namespace cases 
{
  /**
 * \class Anelastic
 * @brief Base class for anelastic cloud simulation cases.
 *
 * Provides initialization of environmental and reference profiles of potential
 * temperature and water vapor, following the anelastic approximation. Designed
 * to be extended by specific cases with defined theta and mixing ratio profiles.
 *
 * @tparam case_ct_params_t Compile-time parameters for the case
 * @tparam n_dims Number of spatial dimensions (2 or 3)
 */
  template<class case_ct_params_t, int n_dims>
  class Anelastic : public CasesCommon<case_ct_params_t, n_dims>
  {
    protected:
    using parent_t = CasesCommon<case_ct_params_t, n_dims>;

    quantity<si::pressure, real_t> p_0 = -10 * si::pascals; 
    quantity<si::mass_density, real_t> rhod_0 = 0 * si::kilograms / si::cubic_meters; 

    /**
 * @brief Liquid water potential temperature at height z.
 *
 * Must be implemented in derived classes.
 * @param z Height (m)
 * @return Potential temperature (K)
 */
    virtual quantity<si::temperature, real_t> th_l(const real_t &z) {throw std::runtime_error("UWLCM: base Anelastic class th_l called");}

    /**
  * @brief Total water mixing ratio at height z.
  *
  * Must be implemented in derived classes.
  * @param z Height (m)
  * @return Water mixing ratio (dimensionless)
  */
    virtual quantity<si::dimensionless, real_t> r_t(const real_t &z) {throw std::runtime_error("UWLCM: base Anelastic class r_t called");}

    /**
 * @brief Initialize environmental profiles of theta and water vapor.
 *
 * Calculates initial temperature, pressure, water vapor, and liquid water
 * mixing ratios using a moist thermodynamic framework.
 *
 * @param profs Reference to profiles structure to populate
 * @param nz Number of vertical levels
 */
    void env_prof(detail::profiles_t &profs, int nz)
    {
      using libcloudphxx::common::moist_air::R_d_over_c_pd;
      using libcloudphxx::common::moist_air::c_pd;
      using libcloudphxx::common::moist_air::R_d;
      using libcloudphxx::common::const_cp::l_tri;
      using libcloudphxx::common::theta_std::p_1000;

      assert(p_0 > 0. * si::pascals);

      // temp profile
      arr_1D_t T(nz);
      real_t dz = (this->Z / si::metres) / (nz-1);

      // input sounding at z=0, for moist air, no liquid water
      T(0) = th_l(0.) / si::kelvins *  pow(p_0 / p_1000<real_t>(),  R_d_over_c_pd<real_t>());
      profs.p_e(0) = p_0 / si::pascals;
      profs.th_e(0) = th_l(0.) / si::kelvins;
      profs.rv_e(0) = r_t(0.);
      profs.rl_e(0) = 0.;

      real_t tt0 = 273.17;
      real_t rv = 461; // specific gas constant for vapor
      real_t ee0 = 611.;
      real_t a = R_d<real_t>() / rv / si::joules * si::kelvins * si::kilograms; // aka epsilon
      real_t b = l_tri<real_t>() / si::joules * si::kilograms / rv / tt0;
      real_t c = l_tri<real_t>() / c_pd<real_t>() / si::kelvins;
      real_t d = l_tri<real_t>() / si::joules * si::kilograms / rv;
      real_t f = R_d_over_c_pd<real_t>(); 

      real_t lwp_env = 0;

      for(int k=1; k<nz; ++k)
      {
        real_t bottom = R_d<real_t>() / si::joules * si::kelvins * si::kilograms * T(k-1) * (1 + 0.61 * profs.rv_e(k-1)); // (p / rho) of moist air at k-1
        real_t rho1 = profs.p_e(k-1) / bottom; // rho at k-1
        profs.p_e(k) = profs.p_e(k-1) - rho1 * 9.81 * dz; // estimate of pre at k (dp = -g * rho * dz)
        real_t thetme = pow(p_1000<real_t>() / si::pascals / profs.p_e(k), f); // 1/Exner
        real_t thi = 1. / (th_l(k * dz) / si::kelvins); // 1/theta_std
        real_t y = b * thetme * tt0 * thi; 
        real_t ees = ee0 * exp(b-y); // saturation vapor pressure (Tetens equation or what?)
        real_t qvs = a * ees / (profs.p_e(k) - ees);  // saturation vapor mixing ratio = R_d / R_v * ees / p_d
        // calculate linearized condensation rate
        real_t cf1 = thetme*thetme*thi*thi;  // T^{-2}
        cf1 *= c * d * profs.p_e(k) / (profs.p_e(k) - ees); // = l_tri^2 / (C_pd * R_v * T^2) * p/p_d
        real_t delta = (r_t(k*dz) - qvs) / (1 + qvs * cf1); // how much supersaturated is the air (divided by sth)
        if(delta < 0.) delta = 0.;
        profs.rv_e(k) = r_t(k*dz) - delta;
        profs.rl_e(k) = delta;
        profs.th_e(k) = th_l(k*dz) / si::kelvins + c * thetme * delta;
        T(k) = profs.th_e(k) * pow(profs.p_e(k) / (p_1000<real_t>() / si::pascals),  f);

        bottom = R_d<real_t>() / si::joules * si::kelvins * si::kilograms * T(k) * (1 + 0.61 * profs.rv_e(k)); // (p / rho) of moist air at k-1
        rho1 = profs.p_e(k) / bottom; // rho at k-1
        lwp_env  += delta * rho1;
      }
      lwp_env = lwp_env * 5  * 1e3;
    }


    /**
     * @brief Initialize reference profiles for theta and dry air density.
     *
     * Computes reference potential temperature and dry air density profiles based
     * on environmental profiles. Used for stability and anelastic calculations.
     *
     * @param profs Reference to profiles structure to populate
     * @param nz Number of vertical levels
     */
    void ref_prof(detail::profiles_t &profs, int nz)
    {
      using libcloudphxx::common::moist_air::R_d_over_c_pd;
      using libcloudphxx::common::moist_air::c_pd;
      using libcloudphxx::common::moist_air::R_d;
      using libcloudphxx::common::theta_std::p_1000;

      real_t dz = (this->Z / si::metres) / (nz-1);

      // compute reference state theta and rhod
      blitz::firstIndex k;
      // calculate average stability
      blitz::Range notopbot(1, nz-2);
      arr_1D_t st(nz);
      st=0;
      st(notopbot) = (profs.th_e(notopbot+1) - profs.th_e(notopbot-1)) / profs.th_e(notopbot);
      real_t st_avg = blitz::sum(st) / (nz-2) / (2.*dz);
      // reference theta
      profs.th_ref = profs.th_e(0) * (1. + 0.608 * profs.rv_e(0)) * exp(st_avg * k * dz);
    //  th_ref = th_e(0) * pow(1 + rv_e(0) / a, f) // calc dry theta at z=0 
    //           * exp(st_avg * k * dz);
      // virtual temp at surface

      real_t T_surf = profs.th_e(0) *  pow(p_0 / p_1000<real_t>(),  R_d_over_c_pd<real_t>());

      real_t T_virt_surf = T_surf * (1. + 0.608 * profs.rv_e(0));
      real_t rho_surf = (p_0 / si::pascals) / T_virt_surf / 287. ; // TODO: R_d instead of 287, its the total, not dry density!
      rhod_0 = rho_surf * si::kilograms / si::cubic_meters;
//      rho_surf /= (1 + rv_e(0)); // turn it into dry air density! TODO: is this correct? TODO2: approp change in the paper

   //   real_t rho_surf = (p_0 / si::pascals) / T_surf / (1. + 29. / 18. * rv_e(0)) / 287. ; // dry air density at the surface TODO: R_d instead of 287

      // real_t cs = 9.81 / (c_pd<real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / T_surf; // this is from Wojtek
       real_t cs = 9.81 / (c_pd<real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / (profs.th_e(0)
           * (1. + 0.608 * profs.rv_e(0))); 
      // rhod profile
      profs.rhod = rho_surf * exp(- st_avg * k * dz) * pow(
               1. - cs * (1 - exp(- st_avg * k * dz)), (1. / R_d_over_c_pd<real_t>()) - 1);


      // theta_std env prof to theta_dry_e
//      for(int k=1; k<nz; ++k)
//        th_e(k) = theta_dry::std2dry<real_t>(th_e(k) * si::kelvins, quantity<si::dimensionless, real_t>(rv_e(k))) / si::kelvins;
    }
  };
};

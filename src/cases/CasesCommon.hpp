#pragma once

#include <iostream>

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>


// simulation parameters container
// TODO: write them to rt_params directly in main()
struct user_params_t
{
  int nt, outfreq, spinup, rng_seed;
  setup::real_t dt, z_rlx_sclr;
  std::string outdir;
  bool serial, th_src, rv_src, uv_src, w_src;
};

namespace setup 
{
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  template<class concurr_t>
  class CasesCommon
  {
    public:
    // some not really common vars, used only in Dycoms but need to be defined everywhere
    const real_t z_i  = 795; //initial inversion height
    const real_t heating_kappa = 85; // m^2/kg
    const real_t F_0 = 70; // w/m^2
    const real_t F_1 = 22; // w/m^2
    const real_t q_i = 8e-3; // kg/kg
    const real_t c_p = 1004; // J / kg / K
  
    const real_t D = 3.75e-6; // large-scale wind horizontal divergence [1/s]
    const real_t rho_i = 1.12; // kg/m^3
  
    const real_t F_sens = 16; //W/m^2, sensible heat flux
    const real_t F_lat = 93; //W/m^2, latent heat flux
    const real_t u_fric = 0.25; // m/s, friction velocity
  
    //aerosol bimodal lognormal dist. 
    const quantity<si::length, real_t>
      mean_rd1 = real_t(.011e-6) * si::metres,
      mean_rd2 = real_t(.06e-6) * si::metres;
    const quantity<si::dimensionless, real_t>
      sdev_rd1 = real_t(1.2),
      sdev_rd2 = real_t(1.7);
    const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
      n1_stp = real_t(125e6*2) / si::cubic_metres, // 125 || 31
      n2_stp = real_t(65e6*2) / si::cubic_metres;  // 65 || 16
    const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
      n_unit_test = real_t(1) / si::cubic_metres;
  
    //aerosol lognormal dist. for GCCN from Jorgen Jensen
    const quantity<si::length, real_t>
      mean_rd3 = real_t(.283e-6) * si::metres;
    const quantity<si::dimensionless, real_t>
      sdev_rd3 = real_t(2.235);
    const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
      n3_stp = real_t(2.216e6) / si::cubic_metres;
  
    //aerosol chemical composition parameters (needed for activation)
    // for lgrngn:
    const quantity<si::dimensionless, real_t> kappa = .61; // ammonium sulphate; CCN-derived value from Table 1 in Petters and Kreidenweis 2007
    const quantity<si::dimensionless, real_t> kappa_gccn = 1.28; // NaCl; CCN-derived value from Table 1 in Petters and Kreidenweis 2007
    // for blk_2m:
    const quantity<si::dimensionless, real_t> chem_b = .55; //ammonium sulphate //chem_b = 1.33; // sodium chloride
  
    //th, rv and surface fluxes relaxation time and height
    const quantity<si::time, real_t> tau_rlx = 300 * si::seconds;
    const quantity<si::length, real_t> z_rlx_vctr = 1 * si::metres;

    // lognormal aerosol distribution
    template <typename T>
    struct log_dry_radii : public libcloudphxx::common::unary_function<T>
    {
      T funval(const T lnrd) const
      {
        return T((
            lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, real_t>(lnrd)) +
            lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
          ) * si::cubic_metres
        );
      }
  
      log_dry_radii *do_clone() const 
      { return new log_dry_radii( *this ); }
    };
  
    // unit test lognormal aerosol distribution
    template <typename T>
    struct log_dry_radii_unit_test : public libcloudphxx::common::unary_function<T>
    {
      T funval(const T lnrd) const
      {
        return T((
            lognormal::n_e(mean_rd1, sdev_rd1, n_unit_test, quantity<si::dimensionless, real_t>(lnrd)) 
          ) * si::cubic_metres
        );
      }
  
      log_dry_radii_unit_test *do_clone() const 
      { return new log_dry_radii_unit_test( *this ); }
    };
  
    // lognormal aerosol distribution with GCCN
    template <typename T>
    struct log_dry_radii_gccn : public libcloudphxx::common::unary_function<T>
    {
      T funval(const T lnrd) const
      {
        return T((
            lognormal::n_e(mean_rd3, sdev_rd3, n3_stp, quantity<si::dimensionless, real_t>(lnrd)) 
          ) * si::cubic_metres
        );
      }
  
      log_dry_radii_gccn *do_clone() const 
      { return new log_dry_radii_gccn( *this ); }
    };

    protected:
  
    // function enforcing cyclic values in horizontal directions
    // 2D version
    template<class arr_t>
    void make_cyclic(arr_t arr,
      typename std::enable_if<arr_t::rank_ == 2>::type* = 0)
    { arr(arr.extent(0) - 1, blitz::Range::all()) = arr(0, blitz::Range::all()); }
  
    // 3D version
    template<class arr_t>
    void make_cyclic(arr_t arr,
      typename std::enable_if<arr_t::rank_ == 3>::type* = 0)
    { 
      arr(arr.extent(0) - 1, blitz::Range::all(), blitz::Range::all()) = 
        arr(0, blitz::Range::all(), blitz::Range::all()); 
      arr(blitz::Range::all(), arr.extent(1) - 1, blitz::Range::all()) = 
        arr(blitz::Range::all(), 0, blitz::Range::all());
    }
  };
};

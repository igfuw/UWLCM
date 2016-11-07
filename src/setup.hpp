#pragma once
#include <blitz/array.h> 
#include <libcloudph++/common/unary_function.hpp>
#include <libcloudph++/common/lognormal.hpp>

namespace setup 
{
  using real_t = float;
  using arr_1D_t = blitz::Array<setup::real_t, 1>;



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
  
  namespace lognormal = libcloudphxx::common::lognormal;

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
};

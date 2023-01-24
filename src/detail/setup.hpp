#pragma once
#include <blitz/array.h> 
#include <libcloudph++/common/unary_function.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <UWLCM/setup.hpp>
#include <chrono>

namespace setup 
{
  using arr_1D_t = blitz::Array<setup::real_t, 1>;

  using clock = std::chrono::system_clock;
  using timer = std::chrono::milliseconds;

/*
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n_unit_test = real_t(1) / si::cubic_metres;
  */
  //aerosol lognormal dist. for GCCN from Jorgen Jensen
  const quantity<si::length, real_t>
    mean_rd3 = real_t(.283e-6) * si::metres;
  const quantity<si::dimensionless, real_t>
    sdev_rd3 = real_t(2.235);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n3_stp = real_t(2.216e6) / si::cubic_metres;
  
  //aerosol chemical composition parameters (needed for activation)
  // for lgrngn:
  const quantity<si::dimensionless, real_t> kappa_gccn = 1.28; // NaCl; CCN-derived value from Table 1 in Petters and Kreidenweis 2007
  // for blk_2m:
  const quantity<si::dimensionless, real_t> chem_b = .55; //ammonium sulphate //chem_b = 1.33; // sodium chloride
  
  namespace lognormal = libcloudphxx::common::lognormal;

  // lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii : public libcloudphxx::common::unary_function<T>
  {
    const quantity<si::length, real_t> mean_rd1, mean_rd2;
    const quantity<si::dimensionless, real_t> sdev_rd1, sdev_rd2;
    const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_stp, n2_stp;

    log_dry_radii(
      quantity<si::length, real_t> mean_rd1,
      quantity<si::length, real_t> mean_rd2,
      quantity<si::dimensionless, real_t> sdev_rd1,
      quantity<si::dimensionless, real_t> sdev_rd2,
      quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_stp,
      quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n2_stp
    ):
    mean_rd1(mean_rd1),
    mean_rd2(mean_rd2),
    sdev_rd1(sdev_rd1),
    sdev_rd2(sdev_rd2),
    n1_stp(n1_stp),
    n2_stp(n2_stp) {}


    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }
  };
  
  // unit test lognormal aerosol distribution
   /*
  template <typename T>
  struct log_dry_radii_unit_test : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(real_t(1e-7) * si::metres, quantity<si::dimensionless, real_t>(1), n_unit_test, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii_unit_test *do_clone() const 
    { return new log_dry_radii_unit_test( *this ); }
  };
*/
  // lognormal aerosol distribution with GCCN
  template <typename T>
  struct log_dry_radii_gccn : public libcloudphxx::common::unary_function<T>
  {
    real_t lnrd_min, lnrd_max, conc_multiplier;

    T funval(const T lnrd) const
    {
      return 
        lnrd < lnrd_min ? 0 :
          lnrd > lnrd_max ? 0 :
            T((
              lognormal::n_e(mean_rd3, sdev_rd3, conc_multiplier * n3_stp, quantity<si::dimensionless, real_t>(lnrd)) 
            ) * si::cubic_metres
            );
    }

    log_dry_radii_gccn(const real_t lnrd_min = 0, const real_t lnrd_max = 1000000, const real_t conc_multiplier = 1): 
      lnrd_min(lnrd_min),
      lnrd_max(lnrd_max),
      conc_multiplier(conc_multiplier)
    {}
  };
};

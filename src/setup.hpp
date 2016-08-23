#pragma once

#include <iostream>

#include <blitz/array.h> 

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

// TODO: relaxation terms still missing

// 8th ICMW case 1 by Wojciech Grabowski)
namespace setup 
{
  using real_t = float;

  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  enum {x, z}; // dimensions
  const quantity<si::pressure, real_t> 
    p_0 = 101780 * si::pascals;
  const quantity<si::length, real_t> 
    z_0  = 0    * si::metres,
    Z    = 1500 * si::metres, // DYCOMS: 1500
    X    = 6400 * si::metres, // DYCOMS: 6400
    Y    = 6400 * si::metres; // DYCOMS: 6400
  const real_t z_i  = 795; //initial inversion height
  const real_t heating_kappa = 85; // m^2/kg
  const real_t F_0 = 70; // w/m^2
  const real_t F_1 = 22; // w/m^2
  const real_t q_i = 8e-3; // kg/kg
  const real_t c_p = 1004; // J / kg / K
  const real_t z_abs = 1250; // [m] height above which absorber works

  const real_t D = 3.75e-6; // large-scale wind horizontal divergence [1/s]
  const real_t rho_i = 1.12; // kg/m^3

  const real_t F_sens = 16; //W/m^2, sensible heat flux
  const real_t F_lat = 93; //W/m^2, latent heat flux
  const real_t u_fric = 0.25; // m/s, friction velocity

  // liquid water potential temperature at height z
  quantity<si::temperature, real_t> th_l(const real_t &z)
  {
    quantity<si::temperature, real_t> ret;
    ret = z < z_i ?
      288.3 * si::kelvins : 
      (295. + pow(z - z_i, 1./3)) * si::kelvins;
    return ret;
  }

  // water mixing ratio at height z
  struct r_t
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
      const quantity<si::dimensionless, real_t> q_t = z < z_i ?
        9.45e-3 : 
        (5. - 3. * (1. - exp((z_i - z)/500.))) * 1e-3;
      return q_t;
    }
    BZ_DECLARE_FUNCTOR(r_t);
  };

  // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
  struct th_dry_fctr
  {
    real_t operator()(const real_t &z) const
    {
      return theta_dry::std2dry<real_t>(th_l(z), r_t()(z)) / si::kelvins;
    }
    BZ_DECLARE_FUNCTOR(th_dry_fctr);
  };

  // westerly wind
  struct u
  {
    real_t operator()(const real_t &z) const
    {
      return 3. + 4.3 * z / 1000.; 
    }
    BZ_DECLARE_FUNCTOR(u);
  };

  // southerly wind
  struct v
  {
    real_t operator()(const real_t &z) const
    {
      return -9. + 5.6 * z / 1000.; 
    }
    BZ_DECLARE_FUNCTOR(v);
  };

  // large-scale vertical wind
  struct w_LS_fctr
  {
    real_t operator()(const real_t &z) const
    {
      return -D * z; 
    }
    BZ_DECLARE_FUNCTOR(w_LS_fctr);
  };

  // density profile as a function of altitude
  // hydrostatic and assuming constant theta (not used now)
  struct rhod_fctr
  {
    real_t operator()(real_t z) const
    {
      quantity<si::pressure, real_t> p = hydrostatic::p(
	z * si::metres, th_dry_fctr()(0.) * si::kelvins, r_t()(0.), z_0, p_0
      );
      
      quantity<si::mass_density, real_t> rhod = theta_std::rhod(
	p, th_dry_fctr()(0.) * si::kelvins, r_t()(0.)
      );

      return rhod / si::kilograms * si::cubic_metres;
    }

    // to make the rhod() functor accept Blitz arrays as arguments
    BZ_DECLARE_FUNCTOR(rhod_fctr);
  };


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
    n1_stp_pristine = real_t(125e6 * 2 * 0.5) / si::cubic_metres, // 125 || 31
    n2_stp_pristine = real_t(65e6 * 2 * 0.5) / si::cubic_metres;  // 65 || 16

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


  // function expecting a libmpdata solver parameters struct as argument
  template <class T>
  void setopts(T &params, int nx, int nz)
  {
    params.dx = (X / si::metres) / (nx-1); 
    params.dz = (Z / si::metres) / (nz-1);
    params.di = params.dx;
    params.dj = params.dz;
  }

  // function expecting a libmpdata++ solver as argument
  template <class concurr_t>
  void intcond(concurr_t &solver, blitz::Array<setup::real_t, 2> &rhod, int rng_seed)
  {
    using ix = typename concurr_t::solver_t::ix;

    // helper ondex placeholders
    blitz::firstIndex i;
    blitz::secondIndex k;

    // dx, dy ensuring 1500x1500 domain
    int 
      nx = solver.advectee().extent(x), 
      nz = solver.advectee().extent(z); 
    real_t 
      dx = (X / si::metres) / (nx-1), 
      dz = (Z / si::metres) / (nz-1); 

    // initial potential temperature & water vapour mixing ratio profiles
    solver.advectee(ix::th) = th_dry_fctr()(k * dz); 
    // randomly prtrb tht
    std::random_device rd;
    auto seed = rd();
    if(rng_seed > 0)
      seed = rng_seed;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(-0.1, 0.1);

    blitz::Array<real_t, 2> prtrb(nx, nz);
    for (int ii = 0; ii < nx; ++ii)
    {
      for (int kk = 0; kk < nz; ++kk)
      {
         prtrb(ii, kk) = dis(gen);
      }
    }
    auto i_r = blitz::Range(0, nx - 1);
    auto k_r = blitz::Range(0, nz - 1);

    // enforce cyclic perturbation
    prtrb(nx - 1, k_r) = prtrb(0, k_r);

    solver.advectee(ix::th)(i_r, k_r) += prtrb(i_r, k_r);

    solver.advectee(ix::rv) = r_t()(k * dz); 

    solver.advectee(ix::u) = 0;
    solver.advectee(ix::u)(i_r, k_r)= setup::u()(k * dz);
    solver.advectee(ix::w) = 0;  
   
    // absorbers
    solver.vab_coefficient() = where(k * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (k * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);

    solver.vab_relaxed_state(0) = solver.advectee(ix::u);
    solver.vab_relaxed_state(1) = 0;

    // density profile
    solver.g_factor() = rhod; 
  }

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

  // pristine lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii_pristine : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd1, sdev_rd1, n1_stp_pristine, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd2, sdev_rd2, n2_stp_pristine, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii_pristine *do_clone() const 
    { return new log_dry_radii_pristine( *this ); }
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


  // calculate the initial environmental theta and rv profiles
  // like in Wojtek's BabyEulag
  void env_prof(blitz::Array<setup::real_t, 2> &th_e, blitz::Array<setup::real_t, 2> &rv_e, blitz::Array<setup::real_t, 2> &th_ref, blitz::Array<setup::real_t, 2> &rhod, int nz)
  {
    using libcloudphxx::common::moist_air::R_d_over_c_pd;
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::moist_air::R_d;
    using libcloudphxx::common::const_cp::l_tri;
    using libcloudphxx::common::theta_std::p_1000;

    blitz::Range all(blitz::Range::all());
    // pressure profile
    blitz::Array<setup::real_t, 1> pre(nz);
    // temperature profile
    blitz::Array<setup::real_t, 1> T(nz);
    setup::real_t dz = (Z / si::metres) / (nz-1);


    r_t rt;
    T(0) = th_l(0.) / si::kelvins *  pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    pre(0) = setup::p_0 / si::pascals;
    th_e(all, 0) = th_l(0.) / si::kelvins;
    rv_e(all, 0) = rt(0.);

    setup::real_t tt0 = 273.17;
    setup::real_t rv = 461;
    setup::real_t ee0 = 611.;
    setup::real_t a = R_d<setup::real_t>() / rv / si::joules * si::kelvins * si::kilograms;
    setup::real_t b = l_tri<setup::real_t>() / si::joules * si::kilograms / rv / tt0;
    setup::real_t c = l_tri<setup::real_t>() / c_pd<setup::real_t>() / si::kelvins;
    setup::real_t d = l_tri<setup::real_t>() / si::joules * si::kilograms / rv;
    setup::real_t f = R_d_over_c_pd<setup::real_t>(); 

    for(int k=1; k<nz; ++k)
    {
      setup::real_t bottom = R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms * T(k-1) * (1 + 0.61 * rv_e(0, k-1));
      setup::real_t rho1 = pre(k-1) / bottom;
      pre(k) = pre(k-1) - rho1 * 9.81 * dz;
      setup::real_t thetme = pow(p_1000<setup::real_t>() / si::pascals / pre(k), f);
      setup::real_t thi = 1. / (th_l(k * dz) / si::kelvins);
      setup::real_t y = b * thetme * tt0 * thi; 
      setup::real_t ees = ee0 * exp(b-y);
      setup::real_t qvs = a * ees / (pre(k) - ees); 
      setup::real_t cf1 = thetme*thetme*thi*thi;
      cf1 *= c * d * pre(k) / (pre(k) - ees);
      setup::real_t delta = (rt(k*dz) - qvs) / (1 + qvs * cf1);
      if(delta < 0.) delta = 0.;
      rv_e(all, k) = rt(k*dz) - delta;
      th_e(all, k) = th_l(k*dz) / si::kelvins + c * thetme * delta;
      T(k) = th_e(0, k) * pow(pre(k) / (p_1000<setup::real_t>() / si::pascals),  f);
    }

    // compute reference state theta and rhod
    blitz::secondIndex k;
    // calculate average stability
    blitz::Range notopbot(1, nz-2);
    blitz::Array<setup::real_t, 1> st(nz);
    st=0;
    st(notopbot) = (th_e(0, notopbot+1) - th_e(0, notopbot-1)) / th_e(0, notopbot);
    setup::real_t st_avg = blitz::sum(st) / (nz-2) / (2.*dz);
    // reference theta
    th_ref = th_e(0,0) * exp(st_avg * k * dz);
    // virtual temp at surface
    using libcloudphxx::common::moist_air::R_d_over_c_pd;
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::moist_air::R_d;
    using libcloudphxx::common::theta_std::p_1000;

    setup::real_t T_surf = th_e(0, 0) *  pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    setup::real_t T_virt_surf = T_surf * (1. + 0.608 * rv_e(0, 0));
    setup::real_t rho_surf = (setup::p_0 / si::pascals) / T_virt_surf / 287. ; // TODO: R_d instead of 287
    setup::real_t cs = 9.81 / (c_pd<setup::real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / T_surf;
    // rhod profile
    rhod = rho_surf * exp(- st_avg * k * dz) * pow(
             1. - cs * (1 - exp(- st_avg * k * dz)), (1. / R_d_over_c_pd<setup::real_t>()) - 1);
  }
};

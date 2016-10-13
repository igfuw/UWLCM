#pragma once

#include <iostream>

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

namespace setup 
{
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;
  namespace moist_air = libcloudphxx::common::moist_air;
  namespace const_cp = libcloudphxx::common::const_cp;

  // RH to rv
  quantity<si::dimensionless, real_t> RH_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
  {
    return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
  }
  // theta std to temperature
  template <class real_t>
  quantity<si::temperature, real_t> th2T(const quantity<si::temperature, real_t> &th, const quantity<si::pressure, real_t> &p) constexpr
  {
    quantity<si::temperature, real_t> T = th * pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    return T;
  }

  const real_t env_RH = 0.2;
  const real_t prtrb_RH = 1.00;

  const quantity<si::temperature, real_t>
    T_0(283. * si::kelvins);  // surface temperature
  const quantity<si::pressure, real_t> 
    p_0 = 85000 * si::pascals;
  const quantity<si::dimensionless, real_t> rv_0(RH_to_rv(env_RH, T_0, p_0));
  // theta (std) at surface
  const quantity<si::temperature, real_t> th_0 = T_0 / pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
  const quantity<si::temperature, real_t> th_0_dry = theta_dry::std2dry<real_t>(th_0, rv_0);
  const real_t S = 1.3e-5; // stability, 1/m

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

  // standard potential temperature at height z
  quantity<si::temperature, real_t> th_std(const real_t &z)
  {
    quantity<si::temperature, real_t> ret;
    ret = th_0 * exp(S * z);

    return ret;
  }


  struct env_rv
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
      return RH_to_rv(env_RH, th2T(th(z), p(z)), p(z));
    }
  BZ_DECLARE_FUNCTOR(env_rv);
  };

  struct prtrb_rv
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
      return RH_to_rv(prtrb_RH, th2T(th(z), p(z)), p(z));
    }
  BZ_DECLARE_FUNCTOR(prtrb_rv);
  };


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
      return theta_dry::std2dry<real_t>(th_std(z), r_t()(z)) / si::kelvins;
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


  template <class T, class U>
  void setopts_hlpr(T &params, const U &user_params)
  {
    params.outdir = user_params.outdir;
    params.outfreq = user_params.outfreq;
    params.spinup = user_params.spinup;
    params.relax_th_rv = user_params.relax_th_rv;
    params.w_src = user_params.w_src;
    params.uv_src = user_params.uv_src;
    params.th_src = user_params.th_src;
    params.rv_src = user_params.rv_src;
    params.prs_tol=1e-6;
    params.dt = user_params.dt;
    params.nt = user_params.nt;
  }

  // function expecting a libmpdata solver parameters struct as argument
  template <class T, class U>
  void setopts(T &params, int nx, int nz, const U &user_params)
  {
    setopts_hlpr(params, user_params);
    params.di = (X / si::metres) / (nx-1); 
    params.dj = (Z / si::metres) / (nz-1);
    params.dz = params.dj;
  }
  template <class T, class U>
  void setopts(T &params, int nx, int ny, int nz, const U &user_params)
  {
    setopts_hlpr(params, user_params);
    params.di = (X / si::metres) / (nx-1); 
    params.dj = (Y / si::metres) / (ny-1);
    params.dk = (Z / si::metres) / (nz-1);
    params.dz = params.dk;
  }


  template <class concurr_t, class index_t>
  void intcond_hlpr(concurr_t &solver, arr_1D_t &rhod, int rng_seed, index_t index)
  {
    using ix = typename concurr_t::solver_t::ix;
    int nz = solver.advectee().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
    real_t dz = (Z / si::metres) / (nz-1); 

    solver.advectee(ix::rv) = r_t()(index * dz); 
    solver.advectee(ix::u)= setup::u()(index * dz);
    solver.advectee(ix::w) = 0;  
   
    // absorbers
    solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);
    solver.vab_relaxed_state(0) = solver.advectee(ix::u);
    solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state

    // density profile
    solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array

    // initial potential temperature
    solver.advectee(ix::th) = th_dry_fctr()(index * dz); 
    // randomly prtrb tht
    std::mt19937 gen(rng_seed);
    std::uniform_real_distribution<> dis(-0.1, 0.1);
    auto rand = std::bind(dis, gen);

    decltype(solver.advectee(ix::th)) prtrb(solver.advectee(ix::th).shape()); // array to store perturbation
    std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
    solver.advectee(ix::th) += prtrb;
  }

  // function enforcing cyclic values in horizontal directions
  // 2D version
  template<int nd, class arr_t>
  void make_cyclic(arr_t arr,
    typename std::enable_if<nd == 2>::type* = 0)
  { arr(arr.extent(0) - 1, blitz::Range::all()) = arr(0, blitz::Range::all()); }

  // 3D version
  template<int nd, class arr_t>
  void make_cyclic(arr_t arr,
    typename std::enable_if<nd == 3>::type* = 0)
  { 
    arr(arr.extent(0) - 1, blitz::Range::all(), blitz::Range::all()) = 
      arr(0, blitz::Range::all(), blitz::Range::all()); 
    arr(blitz::Range::all(), arr.extent(1) - 1, blitz::Range::all()) = 
      arr(blitz::Range::all(), 0, blitz::Range::all());
  }

  // function expecting a libmpdata++ solver as argument
  // 2D version
  template <int nd, class concurr_t>
  void intcond(concurr_t &solver, arr_1D_t &rhod, int rng_seed,
    typename std::enable_if<nd == 2>::type* = 0
  )
  {
    blitz::secondIndex k;
    intcond_hlpr(solver, rhod, rng_seed, k);
    using ix = typename concurr_t::solver_t::ix;
    make_cyclic<2>(solver.advectee(ix::th));
  }

  // 3D version
  template <int nd, class concurr_t>
  void intcond(concurr_t &solver, arr_1D_t &rhod, int rng_seed,
    typename std::enable_if<nd == 3>::type* = 0
  )
  {
    blitz::thirdIndex k;
    intcond_hlpr(solver, rhod, rng_seed, k);
    using ix = typename concurr_t::solver_t::ix;
    make_cyclic<3>(solver.advectee(ix::th));

    int nz = solver.advectee().extent(ix::w);
    real_t dz = (Z / si::metres) / (nz-1); 

    solver.advectee(ix::v)= setup::v()(k * dz);
    solver.vab_relaxed_state(1) = solver.advectee(ix::v);
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
  // alse set w_LS and hgt_fctrs
  // like in Wojtek's BabyEulag
  template<class user_params_t>
  void env_prof(arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &th_ref, arr_1D_t &rhod, arr_1D_t &w_LS, arr_1D_t &hgt_fctr_vctr, arr_1D_t &hgt_fctr_sclr, int nz, const user_params_t &user_params)
  {
    using libcloudphxx::common::moist_air::R_d_over_c_pd;
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::moist_air::R_d;
    using libcloudphxx::common::const_cp::l_tri;
    using libcloudphxx::common::theta_std::p_1000;

    // pressure profile
    arr_1D_t pre(nz);
    // temperature profile
    arr_1D_t T(nz);
    setup::real_t dz = (Z / si::metres) / (nz-1);

    r_t rt;
    T(0) = th_l(0.) / si::kelvins *  pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    pre(0) = setup::p_0 / si::pascals;
    th_e(0) = th_l(0.) / si::kelvins;
    rv_e(0) = rt(0.);

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
      setup::real_t bottom = R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms * T(k-1) * (1 + 0.61 * rv_e(k-1));
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
      rv_e(k) = rt(k*dz) - delta;
      th_e(k) = th_l(k*dz) / si::kelvins + c * thetme * delta;
      T(k) = th_e(k) * pow(pre(k) / (p_1000<setup::real_t>() / si::pascals),  f);
    }

    // compute reference state theta and rhod
    blitz::firstIndex k;
    // calculate average stability
    blitz::Range notopbot(1, nz-2);
    arr_1D_t st(nz);
    st=0;
    st(notopbot) = (th_e(notopbot+1) - th_e(notopbot-1)) / th_e(notopbot);
    setup::real_t st_avg = blitz::sum(st) / (nz-2) / (2.*dz);
    // reference theta
    th_ref = th_e(0) * exp(st_avg * k * dz);
    // virtual temp at surface
    using libcloudphxx::common::moist_air::R_d_over_c_pd;
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::moist_air::R_d;
    using libcloudphxx::common::theta_std::p_1000;

    setup::real_t T_surf = th_e(0) *  pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    setup::real_t T_virt_surf = T_surf * (1. + 0.608 * rv_e(0));
    setup::real_t rho_surf = (setup::p_0 / si::pascals) / T_virt_surf / 287. ; // TODO: R_d instead of 287
    setup::real_t cs = 9.81 / (c_pd<setup::real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / T_surf;
    // rhod profile
    rhod = rho_surf * exp(- st_avg * k * dz) * pow(
             1. - cs * (1 - exp(- st_avg * k * dz)), (1. / R_d_over_c_pd<setup::real_t>()) - 1);

    // subsidence rate
    w_LS = setup::w_LS_fctr()(k * dz);

    // surface sources relaxation factors
    // for vectors
    real_t z_0 = setup::z_rlx_vctr / si::metres;
    hgt_fctr_vctr = exp(- (k-0.5) * dz / z_0); // z=0 at k=1/2
    hgt_fctr_vctr(0) = 1;
    // for scalars
    z_0 = user_params.z_rlx_sclr;
    hgt_fctr_sclr = exp(- (k-0.5) * dz / z_0);
    hgt_fctr_sclr(0) = 1;
  }
};

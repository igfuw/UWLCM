#pragma once
#include "common.hpp"

namespace setup 
{
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;
  namespace moist_air = libcloudphxx::common::moist_air;
  namespace const_cp = libcloudphxx::common::const_cp;


  using libcloudphxx::common::theta_std::p_1000;
  using libcloudphxx::common::moist_air::R_d_over_c_pd;
  using libcloudphxx::common::moist_air::c_pd;
  using libcloudphxx::common::moist_air::R_d;
  using libcloudphxx::common::moist_air::R_v;
  using libcloudphxx::common::const_cp::l_tri;
  using libcloudphxx::common::const_cp::p_vs;
  using libcloudphxx::common::theta_std::p_1000;
  namespace theta_dry = libcloudphxx::common::theta_dry;

  // RH T and p to rv
  quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
  {
    return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
  }

  const quantity<si::temperature, real_t>
    T_0(283. * si::kelvins);  // surface temperature
  const quantity<si::pressure, real_t> 
    p_0 = 85000 * si::pascals;
  const real_t stab = 1.3e-5; // stability, 1/m
  const real_t env_RH = 0.2;
  const real_t prtrb_RH = 1.; //effective value, should be 1.00, but it caused RH in libcloud = 1.01 in the perturbation; TODO: fix it, its caused by wrong initial condition not taking into account rho/rhod differences?
  // theta (std) at surface
  const quantity<si::temperature, real_t> th_0 = T_0 / pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
  const quantity<si::dimensionless, real_t> rv_0(RH_T_p_to_rv(env_RH, T_0, p_0));
//  const quantity<si::temperature, real_t> th_0_dry = theta_dry::std2dry<real_t>(th_0, rv_0);

  const quantity<si::length, real_t> 
    z_0  = 0    * si::metres,
    Z    = 2400 * si::metres, // DYCOMS: 1500
    X    = 3600 * si::metres, // DYCOMS: 6400
    Y    = 3600 * si::metres, // DYCOMS: 6400
    z_prtrb = 800 * si::metres;

  const real_t z_abs = 125000; // [m] height above which absorber works, no absorber

  // T(th_std, p)
  /*
  template <class real_t>
  quantity<si::temperature, real_t> th2T(const quantity<si::temperature, real_t> &th, const quantity<si::pressure, real_t> &p)
  {
    quantity<si::temperature, real_t> T = th * pow(p / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    return T;
  }
*/
  // some more constants copied from env_prof, todo
//  const setup::real_t T_surf = th2T(th_0, p_0) / si::kelvins;
  const setup::real_t rhod_surf = (setup::p_0 / si::pascals) / (T_0 / si::kelvins) /( R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms + rv_0 * R_v<setup::real_t>() / si::joules * si::kelvins * si::kilograms);
  const setup::real_t rhod_surfW = (setup::p_0 / si::pascals) / (T_0 / si::kelvins) /( R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms);
//  const setup::real_t T_virt_surf = (T_0 / si::kelvins) * (1. + 0.608 * (rv_0 / 1. + rv_0)); // T_virt, i.e. with specific humudity
  //const setup::real_t rho_surf = (setup::p_0 / si::pascals) / T_virt_surf / (R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms); 
  const setup::real_t cs = (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared) / (c_pd<setup::real_t>() / si::joules * si::kilograms * si::kelvins) / stab / (T_0 / si::kelvins);


  struct th_std_fctr
  {
    const real_t th_surf;
    th_std_fctr(const real_t th = (th_0 / si::kelvins)) :
      th_surf(th) {}

    real_t operator()(const real_t &z) const
    {
      return (th_surf) * real_t(exp(stab * z));
    }
    BZ_DECLARE_FUNCTOR(th_std_fctr);
  };

  // temperature profile (constant stability atmosphere, Clark Farley 1984)
  real_t T(const real_t &z)
  {
    return (T_0 / si::kelvins) / exp(- stab * z) * (
             1. - cs * (1 - exp(- stab * z)));

//    return (th_std(z) / si::kelvins) * (1. - cs * (1 - exp(- stab * z)));
  }

  // pressure profile (constant stability atmosphere, Clark Farley 1984), dry...
  real_t p(const real_t &z)
  {
    return p_0 / si::pascals * pow( 1. - cs * (1 - exp(- stab * z)), 1. / R_d_over_c_pd<setup::real_t>() );
  }

  // air density profile (constant stability atmosphere, Clark Farley 1984), dry...
  struct rho_fctr
  {
    const real_t rh_surf;
    rho_fctr(const real_t rho = rhod_surf) :
      rh_surf(rho) {}

    real_t operator()(const real_t &z) const
    {
      return rh_surf * exp(- stab * z) * pow(
               1. - cs * (1 - exp(- stab * z)), (1. / R_d_over_c_pd<setup::real_t>()) - 1);
    }
    BZ_DECLARE_FUNCTOR(rho_fctr);
  };


  // rv(RH, th_dry, rhod)
  real_t RH_th_rhod_to_rv(const real_t &RH, const real_t &th_dry, const real_t &rhod)
  {
    real_t T = theta_dry::T(th_dry * si::kelvins, (rhod * si::kilograms / si::cubic_metres)) / si::kelvins;
    return (p_vs(T * si::kelvins) / si::pascals) * RH / (rhod * T * (R_v<real_t>() / si::joules * si::kilograms * si::kelvins));
  }

  // rv(RH, T, rhod)
  real_t RH_T_rhod_to_rv(const real_t &RH, const real_t &T, const real_t &rhod)
  {
    return (p_vs(T * si::kelvins) / si::pascals) * RH / (rhod * T * (R_v<real_t>() / si::joules * si::kilograms * si::kelvins));
  }

  struct RH
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &x, const real_t &z) const
    {
      real_t r = sqrt( pow( x - (X / si::metres / 2.), 2) + pow( z - (z_prtrb / si::metres), 2));
      if(r <= 200.)
        return prtrb_RH;
      else if(r >= 300.)
        return env_RH;
      else // transition layer
        return env_RH + (prtrb_RH - env_RH) * pow( cos(boost::math::constants::pi<real_t>() / 2. * (r - 200) / 100.), 2);
    }
  BZ_DECLARE_FUNCTOR2(RH);
  };

  struct env_rv
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
  //    return RH_T_rhod_to_rv(env_RH, T(z), rhod_fctr()(z));
      return RH_th_rhod_to_rv(env_RH, th_std_fctr()(z) , rho_fctr()(z));
  //    return RH_T_p_to_rv(env_RH, T(z) * si::kelvins, p(z) * si::pascals);
    }
  BZ_DECLARE_FUNCTOR(env_rv);
  };

  struct prtrb_rv
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &x, const real_t &z) const
    {
  //    return RH_T_rhod_to_rv(env_RH, T(z), rhod_fctr()(z));
      return RH_th_rhod_to_rv(RH()(x,z), th_std_fctr()(z) ,rho_fctr()(z));
  //    return RH_T_p_to_rv(env_RH, T(z) * si::kelvins, p(z) * si::pascals);
    }
  BZ_DECLARE_FUNCTOR2(prtrb_rv);
  };


  // dry air density profile 
  /*
  struct rhod_fctr
  {
    real_t operator()(const real_t &z) const
    {
      return libcloudphxx::common::theta_std::rhod<real_t>(p(z) * si::pascals, th_std(z), env_rv()(z)) / si::kilograms * si::cubic_metres;
    }
    BZ_DECLARE_FUNCTOR(rhod_fctr);
  };
*/

/*
  struct prtrb_rv
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &x, const real_t &z) const
    {
      real_t r = sqrt( pow( x - (X / si::metres / 2.), 2) + pow( z - (z_prtrb / si::metres), 2));
      if(r <= 200.)
        return RH_T_p_to_rv(prtrb_RH, T(z) * si::kelvins, p(z) * si::pascals);
      else if(r >= 300.)
        return RH_T_p_to_rv(env_RH, T(z) * si::kelvins, p(z) * si::pascals);
      else // transition layer
      {
        real_t RH = env_RH + (prtrb_RH - env_RH) * pow( cos(boost::math::constants::pi<real_t>() / 2. * (r - 200) / 100.), 2);
        return RH_T_p_to_rv(RH, T(z) * si::kelvins, p(z) * si::pascals);
      }
    }
  BZ_DECLARE_FUNCTOR2(prtrb_rv);
  };
*/

  // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
  /*
  struct th_dry_fctr
  {
    real_t operator()(const real_t &z) const
    {
      return theta_dry::std2dry<real_t>(th_std(z), env_rv()(z)) / si::kelvins;
    }
    BZ_DECLARE_FUNCTOR(th_dry_fctr);
  };
*/

  template <class T, class U>
  void setopts_hlpr(T &params, const U &user_params)
  {
    params.outdir = user_params.outdir;
    params.outfreq = user_params.outfreq;
    params.spinup = user_params.spinup;
    params.w_src = user_params.w_src;
    params.uv_src = user_params.uv_src;
    params.th_src = user_params.th_src;
    params.rv_src = user_params.rv_src;
    params.prs_tol=1e-6;
    params.dt = user_params.dt;
    params.nt = user_params.nt;
    params.buoyancy_wet = true;
    params.subsidence = false;
    params.friction = false;
//    params.n_iters=1;
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
  void intcond_hlpr(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed, index_t index)
  {
    using ix = typename concurr_t::solver_t::ix;
    int nz = solver.advectee().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
    real_t dz = (Z / si::metres) / (nz-1); 
    int nx = solver.advectee().extent(0);  // ix::w is the index of vertical domension both in 2D and 3D
    real_t dx = (X / si::metres) / (nx-1); 

//    solver.advectee(ix::rv) = rv_e(index);
    solver.advectee(ix::rv) = prtrb_rv()(blitz::tensor::i * dx, blitz::tensor::j * dz); 
/*
    for(int x=0; x<nx; ++x)
      for(int z=0; z<nz; ++z)
      {
         if(RH()(x * dx, z * dz) == env_RH)
           solver.advectee(ix::rv)(x,z) = rv_e(z);
         else
           solver.advectee(ix::rv)(x,z) = RH_th_rhod_to_rv(RH()(x * dx, z * dz), th_e(z) ,rhod(z));
      }
  */  
//solver.advectee(ix::rv)(0,0) = rv_0;
//    solver.advectee(ix::rv) = env_rv()(blitz::tensor::j * dz); 
    solver.advectee(ix::u) = 0;
    solver.advectee(ix::w) = 0;  
   
    // absorbers
    solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);
    solver.vab_relaxed_state(0) = 0;
    solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state

    // density profile
    solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
//    solver.g_factor() = rhod_fctr()(blitz::tensor::j * dz);

    // initial potential temperature
//    solver.advectee(ix::th) = th_std_fctr()(index * dz); 
    solver.advectee(ix::th) = th_e(index); 
//solver.advectee(ix::th)(0,0) = th_0_dry / si::kelvins;  
  }

  // function expecting a libmpdata++ solver as argument
  // 2D version
  template <int nd, class concurr_t>
  void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed,
    typename std::enable_if<nd == 2>::type* = 0
  )
  {
    blitz::secondIndex k;
    intcond_hlpr(solver, rhod, th_e, rv_e, rng_seed, k);
    using ix = typename concurr_t::solver_t::ix;
//    make_cyclic<2>(solver.advectee(ix::th));
  }

  // 3D version
  template <int nd, class concurr_t>
  void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed,
    typename std::enable_if<nd == 3>::type* = 0
  )
  {
    blitz::thirdIndex k;
    intcond_hlpr(solver, rhod, th_e, rv_e, rng_seed, k);
    using ix = typename concurr_t::solver_t::ix;
  //  make_cyclic<3>(solver.advectee(ix::th));

    int nz = solver.advectee().extent(ix::w);
    real_t dz = (Z / si::metres) / (nz-1); 

    solver.advectee(ix::v) = 0;
    solver.vab_relaxed_state(1) = 0;
  }

  // calculate the initial environmental theta and rv profiles as Wojtek does it
  template<class user_params_t>
  void env_prof(arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &th_ref, arr_1D_t &pre_ref, arr_1D_t &rhod, arr_1D_t &w_LS, arr_1D_t &hgt_fctr_vctr, arr_1D_t &hgt_fctr_sclr, int nz, const user_params_t &user_params)
  // pre_ref - total pressure
  // th_e - dry potential temp
  // th_ref - dry potential temp refrence profile
  // rhod - dry density profile
  {
    setup::real_t dz = (Z / si::metres) / (nz-1);
    using libcloudphxx::common::moist_air::R_d_over_c_pd;
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::moist_air::R_d;
    using libcloudphxx::common::const_cp::l_tri;
    using libcloudphxx::common::theta_std::p_1000;

    using setup::real_t;
    blitz::firstIndex k;
    // temperature profile
    arr_1D_t T(nz);

    setup::real_t tt0 = 273.17;
    setup::real_t rv = 461;
    setup::real_t ee0 = 611.;
    const real_t gg = 9.81;
    real_t rg = R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms;
    setup::real_t a = R_d<setup::real_t>() / rv / si::joules * si::kelvins * si::kilograms;
//    setup::real_t b = l_tri<setup::real_t>() / si::joules * si::kilograms / rv / tt0;
    setup::real_t c = l_tri<setup::real_t>() / c_pd<setup::real_t>() / si::kelvins;
    setup::real_t d = l_tri<setup::real_t>() / si::joules * si::kilograms / rv;
    setup::real_t cap = R_d_over_c_pd<setup::real_t>(); 
    real_t capi = 1./cap;

    // surface data
    
    real_t tt = T_0 / si::kelvins; // T(0)
    real_t delt = (tt - tt0) / (tt * tt0); 
    real_t esw = ee0*exp(d * delt);
    real_t qvs = a * esw / ((p_0 / si::pascals) -esw);
    rv_e(0) = env_RH * qvs;
    real_t th_e_surf = th_0 / si::kelvins * (1 + a * rv_e(0)); // virtual potential temp
    
    th_e = th_std_fctr(th_e_surf)(k * dz);
    
    pre_ref(0.) = p_0 / si::pascals;
    T(0.) = T_0 / si::kelvins;
    
    for(int k=1; k<nz; ++k)
    {
      real_t zz = k * dz;  
      // predictor
       real_t rhob=pre_ref(k-1) / rg / (T(k-1)*(1.+a*rv_e(k-1)));
       pre_ref(k)=pre_ref(k-1) - gg*rhob*dz;
// iteration for T and qv:
       rv_e(k)=rv_e(k-1);
       T(k)=th_e(k)* pow(pre_ref(k)/1.e5, cap); 
       T(k)=T(k)/(1.+a*rv_e(k));
      
      for(int iter=0; iter<4; ++iter)
      {
        tt=T(k);
        delt=(tt-tt0)/(tt*tt0);
        esw=ee0*exp(d * delt);
        qvs=a * esw /(pre_ref(k)-esw);
        rv_e(k)=env_RH*qvs;
       T(k)=th_e(k)* pow(pre_ref(k)/1.e5, cap);
        T(k)=T(k)/(1.+a*rv_e(k));
      }

      // corrector
       real_t rhon=pre_ref(k) / rg / (T(k)*(1.+a*rv_e(k)));
       pre_ref(k)=pre_ref(k-1) - gg*(rhob+rhon) / 2. *dz;
// iteration for T and qv:
       T(k)=th_e(k)* pow(pre_ref(k)/1.e5, cap);
       T(k)=T(k)/(1.+a*rv_e(k));
      
      for(int iter=0; iter<4; ++iter)
      {
        tt=T(k);
        delt=(tt-tt0)/(tt*tt0);
        esw=ee0*exp(d * delt);
        qvs=a * esw /(pre_ref(k)-esw);
        rv_e(k)=env_RH*qvs;
       T(k)=th_e(k)* pow(pre_ref(k)/1.e5, cap);
        T(k)=T(k)/(1.+a*rv_e(k));
      }

    }

    th_ref = th_std_fctr()(k * dz);
//    rhod = rho_fctr(rhod_surfW)(k * dz); // this way rhod is total density profile, not dry
    rhod = rho_fctr(rhod_surf)(k * dz); // rhod is dry density profile?

    std::cout << "th_v_e: " << th_e << std::endl;
    std::cout << "rv_e: " << rv_e << std::endl;
    std::cout << "T_e: " << T << std::endl;
    std::cout << "th ref: " << th_ref << std::endl;
    std::cout << "rho_ref: " << rhod << std::endl;

    // turn virtual potential temperature env profile into env profile of standard potential temp
    th_e = th_e / (1. + a * rv_e);
    std::cout << "th_e: " << th_e << std::endl;


    // calc rhod using our formulas...; gives the same rhod
    /*
    for(int k=1; k<nz; ++k)
    {
      quantity<si::dimensionless, real_t> si_rv_e = rv_e(k);
      rhod(k) = libcloudphxx::common::theta_std::rhod(pre_ref(k) * si::pascals, th_e(k) * si::kelvins, si_rv_e) / si::kilograms * si::cubic_metres; 
    }
    std::cout << "rhod(p,th_std,rv_e): " << rhod << std::endl;
    */

    // turn standard potential temp into dry potential temp
    for(int k=1; k<nz; ++k)
    {
      quantity<si::dimensionless, real_t> si_rv_e = rv_e(k);
      th_e(k) = libcloudphxx::common::theta_dry::std2dry(th_e(k) * si::kelvins, si_rv_e) / si::kelvins; 
    }
    std::cout << "th_e_dry: " << th_e << std::endl;

    
    // adjust rv_e according to prtrb_rv later...
    /*
    for(int z=0; z<nz; ++z)
       rv_e(z) = RH_th_rhod_to_rv(env_RH, th_e(z) ,rhod(z));
    std::cout << "rv_e: " << rv_e << std::endl;
*/
  }

  // calculate the initial environmental theta and rv profiles
   /*
  template<class user_params_t>
  void env_prof(arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &th_ref, arr_1D_t &pre_ref, arr_1D_t &rhod, arr_1D_t &w_LS, arr_1D_t &hgt_fctr_vctr, arr_1D_t &hgt_fctr_sclr, int nz, const user_params_t &user_params)
  {
    setup::real_t dz = (Z / si::metres) / (nz-1);
    blitz::firstIndex k;
    th_ref = th_std_fctr()(k * dz);
    th_e = th_ref;
    rv_e = env_rv()(k * dz);
    rhod = rho_fctr()(k * dz);

    std::cout << "th ref: " << th_ref << std::endl;
    std::cout << "th e: " << th_e << std::endl;
    std::cout << "rv e: " << rv_e << std::endl;
    std::cout << "rho_ref: " << rhod << std::endl;

    // subsidence rate
    //w_LS = setup::w_LS_fctr()(k * dz);

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
*/
};

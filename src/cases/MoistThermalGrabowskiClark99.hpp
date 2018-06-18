// In fact its the setup from the GMD 2017 paper Grabowski, Dziekan, Pawlowska on Twomey activation of SDs,
// it's very similar to the Grabowski Clark 99 setup

#pragma once
#include "CasesCommon.hpp"

namespace setup 
{
  namespace moist_thermal
  {
    namespace moist_air = libcloudphxx::common::moist_air;
    namespace const_cp = libcloudphxx::common::const_cp;
    namespace theta_std = libcloudphxx::common::theta_std;

    using libcloudphxx::common::theta_std::p_1000;
    using libcloudphxx::common::moist_air::R_d_over_c_pd;
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::moist_air::R_d;
    using libcloudphxx::common::moist_air::R_v;
    using libcloudphxx::common::const_cp::l_tri;
    using libcloudphxx::common::const_cp::p_vs;
    using libcloudphxx::common::theta_std::p_1000;


    // RH T and p to rv
    quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
    }
    
    const quantity<si::temperature, real_t> T_0(283. * si::kelvins);  // surface temperature
    const quantity<si::pressure, real_t> p_0(85000 * si::pascals); // total surface temperature
    const real_t stab = 1.3e-5; // stability, 1/m
    const real_t env_RH = 0.2;
    const real_t prtrb_RH = 1. - 1e-10;
    const quantity<si::temperature, real_t> th_std_0 = T_0 / pow(p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    const quantity<si::dimensionless, real_t> rv_0 = RH_T_p_to_rv(env_RH, T_0, p_0);
    const quantity<si::dimensionless, real_t> qv_0 = rv_0 / (1. + rv_0); // specific humidity at surface
    const quantity<si::length, real_t> 
//     z_0  ( 0    * si::metres),
     Z    ( 2400 * si::metres), // DYCOMS: 1500
     X    ( 3600 * si::metres), // DYCOMS: 6400
     Y    ( 3600 * si::metres), // DYCOMS: 6400
     z_prtrb ( 800 * si::metres);
//    const setup::real_t rhod_surf = (p_0 / si::pascals) / (T_0 / si::kelvins) /( R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms + rv_0 * R_v<setup::real_t>() / si::joules * si::kelvins * si::kilograms);
 //   const setup::real_t rhod_surf = (p_0 / si::pascals) / (T_0 / si::kelvins) /( R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms) / (1+0.608 * qv_0);
    const setup::real_t rhod_surf = theta_std::rhod(p_0, th_std_0, rv_0) * si::cubic_metres / si::kilograms;
    //const setup::real_t rhod_surfW = (p_0 / si::pascals) / (T_0 / si::kelvins) /( R_d<setup::real_t>() / si::joules * si::kelvins * si::kilograms);
    const setup::real_t cs = (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared) / (c_pd<setup::real_t>() / si::joules * si::kilograms * si::kelvins) / stab / (T_0 / si::kelvins);

    const real_t z_abs = 125000; // [m] height above which absorber works, no absorber

///WARNING: these functors, taken from Clark Farley 1984, are for dry air!!

    struct th_std_fctr
    {
      const real_t th_surf;
      th_std_fctr(const real_t th) :
        th_surf(th) {}
    
      real_t operator()(const real_t &z) const
      {
        return (th_surf) * real_t(exp(stab * z));
      }
      BZ_DECLARE_FUNCTOR(th_std_fctr);
    };
  
/*  
    // temperature profile (constant stability atmosphere, Clark Farley 1984)
    real_t T(const real_t &z)
    {
      return (T_0 / si::kelvins) / exp(- stab * z) * (
               1. - cs * (1 - exp(- stab * z)));
    }
    
    // pressure profile (constant stability atmosphere, Clark Farley 1984), dry...
    real_t p(const real_t &z)
    {
      return p_0 / si::pascals * pow( 1. - cs * (1 - exp(- stab * z)), 1. / R_d_over_c_pd<setup::real_t>() );
    }
    
*/
    // air density profile (constant stability atmosphere, Clark Farley 1984), dry...
    struct rho_fctr
    {
      const real_t rh_surf;
      rho_fctr(const real_t &rho) :
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
      quantity<si::dimensionless, real_t> operator()(const real_t &r) const // r - distance from the center of the perturbation
      {
        if(r <= 250.)
          return prtrb_RH;
        else if(r >= 350.)
          return env_RH;
        else // transition layer
          return env_RH + (prtrb_RH - env_RH) * pow( cos(boost::math::constants::pi<real_t>() / 2. * (r - 250) / 100.), 2);
      }
    BZ_DECLARE_FUNCTOR(RH);
    };
    
/*
    struct env_rv
    {
      quantity<si::dimensionless, real_t> operator()(const real_t &z) const
      {
    //    return RH_T_rhod_to_rv(env_RH, T(z), rhod_fctr()(z));
        return RH_th_rhod_to_rv(env_RH, th_std_fctr(th_std_0 / si::kelvins)(z) , rho_fctr(rhod_surf)(z));
    //    return RH_T_p_to_rv(env_RH, T(z) * si::kelvins, p(z) * si::pascals);
      }
    BZ_DECLARE_FUNCTOR(env_rv);
    };
  */  

    struct prtrb_rv 
    {
      arr_1D_t &_th_std, &_rhod;
      real_t dz;
      prtrb_rv(arr_1D_t _th_std, arr_1D_t _rhod, real_t dz): _th_std(_th_std), _rhod(_rhod), dz(dz) {}

      quantity<si::dimensionless, real_t> operator()(const real_t &r, const real_t &z) const
      {
        return RH_th_rhod_to_rv(RH()(r), this->_th_std(z/this->dz) , this->_rhod(z/this->dz));
    //    return RH_T_p_to_rv(env_RH, T(z) * si::kelvins, p(z) * si::pascals);
    //    return RH_T_rhod_to_rv(env_RH, T(z), rhod_fctr()(z));
      }
    BZ_DECLARE_FUNCTOR2(prtrb_rv);
    };

    // its in fact the moist thermal from our 2017 GMD paper on Twomey SDs? differences: kappa=1.28, i.e. sea salt aerosol
    template<class concurr_t>
    class MoistThermalGrabowskiClark99 : public CasesCommon<concurr_t>
    {

      protected:
    
      void setopts_hlpr(typename concurr_t::solver_t::rt_params_t &params, const user_params_t &user_params)
      {
        params.outdir = user_params.outdir;
        params.outfreq = user_params.outfreq;
        params.spinup = user_params.spinup;
        params.w_src = true;
        params.uv_src = false;
        params.th_src = false;
        params.rv_src = false;
        params.dt = user_params.dt;
        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = false;
        params.friction = false;
    //    params.n_iters=1;
      }
    
      template <class index_t>
      void intcond_hlpr(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed, index_t index)
      {
        using ix = typename concurr_t::solver_t::ix;
        int nz = solver.advectee().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (Z / si::metres) / (nz-1); 
        int nx = solver.advectee().extent(0);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dx = (X / si::metres) / (nx-1); 
    
       // solver.advectee(ix::rv) = rv_e(index);

        solver.advectee(ix::u) = 0;
        solver.advectee(ix::w) = 0;  
       
        // absorbers
        solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);
        solver.vab_relaxed_state(0) = 0;
        solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
    
        // density profile
        solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
    
        // initial potential temperature
        solver.advectee(ix::th) = th_e(index); 
//libcloudphxx::common::theta_dry::std2dry(th_std_0, rv_0) / si::kelvins; 
      }
    
    
      public:
      // calculate the initial environmental theta and rv profiles as Wojtek does it
      // i.e. for stable virtual standard potential temperature
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
        //rv_e(0) = env_RH * qvs;
        rv_e(0) = rv_0;// env_RH * qvs;
        real_t th_e_surf = th_std_0 / si::kelvins * (1 + a * rv_e(0)); // virtual potential temp
        
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
          rv_e(k) =  RH_T_p_to_rv(env_RH, T(k) * si::kelvins, pre_ref(k) * si::pascals); // cheating!
    
        }
    
        //th_ref = th_std_fctr(th_std_0 / si::kelvins)(k * dz);
        rhod = rho_fctr(rhod_surf)(k * dz); // rhod is dry density profile?
    
        // turn virtual potential temperature env profile into env profile of standard potential temp
        th_e = th_e / (1. + a * rv_e);
    
        // turn standard potential temp into dry potential temp
        for(int k=0; k<nz; ++k)
        {
          quantity<si::dimensionless, real_t> si_rv_e = rv_e(k);
          th_e(k) = libcloudphxx::common::theta_dry::std2dry(th_e(k) * si::kelvins, si_rv_e) / si::kelvins; 
        }
        th_ref = th_e;//th_std_fctr(th_std_0 / si::kelvins)(k * dz);
      }

      // ctor
      MoistThermalGrabowskiClark99()
      {
        this->kappa = 1.28; // NaCl aerosol
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

    // 2d/3d children
    template<class concurr_t>
    class MoistThermalGrabowskiClark99_2d : public MoistThermalGrabowskiClark99<concurr_t>
    {
      public:
      // function expecting a libmpdata solver parameters struct as argument
      void setopts(typename concurr_t::solver_t::rt_params_t &params, int nx, int nz, const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nx-1); 
        params.dj = (Z / si::metres) / (nz-1);
        params.dz = params.dj;
      }

      // function expecting a libmpdata++ solver as argument
      void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(solver, rhod, th_e, rv_e, rng_seed, k);

        using ix = typename concurr_t::solver_t::ix;
        int nz = solver.advectee().extent(ix::w); 
        real_t dz = (Z / si::metres) / (nz-1); 
        int nx = solver.advectee().extent(0); 
        real_t dx = (X / si::metres) / (nx-1); 
        solver.advectee(ix::rv) = prtrb_rv(th_e, rhod, dz)(
          sqrt(
            pow(blitz::tensor::i * dx - (X / si::metres / 2.), 2) + 
            pow(blitz::tensor::j * dz - (z_prtrb / si::metres), 2)
          ),
          blitz::tensor::j * dz
        );

      }
    };

    template<class concurr_t>
    class MoistThermalGrabowskiClark99_3d : public MoistThermalGrabowskiClark99<concurr_t>
    {
      public:
      // function expecting a libmpdata solver parameters struct as argument
      void setopts(typename concurr_t::solver_t::rt_params_t &params, int nx, int ny, int nz, const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nx-1); 
        params.dj = (Y / si::metres) / (ny-1);
        params.dk = (Z / si::metres) / (nz-1);
        params.dz = params.dk;
      }

      // function expecting a libmpdata++ solver as argument
      void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(solver, rhod, th_e, rv_e, rng_seed, k);
        using ix = typename concurr_t::solver_t::ix;
        int nz = solver.advectee().extent(2); 
        real_t dz = (Z / si::metres) / (nz-1); 
        int nx = solver.advectee().extent(0); 
        real_t dx = (X / si::metres) / (nx-1); 
        int ny = solver.advectee().extent(1); 
        real_t dy = (Y / si::metres) / (ny-1); 
        solver.advectee(ix::rv) = prtrb_rv(th_e, rhod, dz)(
          sqrt(
            pow(blitz::tensor::i * dx - (X / si::metres / 2.), 2) + 
            pow(blitz::tensor::j * dy - (Y / si::metres / 2.), 2) + 
            pow(blitz::tensor::k * dz - (z_prtrb / si::metres), 2)
          ),
          blitz::tensor::k * dz
        );
    
        solver.advectee(ix::v) = 0;
        solver.vab_relaxed_state(1) = 0;
      }

/*
      // TODO: make it work in 3d?
      MoistThermalGrabowskiClark99_3d()
      {
        throw std::runtime_error("Moist Thermal doesn't work in 3d yet");
      }
*/
    };
  };
};

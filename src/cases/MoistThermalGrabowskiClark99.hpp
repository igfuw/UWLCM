// In fact its the setup from the GMD 2017 paper Grabowski, Dziekan, Pawlowska on Twomey activation of SDs,
// it's very similar to the Grabowski Clark 99 setup

#pragma once
#include "CasesCommon.hpp"

namespace detail
{
  struct calc_p_v
  {
    setup::real_t operator()(setup::real_t p, setup::real_t rv) const
    {return libcloudphxx::common::moist_air::p_v<setup::real_t>(p * si::pascals, rv)  / si::pascals;}
    BZ_DECLARE_FUNCTOR2(calc_p_v)
  };
};

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


    // RH T and p to rv assuming RH = r_v / r_vs
    inline quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return  RH * const_cp::r_vs<real_t>(T, p);
    }

/*
    // RH T and p to rv assuming RH = p_v / p_vs
    quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
    }

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
*/
    
    const quantity<si::temperature, real_t> T_0(283. * si::kelvins);  // surface temperature
    const quantity<si::pressure, real_t> p_0(85000 * si::pascals); // total surface temperature
    const real_t stab = 1.3e-5; // stability, 1/m
    const real_t env_RH = 0.2;
    const real_t prtrb_RH = 1. - 1e-10;
    const quantity<si::temperature, real_t> th_std_0 = T_0 / pow(p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
    const quantity<si::dimensionless, real_t> rv_0 = RH_T_p_to_rv(env_RH, T_0, p_0);
    const quantity<si::dimensionless, real_t> qv_0 = rv_0 / (1. + rv_0); // specific humidity at surface
    const quantity<si::length, real_t> 
     Z    ( 2400 * si::metres), 
     X    ( 3600 * si::metres), 
     Y    ( 3600 * si::metres), 
     z_prtrb ( 800 * si::metres);
    const setup::real_t rhod_surf = theta_std::rhod(p_0, th_std_0, rv_0) * si::cubic_metres / si::kilograms;
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
    
    struct prtrb_rv 
    {
      arr_1D_t &_T, &_p;
      real_t dz;
      prtrb_rv(arr_1D_t &_T, arr_1D_t &_p, real_t dz): _T(_T), _p(_p), dz(dz) {}

      quantity<si::dimensionless, real_t> operator()(const real_t &r, const real_t &z) const
      {
        return RH_T_p_to_rv(RH()(r), this->_T(z/this->dz) * si::kelvins , this->_p(z/this->dz) * si::pascals);
      }
      BZ_DECLARE_FUNCTOR2(prtrb_rv);
    };

    // its in fact the moist thermal from our 2017 GMD paper on Twomey SDs? differences: kappa=1.28, i.e. sea salt aerosol
    template<class case_ct_params_t, int n_dims>
    class MoistThermalGrabowskiClark99Common : public CasesCommon<case_ct_params_t, n_dims>
    {

      protected:
      using parent_t = CasesCommon<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;
    
      void setopts_hlpr(rt_params_t &params, const user_params_t &user_params)
      {
        params.outdir = user_params.outdir;
        params.outfreq = user_params.outfreq;
        params.spinup = user_params.spinup;
        params.w_src = true;
        params.uv_src = false;
        params.th_src = false;
        params.rv_src = false;
        params.rc_src = false;
        params.rr_src = false;
        params.nc_src = false;
        params.nr_src = false;
        params.dt = user_params.dt;
        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = false;
        params.friction = false;
        params.coriolis = false;
        params.radiation = false;
    //    params.n_iters=1;

        this->setopts_sgs(params);
      }
    
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &solver,
                        arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, int rng_seed, index_t index)
      {
        int nz = solver.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (Z / si::metres) / (nz-1); 
    
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
      }
    
    
      public:
      // calculate the initial environmental theta and rv profiles as Wojtek does it
      // i.e. for stable virtual standard potential temperature
      void set_profs(profiles_t &profs, int nz, const user_params_t &user_params)
      // pre_ref - total pressure
      // th_e - dry potential temp
      // th_ref - dry potential temp refrence profsile
      // rhod - dry density profsile
      {

        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::const_cp::l_tri;
        using libcloudphxx::common::theta_std::p_1000;
        using setup::real_t;

        parent_t::set_profs(profs, nz, user_params);

        real_t dz = (Z / si::metres) / (nz-1);
        blitz::firstIndex k;
        // temperature and total pressure profiles
        arr_1D_t T(nz), pre_ref(nz);
    
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
        profs.rv_e(0) = rv_0;// env_RH * qvs;
        profs.rl_e = 0.;
        real_t th_e_surf = th_std_0 / si::kelvins * (1 + a * profs.rv_e(0)); // virtual potential temp

        profs.th_e = th_std_fctr(th_e_surf)(k * dz);
        
        pre_ref(0.) = p_0 / si::pascals;
        profs.p_e(0) = pre_ref(0);
        T(0.) = T_0 / si::kelvins;
        
        for(int k=1; k<nz; ++k)
        {
          real_t zz = k * dz;  
          // predictor
           real_t rhob=pre_ref(k-1) / rg / (T(k-1)*(1.+a*profs.rv_e(k-1))); // density of air at k-1
           pre_ref(k)=pre_ref(k-1) - gg*rhob*dz; // estimate of pre at k (dp = -g * rho * dz)
    // iteration for T and qv:
           profs.rv_e(k)=profs.rv_e(k-1);
           T(k)=profs.th_e(k)* pow(pre_ref(k)/1.e5, cap); 
           T(k)=T(k)/(1.+a*profs.rv_e(k));
          
          for(int iter=0; iter<4; ++iter)
          {
            tt=T(k);
            delt=(tt-tt0)/(tt*tt0);
            esw=ee0*exp(d * delt);
            qvs=a * esw /(pre_ref(k)-esw);
            profs.rv_e(k)=env_RH*qvs;
           T(k)=profs.th_e(k)* pow(pre_ref(k)/1.e5, cap);
            T(k)=T(k)/(1.+a*profs.rv_e(k));
          }
    
          // corrector
           real_t rhon=pre_ref(k) / rg / (T(k)*(1.+a*profs.rv_e(k)));
           pre_ref(k)=pre_ref(k-1) - gg*(rhob+rhon) / 2. *dz;
    // iteration for T and qv:
           T(k)=profs.th_e(k)* pow(pre_ref(k)/1.e5, cap);
           T(k)=T(k)/(1.+a*profs.rv_e(k));
          
          for(int iter=0; iter<4; ++iter)
          {
            tt=T(k);
            delt=(tt-tt0)/(tt*tt0);
            esw=ee0*exp(d * delt);
            qvs=a * esw /(pre_ref(k)-esw);
            profs.rv_e(k)=env_RH*qvs;
            T(k)=profs.th_e(k)* pow(pre_ref(k)/1.e5, cap);
            T(k)=T(k)/(1.+a*profs.rv_e(k));
          }
          //rv_e(k) =  RH_T_p_to_rv(env_RH, T(k) * si::kelvins, pre_ref(k) * si::pascals); // cheating!
          profs.p_e(k) = pre_ref(k);
        }
    
        //th_ref = th_std_fctr(th_std_0 / si::kelvins)(k * dz);
        profs.rhod = rho_fctr(rhod_surf)(k * dz); // rhod is dry density profsile?
    
        // turn virtual potential temperature env profsile into env profsile of standard potential temp
        profs.th_e = profs.th_e / (1. + a * profs.rv_e);
    
        // turn standard potential temp into dry potential temp
/*
        for(int k=0; k<nz; ++k)
        {
          quantity<si::dimensionless, real_t> si_rv_e = rv_e(k);
          th_e(k) = libcloudphxx::common::theta_dry::std2dry(th_e(k) * si::kelvins, si_rv_e) / si::kelvins; 
          real_t p_d = pre_ref(k) - libcloudphxx::common::moist_air::p_v<real_t>(pre_ref(k) * si::pascals, rv_e(k))  / si::pascals;
        }
*/
        profs.th_ref = profs.th_e;//th_std_fctr(th_std_0 / si::kelvins)(k * dz);

        profs.w_LS = 0.; // no subsidence
        profs.th_LS = 0.; // no large-scale horizontal advection
        profs.rv_LS = 0.; 
      }

      // ctor
      MoistThermalGrabowskiClark99Common()
      {
        this->kappa = 1.28; // NaCl aerosol
        this->Z = Z;
      }
    };

    // 2d/3d children
    template<class case_ct_params_t, int n_dims>
    class MoistThermalGrabowskiClark99;

    template<class case_ct_params_t>
    class MoistThermalGrabowskiClark99<case_ct_params_t, 2> : public MoistThermalGrabowskiClark99Common<case_ct_params_t, 2>
    {
      using parent_t = MoistThermalGrabowskiClark99Common<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // function expecting a libmpdata solver parameters struct as argument
      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      // function expecting a libmpdata++ solver as argument
      void intcond(typename parent_t::concurr_any_t &solver,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(solver, rhod, th_e, rv_e, rl_e, rng_seed, k);

//        arr_1D_t p_d_e(p_e - detail::calc_p_v()(p_e, rv_e));
        arr_1D_t T(th_e * pow(p_e / 1.e5, R_d_over_c_pd<setup::real_t>()));

        int nz = solver.advectee_global().extent(ix::w); 
        real_t dz = (Z / si::metres) / (nz-1); 
        int nx = solver.advectee_global().extent(0); 
        real_t dx = (X / si::metres) / (nx-1); 
        solver.advectee(ix::rv) = prtrb_rv(T, p_e, dz)(
          sqrt(
            pow(blitz::tensor::i * dx - (X / si::metres / 2.), 2) + 
            pow(blitz::tensor::j * dz - (z_prtrb / si::metres), 2)
          ),
          blitz::tensor::j * dz
        );
     
      }

      public:
      MoistThermalGrabowskiClark99()
      {
        this->X = X;
      }
    };

    template<class case_ct_params_t>
    class MoistThermalGrabowskiClark99<case_ct_params_t, 3> : public MoistThermalGrabowskiClark99Common<case_ct_params_t, 3>
    {
      using parent_t = MoistThermalGrabowskiClark99Common<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // function expecting a libmpdata solver parameters struct as argument
      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Y / si::metres) / (nps[1]-1);
        params.dk = (Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      // function expecting a libmpdata++ solver as argument
      void intcond(typename parent_t::concurr_any_t &solver,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(solver, rhod, th_e, rv_e, rl_e, rng_seed, k);

//        arr_1D_t p_d_e(p_e - detail::calc_p_v()(p_e, rv_e));
        arr_1D_t T(th_e * pow(p_e / 1.e5, R_d_over_c_pd<setup::real_t>()));

        int nz = solver.advectee_global().extent(2); 
        real_t dz = (Z / si::metres) / (nz-1); 
        int nx = solver.advectee_global().extent(0); 
        real_t dx = (X / si::metres) / (nx-1); 
        int ny = solver.advectee_global().extent(1); 
        real_t dy = (Y / si::metres) / (ny-1); 
        solver.advectee(ix::rv) = prtrb_rv(T, p_e, dz)(
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

      public:
      MoistThermalGrabowskiClark99()
      {
        this->X = X;
        this->Y = Y;
      }
    };
  };
};

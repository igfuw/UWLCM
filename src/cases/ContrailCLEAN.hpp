// Contrail setup for NASA CLEAN experiment

#pragma once
#include "Anelastic.hpp"

namespace cases 
{
  namespace contrail_CLEAN
  {
    namespace moist_air = libcloudphxx::common::moist_air;
    namespace const_cp = libcloudphxx::common::const_cp;
    namespace theta_std = libcloudphxx::common::theta_std;

    // RH T and p to rv assuming RH = p_v / p_vs
    inline quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
    }
 
    const quantity<si::temperature, real_t> T_0((273.15 - 50) * si::kelvins);  // surface temperature
    const quantity<si::pressure, real_t> p_0(26500 * si::pascals);
    const quantity<si::dimensionless, real_t> RH_0(0.5);
    const quantity<si::length, real_t> 
     Z_def    ( 100 * si::metres), 
     X_def    ( 500 * si::metres), 
     Y_def    ( 500 * si::metres);

    const real_t z_abs = 20; // velocity absorbers working z_abs meters from top bottom edges; 
    // const real_t x_abs = 1000; // velocity absorbers working 1000 meters right edge TODO: left edge as well?

    inline quantity<si::temperature, real_t> th_l_CLEAN(const real_t &z)
    {
      return T_0 / theta_std::exner(p_0);
    }

    inline quantity<si::dimensionless, real_t> r_t_CLEAN(const real_t &z)
    {
      return RH_T_p_to_rv(RH_0, T_0, p_0);
    }


    template<class case_ct_params_t, int n_dims>
    class ContrailCLEANCommon : public Anelastic<case_ct_params_t, n_dims>
    {

      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      real_t spinup_time; // time to turn on large-scale sources after spinup

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return th_l_CLEAN(z);
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return r_t_CLEAN(z);
      }
    
      void setopts_hlpr(rt_params_t &params, const user_params_t &user_params)
      {
        params.w_src = true;
        params.uv_src = false;
        params.th_src = true;
        params.rv_src = true;
        params.rc_src = false;
        params.rr_src = false;
        params.nc_src = false;
        params.nr_src = false;
        // params.buoyancy_wet = false;
        params.buoyancy_wet = true;
        params.subsidence = subs_t::none;
        params.vel_subsidence = false;
        params.friction = false;
        params.coriolis = false;
        params.radiation = false;
    //    params.n_iters=1;

        spinup_time = user_params.spinup * user_params.dt; // save spinup time to turn on large-scale sources after spinup

        this->setopts_sgs(params);
      }
    
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr,
                        arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, int rng_seed, index_t index)
      {
        // int nz = concurr.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        int nz = rhod.extent(0)-1;
        real_t dz = (this->Z / si::metres) / (nz-1); 
    
        // zero
        concurr.advectee(ix::u) = 0;
        concurr.advectee(ix::w) = 0;  
        concurr.vab_relaxed_state(0) = 0; // horizontal velocity (x) relaxed state
        concurr.vab_relaxed_state(ix::w) = 0; // vertical velocity relaxed state
        concurr.vab_coefficient() = 0; // no velocity absorber by default

        // gravity wave absorber
        double gw_abs_mag = 1.;///100; // absorber strength
        std::conditional_t<n_dims==2, blitz::secondIndex, blitz::thirdIndex> k; // vertical index for iterating over vertical levels; works both for 2D and 3D
        // bottom edge
        concurr.vab_coefficient() += where(k * dz <= z_abs,  gw_abs_mag * (z_abs - k * dz) / (z_abs), 0);
        // top
        concurr.vab_coefficient() += where(k * dz >= this->Z / si::metres - z_abs,  gw_abs_mag * (k * dz - (this->Z / si::metres - z_abs)) / (z_abs), 0);
    
        // density profile
        concurr.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
    
        // initial water vapor mixing ratio
        concurr.advectee(ix::rv) = rv_e(index); 
        // initial potential temperature
        concurr.advectee(ix::th) = th_e(index); 

        // randomly prtrb tht and rv
        // TODO: make it a function cause its used in many cases
        // TODO2: dont use advectee_global
        // NOTE: all processes do this, but ultimately only perturbation calculated by MPI rank 0 is used
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-0.1, 0.1);
          auto rand = std::bind(dis, gen);
  
          auto th_global = concurr.advectee_global(ix::th);
          decltype(concurr.advectee(ix::th)) prtrb(th_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          // prtrb = where(index * dz >= 1600., 0., prtrb); // no perturbation above 1.6km
          th_global += prtrb;
          this->make_cyclic(th_global);
          concurr.advectee_global_set(th_global, ix::th);
        }
        {
          std::mt19937 gen(rng_seed+1); // different seed than in th. NOTE: if the same instance of gen is used in th and rv, for some reason it gives the same sequence in rv as in th despite being advanced in th prtrb
          std::uniform_real_distribution<> dis(-0.025e-3, 0.025e-3);
          auto rand = std::bind(dis, gen);
  
          auto rv_global = concurr.advectee_global(ix::rv);
          decltype(concurr.advectee(ix::rv)) prtrb(rv_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          // prtrb = where(index * dz >= 1600., 0., prtrb); // no perturbation above 1.6km
          rv_global += prtrb;
          this->make_cyclic(rv_global);
          concurr.advectee_global_set(rv_global, ix::rv);
        }
      }
    
    
      public:
      // calculate the initial environmental theta and rv profiles as Wojtek does it
      // i.e. for stable virtual standard potential temperature
      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
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

        blitz::firstIndex k;
        real_t dz = (this->Z / si::metres) / (nz-1);

        parent_t::set_profs(profs, nz, user_params);
        parent_t::env_prof(profs, nz);
        parent_t::ref_prof(profs, nz);
  
        // subsidence rate
        profs.w_LS = 0;
      }

      // ctor
      ContrailCLEANCommon(const real_t _X, const real_t _Y, const real_t _Z)
      {
        init();
        // this->kappa = 1.28; // NaCl aerosol

        this->X = _X < 0 ? X_def : _X * si::meters;
        if constexpr (n_dims == 3)
          this->Y = _Y < 0 ? Y_def : _Y * si::meters;
        this->Z = _Z < 0 ? Z_def : _Z * si::meters;
      }

      void init()
      {
        this->p_0 = p_0;
      }
    };

    // 2d/3d children
    template<class case_ct_params_t, int n_dims>
    class ContrailCLEAN;

    template<class case_ct_params_t>
    class ContrailCLEAN<case_ct_params_t, 2> : public ContrailCLEANCommon<case_ct_params_t, 2>
    {
      using parent_t = ContrailCLEANCommon<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // function expecting a libmpdata concurr parameters struct as argument
      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      // function expecting a libmpdata++ concurr as argument
      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed, const int nps[2]) override
      {
        blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, th_e, rv_e, rl_e, rng_seed, k);

        // int nx = nps[0];
        // real_t dx = (this->X / si::metres) / (nx-1); 
        // int nz = nps[1];
        // real_t dz = (this->Z / si::metres) / (nz-1); 
      }

      virtual void update_rv_LS(blitz::Array<real_t, 2> rv_LS, const real_t t, const real_t dx, const real_t dy, const real_t dz)
      {
        if(t==0) rv_LS=0;
        else if (t > this->spinup_time)
        {
          blitz::firstIndex i;
          blitz::secondIndex k;

          // define a large-scale source of water vapor to represent the exhaust
          // dx is in meters, so i*dx is x position of given cell
          // t is in seconds
          // rv_LS is in [1/s] 
          rv_LS = where(
            (i * dx > 100 && i*dx < 400 && k * dz > 20 && k*dz < 80), // condition
            1e-4, // if true
            0     // if false
          );
        }
      };

      // same as above, but for potential temperature source [K/s]
      // virtual void update_th_LS(blitz::Array<real_t, 2> th_LS, const real_t t, const real_t dx, const real_t dy, const real_t dz)
      // {}

      public:
      ContrailCLEAN(const real_t _X=-1, const real_t _Y=-1, const real_t _Z=-1):
        parent_t(_X, _Y, _Z)
      {}
    };

    template<class case_ct_params_t>
    class ContrailCLEAN<case_ct_params_t, 3> : public ContrailCLEANCommon<case_ct_params_t, 3>
    {
      using parent_t = ContrailCLEANCommon<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // function expecting a libmpdata concurr parameters struct as argument
      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Y / si::metres) / (nps[1]-1);
        params.dk = (this->Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      // function expecting a libmpdata++ concurr as argument
      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed, const int nps[3]) override
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(concurr, rhod, th_e, rv_e, rl_e, rng_seed, k);

        concurr.advectee(ix::v) = 0;
        concurr.vab_relaxed_state(1) = 0;
      }

      public:
      ContrailCLEAN(const real_t _X=-1, const real_t _Y=-1, const real_t _Z=-1):
        parent_t(_X, _Y, _Z)
      {}
    };
  };
};

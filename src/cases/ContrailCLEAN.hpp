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
    quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
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

    const real_t z_abs = 100000; // no absorber

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
        params.th_src = false;
        params.rv_src = false;
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

        this->setopts_sgs(params);
      }
    
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr,
                        arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, int rng_seed, index_t index)
      {
        int nz = concurr.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (this->Z / si::metres) / (nz-1); 
    
        // concurr.advectee(ix::u) = 0;
        concurr.advectee(ix::w) = 0;  
       
        // absorbers; TODO: they are needed, because gravity waves are reflected from boundaries?
        concurr.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (this->Z / si::metres - z_abs)), 2), 0);
        concurr.vab_relaxed_state(0) = 0;
        concurr.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
    
        // density profile
        concurr.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
    
        concurr.advectee(ix::rv) = rv_e(index); 
        // initial potential temperature
        concurr.advectee(ix::th) = th_e(index); 

        std::cerr << "th_e:" << th_e << std::endl;
        std::cerr << "rv_e:" << rv_e << std::endl;
        std::cerr << "rhod:" << rhod << std::endl;
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
        // large-scale horizontal advection
        profs.th_LS = 0;
        profs.rv_LS = 0;
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
        // this->mean_rd1 = real_t(.03e-6) * si::metres,
        // this->mean_rd2 = real_t(.14e-6) * si::metres;
        // this->sdev_rd1 = real_t(1.28),
        // this->sdev_rd2 = real_t(1.75);
        // this->n1_stp = real_t(90e6) / si::cubic_metres, // 125 || 31
        // this->n2_stp = real_t(15e6) / si::cubic_metres;  // 65 || 16
        // this->ForceParameters.coriolis_parameter = 0.376e-4; // [1/s] 
        // this->z_rlx = z_rlx;
        // this->gccn_max_height = gccn_max_height;
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
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, th_e, rv_e, rl_e, rng_seed, k);

        int nx = concurr.advectee_global().extent(0);
        real_t dx = (this->X / si::metres) / (nx-1); 
        int nz = concurr.advectee_global().extent(ix::w);
        real_t dz = (this->Z / si::metres) / (nz-1); 

        blitz::firstIndex i;

        int mid_k = nz / 2;

        concurr.advectee(ix::u) = 0 
          + where(
            // if
            // i == 0 && k > mid_k && k <= mid_k + 1, 
            // (i == 0 || i == 1) && k == mid_k, 
            i == 0 && k == mid_k, 
            // k == mid_k, 
            // then
            1,
            // 67.5, 
            // else
            0
          );
        // concurr.advectee(ix::w) = 0 
        //   + where(
        //     // if
        //     // i == 0 && k > mid_k && k <= mid_k + 1, 
        //     // (i == 0 || i == 1) && k == mid_k, 
        //     i == 0 && k == mid_k, 
        //     // k == mid_k, 
        //     // then
        //     -67.5, 
        //     // else
        //     0
        //   );

      }

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
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
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

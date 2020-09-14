#pragma once

// like thermal test from limpdataxx
#include "CasesCommon.hpp"

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

namespace setup 
{
  namespace dry_thermal
  {
    const quantity<si::length, real_t> 
      Z    = 2000 * si::metres,
      Y    = 2000 * si::metres,
      X    = 2000 * si::metres;
  
    const real_t z_abs = 100000; // no absorber

    template<class case_ct_params_t, int n_dims>
    class DryThermalCommon : public CasesCommon<case_ct_params_t, n_dims>
    {
      protected:
      
      using parent_t = CasesCommon<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts_hlpr(rt_params_t &params, const user_params_t &user_params)
      {
        params.outdir = user_params.outdir;
        params.outfreq = user_params.outfreq;
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
        params.buoyancy_wet = false;
        params.subsidence = false;
        params.friction = false;
        params.coriolis = false;
        params.radiation = false;

        this->setopts_sgs(params);
      }
  
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &solver, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        int nz = solver.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (Z / si::metres) / (nz-1); 
        int nx = solver.advectee_global().extent(0);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dx = (X / si::metres) / (nx-1); 
    
    //    solver.advectee(ix::rv) = r_t()(index * dz); 
        solver.advectee(ix::u)= 0;// setup::u()(index * dz);
        solver.advectee(ix::w) = 0;  
       
        // absorbers
        solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);
        solver.vab_relaxed_state(0) = solver.advectee(ix::u);
        solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
    
        // density profile
        solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
    
        // initial potential temperature
        real_t r0 = 250;
        solver.advectee(ix::rv) = 1e-3; // some rv, but no actual wet physics
        solver.advectee(ix::th) = 300. + where(
          // if
          pow(blitz::tensor::i * dx - 4    * r0 , 2) + 
          pow(blitz::tensor::j * dz - 1.04 * r0 , 2) <= pow(r0, 2), 
          // then
          .5, 
          // else
          0
        );
      }
  
      // calculate the initial environmental theta and rv profiles
      // alse set w_LS and hgt_fctrs
      // like in Wojtek's BabyEulag
      void set_profs(profiles_t &profs, int nz, const user_params_t &user_params)
      {
        using libcloudphxx::common::theta_std::p_1000;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d_over_c_pd;

        parent_t::set_profs(profs, nz, user_params);
       
        profs.rhod = 1;
        profs.th_e = 300;
        profs.rv_e = 0; // doesnt matter in dry case, just to have consistent output between runs

        const quantity<si::temperature, real_t> T(
          libcloudphxx::common::theta_dry::T(
            quantity<si::temperature, real_t>(300 * si::kelvins),
            quantity<si::mass_density, real_t>(1 * si::kilograms / si::cubic_metres)
          )
        );

        const quantity<si::pressure, real_t> p(
          libcloudphxx::common::theta_dry::p(
            quantity<si::mass_density, real_t>(1 * si::kilograms / si::cubic_metres),
            quantity<si::dimensionless, real_t>(0.),
            T
          )
        );
        profs.p_e = real_t(p / si::pascals); // total env pressure

        profs.th_ref = 300;
        profs.w_LS = 0.;  // no subsidence
        profs.th_LS = 0.; // no large-scale horizontal advection
        profs.rv_LS = 0.; 
      }

      public:
      DryThermalCommon()
      {
        this->Z = Z;
      }
    };
    
    // 2d/3d children
    template<class case_ct_params_t, int n_dims>
    class DryThermal;

    template<class case_ct_params_t>
    class DryThermal<case_ct_params_t, 2> : public DryThermalCommon<case_ct_params_t, 2>
    {
      using parent_t = DryThermalCommon<case_ct_params_t, 2>;
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
        this->intcond_hlpr(solver, rhod, rng_seed, k);
      }

      public:
      DryThermal()
      {
        this->X = X;
      }
    };

    template<class case_ct_params_t>
    class DryThermal<case_ct_params_t, 3> : public DryThermalCommon<case_ct_params_t, 3>
    {
      using parent_t = DryThermalCommon<case_ct_params_t, 3>;
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
        this->intcond_hlpr(solver, rhod, rng_seed, k);
    
        solver.advectee(ix::v) = 0;
        solver.vab_relaxed_state(1) = 0;
      }

      public:
      // TODO: make it work in 3d?
      DryThermal()
      {
        throw std::runtime_error("Dry Thermal doesn't work in 3d yet");
      }
    };
  };
};


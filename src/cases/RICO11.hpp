// RICO trade cumulus based on 
// http://projects.knmi.nl/rico/setup3d.html
// https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011MS000056

#pragma once
#include <random>
#include "Anelastic.hpp"
#include "detail/formulas.hpp"

namespace setup 
{
  namespace rico
  {
    namespace hydrostatic = libcloudphxx::common::hydrostatic;
    namespace theta_std   = libcloudphxx::common::theta_std;
    namespace theta_dry   = libcloudphxx::common::theta_dry;
    namespace lognormal   = libcloudphxx::common::lognormal;
    namespace const_cp    = libcloudphxx::common::const_cp;

  
    const quantity<si::pressure, real_t> 
      p_0 = 101540 * si::pascals;
    const quantity<si::temperature, real_t> 
      T_SST = real_t(299.8) * si::kelvins;
    const quantity<si::length, real_t> 
      z_0  = 0    * si::metres,
      Z    = 4000 * si::metres, 
      X    = 12800 * si::metres, 
      Y    = 12800 * si::metres; 
    const real_t z_abs = 3000;
//    const real_t z_i = 795; //initial inversion height
    const quantity<si::length, real_t> z_rlx = 100 * si::metres;

    inline quantity<si::temperature, real_t> th_l_rico(const real_t &z)
    {
      quantity<si::temperature, real_t> ret;
      ret = z < 740. ?
        297.9 * si::kelvins : 
        (297.9 + (317. - 297.9)/(4000. - 740) * (z-740)) * si::kelvins;
      return ret;
    }

    inline quantity<si::dimensionless, real_t> r_t_rico(const real_t &z)
    {
      const quantity<si::dimensionless, real_t> q_t = z < 740 ?
        (16 + (13.8 - 16) / 740. * z) * 1e-3 : 
        z < 3260 ?
          (13.8 + (2.4 - 13.8) / (3260 - 740) * (z-740)) * 1e-3 :
          (2.4 + (1.8 - 2.4) / (4000 - 3260) * (z-3260)) * 1e-3;
      return q_t;
    }

    template<class case_ct_params_t, int n_dims>
    class Rico11Common : public Anelastic<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return th_l_rico(z);
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return r_t_rico(z);
      }

      // water mixing ratio at height z
      struct r_t_fctr
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return r_t_rico(z);
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };

      // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return th_l_rico(z) / si::kelvins;
        }
        BZ_DECLARE_FUNCTOR(th_std_fctr);
      };
    
      // westerly wind
      struct u
      {
        real_t operator()(const real_t &z) const
        {
          return -9.9 + 2e-3 * z; 
        }
        BZ_DECLARE_FUNCTOR(u);
      };
    
      // large-scale vertical wind
      struct w_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          real_t sub_vel = z < 2260 ?
            -(0.005 / 2260) * z :
            -0.005;
          return sub_vel; 
        }
        BZ_DECLARE_FUNCTOR(w_LS_fctr);
      };
    
      // large-scale horizontal advection of th + radiative cooling [K/s]
      struct th_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return -2.5 / 86400;
        }
        BZ_DECLARE_FUNCTOR(th_LS_fctr);
      };
    
      // large-scale horizontal advection of rv [1/s]
      struct rv_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          real_t rv_LS = z < 2980 ?
            -1. / 86400 + (1.3456 / 86400) * z / 2980 :
            4e-6;
          return rv_LS * 1e-3; 
        }
        BZ_DECLARE_FUNCTOR(rv_LS_fctr);
      };
    
      // density profile as a function of altitude
      // hydrostatic and assuming constant theta (not used now)
  /*    struct rhod_fctr
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
  */

      template<bool enable_sgs = case_ct_params_t::enable_sgs>
      void setopts_sgs(rt_params_t &params,
                       typename std::enable_if<!enable_sgs>::type* = 0) 
      {
        parent_t::setopts_sgs(params);
      }

      template<bool enable_sgs = case_ct_params_t::enable_sgs>
      void setopts_sgs(rt_params_t &params,
                       typename std::enable_if<enable_sgs>::type* = 0) 
      {
        parent_t::setopts_sgs(params);
        params.cdrag = 0.001229; // NOTE: in SMG simulations, we do not apply the height correction to drag coefficient (i.e. it is assumed that U ground is at 20 m)
      }
  
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
        params.rc_src = user_params.rc_src;
        params.rr_src = user_params.rr_src;
        params.nc_src = user_params.nc_src;
        params.nr_src = user_params.nr_src;
        params.dt = user_params.dt;
        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = true;
        params.vel_subsidence = false;
        params.friction = true;
        params.coriolis = true;
        params.radiation = false;

        this->setopts_sgs(params);
      }
  
  
  
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &solver, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        int nz = solver.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (Z / si::metres) / (nz-1); 
  
        solver.advectee(ix::rv) = r_t_fctr{}(index * dz); 
        solver.advectee(ix::u)= u{}(index * dz);
        solver.advectee(ix::w) = 0;  
       
        // absorbers
        solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 1020 * (index * dz - z_abs) / (Z / si::metres - z_abs), 0);
        solver.vab_relaxed_state(0) = solver.advectee(ix::u);
        solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
  
        // density profile
        solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
  
        // initial potential temperature
        solver.advectee(ix::th) = th_std_fctr{}(index * dz); 

        // randomly prtrb tht and rv
        // NOTE: all processes do this, but ultimately only perturbation calculated by MPI rank 0 is used
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-0.1, 0.1);
          auto rand = std::bind(dis, gen);
  
          auto th_global = solver.advectee_global(ix::th);
          decltype(solver.advectee(ix::th)) prtrb(th_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          th_global += prtrb;
          this->make_cyclic(th_global);
          solver.advectee_global_set(th_global, ix::th);
        }
        {
          std::mt19937 gen(rng_seed+1); // different seed than in th. NOTE: if the same instance of gen is used in th and rv, for some reason it gives the same sequence in rv as in th despite being advanced in th prtrb
          std::uniform_real_distribution<> dis(-0.025e-3, 0.025e-3);
          auto rand = std::bind(dis, gen);
  
          auto rv_global = solver.advectee_global(ix::rv);
          decltype(solver.advectee(ix::rv)) prtrb(rv_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          rv_global += prtrb;
          this->make_cyclic(rv_global);
          solver.advectee_global_set(rv_global, ix::rv);
        }
      }
  
  
  
      // calculate the initial environmental theta and rv profiles
      // like in Wojtek's BabyEulag
      // alse set w_LS and hgt_fctrs
      // TODO: same in DYCOMS (and others?), move to a common function
      void set_profs(profiles_t &profs, int nz, const user_params_t &user_params)
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
        profs.w_LS = w_LS_fctr()(k * dz);
        // large-scale horizontal advection
        profs.th_LS = th_LS_fctr()(k * dz);
        profs.rv_LS = rv_LS_fctr()(k * dz);
      }

      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                 blitz::Array<real_t, n_dims> th_ground,    // value of th on the ground
                                 blitz::Array<real_t, n_dims> U_ground,     // magnitude of horizontal ground wind
                                 const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        static const real_t th_0 = (T_SST / si::kelvins) / theta_std::exner(p_0);
        surf_flux_sens =   formulas::surf_flux_coeff_scaling<real_t>(U_ground_z, 20) * real_t(0.001094) * U_ground * (th_ground - th_0) * (this->rhod_0 / si::kilograms * si::cubic_meters) * theta_std::exner(p_0); // [K kg / (m^2 s)]; *= -1 because gradient is taken later and negative gradient of upward flux means inflow
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        static const real_t rsat_0 = const_cp::r_vs(T_SST, p_0); // if we wanted to use the Tetens formula, this would need to be changed
        surf_flux_lat =   formulas::surf_flux_coeff_scaling<real_t>(U_ground_z, 20) * real_t(0.001133) * U_ground * (rt_ground - rsat_0) * (this->rhod_0 / si::kilograms * si::cubic_meters); // [kg / (m^2 s)]
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind
                               const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        surf_flux_uv = - formulas::surf_flux_coeff_scaling<real_t>(U_ground_z, 20) * real_t(0.001229) * U_ground * uv_ground * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters); // [ kg m/s / (m^2 s) ]
      }

      // ctor
      Rico11Common()
      {
        this->p_0 = p_0;
        this->mean_rd1 = real_t(.03e-6) * si::metres,
        this->mean_rd2 = real_t(.14e-6) * si::metres;
        this->sdev_rd1 = real_t(1.28),
        this->sdev_rd2 = real_t(1.75);
        this->n1_stp = real_t(90e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(15e6) / si::cubic_metres;  // 65 || 16
        this->ForceParameters.coriolis_parameter = 0.449e-4; // [1/s] @ 18.0 deg N
        this->X = X;
        this->Z = Z;
        this->z_rlx = z_rlx;
      }
    };
    
    template<class case_ct_params_t, int n_dims>
    class Rico11;

    template<class case_ct_params_t>
    class Rico11<case_ct_params_t, 2> : public Rico11Common<case_ct_params_t, 2>
    {
      using parent_t = Rico11Common<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      void intcond(typename parent_t::concurr_any_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);

        auto th_global = solver.advectee_global(ix::th);
        this->make_cyclic(th_global);
        solver.advectee_global_set(th_global, ix::th);

        auto rv_global = solver.advectee_global(ix::rv);
        this->make_cyclic(rv_global);
        solver.advectee_global_set(rv_global, ix::rv);
      }

      public:
      Rico11()
      {
        this->X = X;
      }
    };

    template<class case_ct_params_t>
    class Rico11<case_ct_params_t, 3> : public Rico11Common<case_ct_params_t, 3>
    {
      using parent_t = Rico11Common<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // southerly wind
      struct v
      {
        real_t operator()(const real_t &z) const
        {
          return -3.8;
        }
        BZ_DECLARE_FUNCTOR(v);
      };

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Y / si::metres) / (nps[1]-1);
        params.dk = (Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      void intcond(typename parent_t::concurr_any_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);

        auto th_global = solver.advectee_global(ix::th);
        this->make_cyclic(th_global);
        solver.advectee_global_set(th_global, ix::th);

        auto rv_global = solver.advectee_global(ix::rv);
        this->make_cyclic(rv_global);
        solver.advectee_global_set(rv_global, ix::rv);
  
        int nz = solver.advectee_global().extent(ix::w);
        real_t dz = (Z / si::metres) / (nz-1); 
  
        solver.advectee(ix::v)= v()(k * dz);
        solver.vab_relaxed_state(1) = solver.advectee(ix::v);
      }

      void set_profs(profiles_t &profs, int nz, const user_params_t &user_params)
      {
        parent_t::set_profs(profs, nz, user_params);
        // geostrophic wind equal to the initial velocity profile
        blitz::firstIndex k;
        typename parent_t::u u;
        real_t dz = (Z / si::metres) / (nz-1);
        profs.geostr[0] = u(k * dz);
        profs.geostr[1] = v()(k * dz);
      }

      public:
      Rico11()
      {
        this->X = X;
        this->Y = Y;
      }
    };
  };
};

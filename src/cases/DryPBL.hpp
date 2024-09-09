// Dry planetary boundary layer simulation
// based on the pbl test from libmpdata++
// but here it's anelastic, not Boussinesq

#pragma once
#include <random>
#include "Anelastic.hpp"
#include "detail/formulas.hpp"

namespace cases 
{
  namespace pbl
  {
    namespace hydrostatic = libcloudphxx::common::hydrostatic;
    namespace theta_std   = libcloudphxx::common::theta_std;
    namespace theta_dry   = libcloudphxx::common::theta_dry;
    namespace lognormal   = libcloudphxx::common::lognormal;
    namespace const_cp    = libcloudphxx::common::const_cp;

  
    const quantity<si::pressure, real_t> 
      p_0 = 101540 * si::pascals;
    const quantity<si::length, real_t> 
      Z_def    = 1500 * si::metres, 
      X_def    = 3200 * si::metres, 
      Y_def    = 3200 * si::metres; 
    const real_t z_abs = 1000; //[m]
//    const real_t z_i = 795; //initial inversion height
    const quantity<si::length, real_t> z_rlx = 25 * si::metres;
    const real_t mixed_length = 500; // [m]
    const quantity<si::dimensionless, real_t> st = 1e-5;

    inline quantity<si::temperature, real_t> th_l_pbl(const real_t &z)
    {
      quantity<si::temperature, real_t> ret = real_t(300) * std::max(real_t(1), real_t(1 + (z - mixed_length) * st)) * si::kelvins; // tht_e from pbl in libmpdata++, tht there is 0

      return ret;
    }

    inline quantity<si::dimensionless, real_t> r_t_pbl(const real_t &z)
    {
      const quantity<si::dimensionless, real_t> q_t = 0;
      return q_t;
    }

    template<class case_ct_params_t, int n_dims>
    class DryPBLCommon : public Anelastic<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return th_l_pbl(z);
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return r_t_pbl(z);
      }

      // water mixing ratio at height z
      struct r_t_fctr
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return r_t_pbl(z);
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };

      // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return th_l_pbl(z) / si::kelvins;
        }
        BZ_DECLARE_FUNCTOR(th_std_fctr);
      };
    
      // large-scale vertical wind
      struct w_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          real_t sub_vel = 0;
          return sub_vel; 
        }
        BZ_DECLARE_FUNCTOR(w_LS_fctr);
      };
    
      // large-scale horizontal advection of th + radiative cooling [K/s]
      struct th_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return 0;
        }
        BZ_DECLARE_FUNCTOR(th_LS_fctr);
      };
    
      // large-scale horizontal advection of rv [1/s]
      struct rv_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          real_t rv_LS = 0;
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
        params.cdrag = 0.1; // NOTE: in SMG simulations, we do not apply the height correction to drag coefficient (i.e. it is assumed that U ground is at 20 m)
      }
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
        params.buoyancy_wet = false;
        params.subsidence = subs_t::none;
        params.vel_subsidence = false;
        params.friction = true;
        params.coriolis = false;
        params.radiation = false;

        params.rv_src = false;
        params.rc_src = false;
        params.rr_src = false;
        params.nc_src = false;
        params.nr_src = false;

        params.user_params.relax_th_rv = true;

        this->setopts_sgs(params);
      }
  
  
  
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        int nz = concurr.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (this->Z / si::metres) / (nz-1); 
  
        concurr.advectee(ix::rv) = r_t_fctr{}(index * dz); 
        concurr.advectee(ix::u) = 0;
        concurr.advectee(ix::w) = 0;  
       
        // absorbers
        concurr.vab_coefficient() = where(index * dz >= z_abs,  1. / 1020 * (index * dz - z_abs) / (this->Z / si::metres - z_abs), 0);
        concurr.vab_relaxed_state(0) = 0;
        concurr.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
  
        // density profile
        concurr.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
  
        // initial potential temperature
        concurr.advectee(ix::th) = th_std_fctr{}(index * dz); 

        // randomly prtrb tht and w
        // NOTE: all processes do this, but ultimately only perturbation calculated by MPI rank 0 is used
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-0.0005, 0.0005);
          auto rand = std::bind(dis, gen);
  
          auto th_global = concurr.advectee_global(ix::th);
          decltype(concurr.advectee(ix::th)) prtrb(th_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          prtrb = where(index * dz >= mixed_length, 0, prtrb * (1. - (index * dz / mixed_length)));
          th_global += prtrb;
          this->make_cyclic(th_global);
          concurr.advectee_global_set(th_global, ix::th);
        }
        {
          std::mt19937 gen(rng_seed); // same seed as in th on purpose
          std::uniform_real_distribution<> dis(-0.1, 0.1);
          auto rand = std::bind(dis, gen);
  
          auto w_global = concurr.advectee_global(ix::w);
          decltype(concurr.advectee(ix::w)) prtrb(w_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          prtrb = where(index * dz >= mixed_length, 0, prtrb * (1. - (index * dz / mixed_length)));
          w_global += prtrb;
          this->make_cyclic(w_global);
          concurr.advectee_global_set(w_global, ix::w);
        }
      }
  
      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
      {
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

        //nudging of th and rv
        profs.relax_th_rv_coeff = where(k * dz >= z_abs, 1. / 1020 * (k * dz - z_abs) / (this->Z / si::metres - z_abs), 0);
      }

      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                 blitz::Array<real_t, n_dims> th_ground,    // value of th on the ground
                                 blitz::Array<real_t, n_dims> U_ground,     // magnitude of horizontal ground wind
                                 const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        surf_flux_sens = .01 * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters) * theta_std::exner(p_0); // [K kg / (m^2 s)]; -1 because negative gradient of upward flux means inflow

      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        surf_flux_lat = 0;
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind
                               const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy, const real_t &uv_mean) override
      {
        surf_flux_uv = - real_t(0.1) * U_ground * (uv_ground + uv_mean) * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters); // [ kg m/s / (m^2 s) ]
      }

      // ctor
      DryPBLCommon()
      {
        this->p_0 = p_0;
        this->z_rlx = z_rlx;
        // TODO: no aerosols!
      }
    };
    
    template<class case_ct_params_t, int n_dims>
    class DryPBL;

    template<class case_ct_params_t>
    class DryPBL<case_ct_params_t, 2> : public DryPBLCommon<case_ct_params_t, 2>
    {
      using parent_t = DryPBLCommon<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      void intcond(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);

        auto th_global = concurr.advectee_global(ix::th);
        this->make_cyclic(th_global);
        concurr.advectee_global_set(th_global, ix::th);

        auto rv_global = concurr.advectee_global(ix::rv);
        this->make_cyclic(rv_global);
        concurr.advectee_global_set(rv_global, ix::rv);
      }

      public:
      DryPBL(const real_t _X=-1, const real_t _Y=-1, const real_t _Z=-1)
      {
        this->X = _X < 0 ? X_def : _X * si::meters;
        this->Z = _Z < 0 ? Z_def : _Z * si::meters;
      }
    };

    template<class case_ct_params_t>
    class DryPBL<case_ct_params_t, 3> : public DryPBLCommon<case_ct_params_t, 3>
    {
      using parent_t = DryPBLCommon<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Y / si::metres) / (nps[1]-1);
        params.dk = (this->Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      void intcond(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);

        auto th_global = concurr.advectee_global(ix::th);
        this->make_cyclic(th_global);
        concurr.advectee_global_set(th_global, ix::th);

        auto rv_global = concurr.advectee_global(ix::rv);
        this->make_cyclic(rv_global);
        concurr.advectee_global_set(rv_global, ix::rv);
  
        int nz = concurr.advectee_global().extent(ix::w);
        real_t dz = (this->Z / si::metres) / (nz-1); 
  
        concurr.advectee(ix::v)= 0;
        concurr.vab_relaxed_state(1) = 0;
      }

      public:
      DryPBL(const real_t _X=-1, const real_t _Y=-1, const real_t _Z=-1)
      {
        this->X = _X < 0 ? X_def : _X * si::meters;
        this->Y = _Y < 0 ? Y_def : _Y * si::meters;
        this->Z = _Z < 0 ? Z_def : _Z * si::meters;
      }
    };
  };
};

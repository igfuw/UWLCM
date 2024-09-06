// BOMEX trade cumulus based on Siebesma et al. 2003
//
// TODO:  subsidence should be proportional to the gradient of HORIZONTALLLY AVERAGED fields, but
//        it is the local gradient (same goes for SD subsidence) 

#pragma once
#include <random>
#include "Anelastic.hpp"
#include "detail/formulas.hpp"

namespace cases 
{
  namespace bomex
  {
    namespace hydrostatic = libcloudphxx::common::hydrostatic;
    namespace theta_std   = libcloudphxx::common::theta_std;
    namespace theta_dry   = libcloudphxx::common::theta_dry;
    namespace lognormal   = libcloudphxx::common::lognormal;
    namespace const_cp    = libcloudphxx::common::const_cp;

  
    const quantity<si::pressure, real_t> 
      p_0 = 101500 * si::pascals;
    const quantity<si::length, real_t> 
      z_0      = 0      * si::metres,
      Z_def    = 75*40  * si::metres, 
      X_def    = 64*100 * si::metres, 
      Y_def    = 64*100 * si::metres; 
    const real_t z_abs = 2500; // no word about absorber in Siebesma (2003)
    const quantity<si::length, real_t> z_rlx = 100 * si::metres;
    const quantity<si::length, real_t> gccn_max_height = 450 * si::metres; // below cloud base

    inline real_t linint(const real_t a, const real_t b, const real_t fa, const real_t fb, const real_t x)
    {
      return((x - a) / (b - a) * (fb - fa) + fa); 
    }

    inline real_t init_prof(const std::vector<real_t> &lvls, const std::vector<real_t> &vals, const real_t z)
    {
      for(int i=0; i<lvls.size(); ++i)
      {
        if(z < lvls.at(i)) return linint(lvls.at(i-1), lvls.at(i), vals.at(i-1), vals.at(i), z); 
      }
      return vals.back();
    }

    inline quantity<si::velocity, real_t> u_bomex(const real_t &z)
    {
      return init_prof({0, 700, 3000}, {-8.75, -8.75, -4.61}, z) * si::meters / si::seconds;
    }
    inline quantity<si::velocity, real_t> v_bomex(const real_t &z)
    {
      return real_t(0) * si::meters / si::seconds;
    }

    inline quantity<si::temperature, real_t> th_l_bomex(const real_t &z)
    {
      return init_prof({0, 520, 1480, 2000, 3000}, {298.7, 298.7, 302.4, 308.2, 311.85}, z) * si::kelvins;
    }

    inline quantity<si::dimensionless, real_t> r_t_bomex(const real_t &z)
    {
      // q_t is specific humidity
      //const quantity<si::dimensionless, real_t> q_t =
      real_t q_t =
        init_prof({0, 520, 1480, 2000, 3000}, {17., 16.3, 10.7, 4.2, 3.0}, z) * 1e-3;
      // turn to mixing ratio
      q_t = q_t / (1. - q_t);
      return q_t;
    }

    template<class case_ct_params_t, int n_dims>
    class Bomex03Common : public Anelastic<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return th_l_bomex(z);
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return r_t_bomex(z);
      }

      // water mixing ratio at height z
      struct r_t_fctr
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return r_t_bomex(z);
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };

      // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return th_l_bomex(z) / si::kelvins;
        }
        BZ_DECLARE_FUNCTOR(th_std_fctr);
      };
    
      // westerly wind
      struct u_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        u_t() : hori_vel_t(&u_bomex) {}

        BZ_DECLARE_FUNCTOR(u_t);
      };

      u_t u;
    
      // large-scale vertical wind
      struct w_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return init_prof({0, 1500, 2100, 3000}, {0, -0.65e-2, 0, 0}, z); 
        }
        BZ_DECLARE_FUNCTOR(w_LS_fctr);
      };
    
      // large-scale horizontal advection of th + radiative cooling [K/s]
      struct th_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return init_prof({0, 1500, 3000}, {-2, -2, 0}, z) / (24. * 3600.); // [K/s]
        }
        BZ_DECLARE_FUNCTOR(th_LS_fctr);
      };
    
      // large-scale horizontal advection of rv [1/s]
      struct rv_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          // assume dr/dt = dq/dt
          // in reality, dr/dq = (1-q)^{-1} + q * (1-q)^{-2}
          return init_prof({0, 300, 500, 3000}, {-1.2e-8, -1.2e-8, 0, 0}, z);
        }
        BZ_DECLARE_FUNCTOR(rv_LS_fctr);
      };

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
        params.fricvelsq = 0.0784; // 0.28^2
      }
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
        params.buoyancy_wet = true;
        params.subsidence = true;
        params.vel_subsidence = true;
        params.friction = true;
        params.coriolis = true;
        params.radiation = false; // that's Dycoms-style radiation

        this->setopts_sgs(params);
      }
  
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        int nz = concurr.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (this->Z / si::metres) / (nz-1); 
  
        concurr.advectee(ix::rv) = r_t_fctr{}(index * dz); 
        concurr.advectee(ix::u)= u(index * dz);
        concurr.advectee(ix::w) = 0;  
       
        // absorbers
        concurr.vab_coefficient() = where(index * dz >= z_abs,  1. / 1020 * (index * dz - z_abs) / (this->Z / si::metres - z_abs), 0);
        concurr.vab_relaxed_state(0) = concurr.advectee(ix::u);
        concurr.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
  
        // density profile
        concurr.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
  
        // initial potential temperature
        concurr.advectee(ix::th) = th_std_fctr{}(index * dz); 

        // randomly prtrb tht and rv
        // NOTE: all processes do this, but ultimately only perturbation calculated by MPI rank 0 is used
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-0.1, 0.1);
          auto rand = std::bind(dis, gen);
  
          auto th_global = concurr.advectee_global(ix::th);
          decltype(concurr.advectee(ix::th)) prtrb(th_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          prtrb = where(index * dz >= 1600., 0., prtrb); // no perturbation above 1.6km
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
          prtrb = where(index * dz >= 1600., 0., prtrb); // no perturbation above 1.6km
          rv_global += prtrb;
          this->make_cyclic(rv_global);
          concurr.advectee_global_set(rv_global, ix::rv);
        }
      }
  
  
  
      // calculate the initial environmental theta and rv profiles
      // like in Wojtek's BabyEulag
      // alse set w_LS and hgt_fctrs
      // TODO: same in DYCOMS (and others?), move to a common function
      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
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
        if(timestep == 0)
          surf_flux_sens = 8e-3 * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters) * theta_std::exner(p_0); // [K kg / (m^2 s)]; -1 because negative gradient of upward flux means inflow
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        if(timestep == 0)
          surf_flux_lat = 5.2e-5 * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters); // [m/s]
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground (without mean)
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind (total, including mean)
                               const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy, const real_t &uv_mean) override
      {
        surf_flux_uv = where(U_ground < 1e-4,
            - 0.0784 * (uv_ground + uv_mean) / real_t(1e-4) * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters), // 0.0784 m^2 / s^2 is the square of friction velocity = 0.28 m / s
            - 0.0784 * (uv_ground + uv_mean) / U_ground * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters)
          );
      }

/*
      void update_rv_LS(blitz::Array<real_t, 1> rv_LS,
                        int timestep, const real_t dt, real_t dz)
      {
        blitz::firstIndex k;
        rv_LS = rv_LS_var_fctr(timestep * dt)(k * dz);
      };

      void update_th_LS(blitz::Array<real_t, 1> th_LS,
                        int timestep, const real_t dt, real_t dz)
      {
        blitz::firstIndex k;
        th_LS = th_LS_var_fctr(timestep * dt)(k * dz);
      };
      */

      void init()
      {
        // Siebesma 2003 does not specify aerosol
        // use one from RICO (vanZanten 2011) as in Sato, Shima et al. 2017
        this->p_0 = p_0;
        this->mean_rd1 = real_t(.03e-6) * si::metres,
        this->mean_rd2 = real_t(.14e-6) * si::metres;
        this->sdev_rd1 = real_t(1.28),
        this->sdev_rd2 = real_t(1.75);
        this->n1_stp = real_t(90e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(15e6) / si::cubic_metres;  // 65 || 16
        this->ForceParameters.coriolis_parameter = 0.376e-4; // [1/s] 
        this->z_rlx = z_rlx;
        this->gccn_max_height = gccn_max_height;
      }

      public:

      // ctor
      Bomex03Common(const real_t _X, const real_t _Y, const real_t _Z, const bool window)
      {
        init();

        this->X = _X < 0 ? X_def : _X * si::meters;
        if(n_dims == 3)
          this->Y = _Y < 0 ? Y_def : _Y * si::meters;
        this->Z = _Z < 0 ? Z_def : _Z * si::meters;

        u.init(window, this->Z);

        this->ForceParameters.uv_mean[0] = u.mean_vel;
      }
    };
    
    template<class case_ct_params_t, int n_dims>
    class Bomex03;


    // TODO: stuff from 2D can be moved to common? 
    template<class case_ct_params_t>
    class Bomex03<case_ct_params_t, 2> : public Bomex03Common<case_ct_params_t, 2>
    {
      using parent_t = Bomex03Common<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;
      using parent_t::parent_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      void intcond(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
       // blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, blitz::secondIndex{});
      }
    };

    template<class case_ct_params_t>
    class Bomex03<case_ct_params_t, 3> : public Bomex03Common<case_ct_params_t, 3>
    {
      using parent_t = Bomex03Common<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // southerly wind
      struct v_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        v_t() : hori_vel_t(&v_bomex) {}

        BZ_DECLARE_FUNCTOR(v_t);
      };

      v_t v;

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
  
        int nz = concurr.advectee_global().extent(ix::w);
        real_t dz = (this->Z / si::metres) / (nz-1); 
  
        concurr.advectee(ix::v)= v(k * dz);
        concurr.vab_relaxed_state(1) = concurr.advectee(ix::v);
      }

      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
      {
        parent_t::set_profs(profs, nz, user_params);
        // geostrophic wind used in Coriolis force
        blitz::firstIndex k;
        real_t dz = (this->Z / si::metres) / (nz-1);
        profs.geostr[0] = -10 + 1.8e-3 * k * dz;
        profs.geostr[1] = 0;
      }

      public:
      Bomex03(const real_t _X, const real_t _Y, const real_t _Z, const bool window):
        parent_t(_X, _Y, _Z, window)
        {
          v.init(window, this->Z);
          this->ForceParameters.uv_mean[1] = v.mean_vel;
        }
    };
  };
};

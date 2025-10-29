// Daytime convective development over land based on Grabowski et al. (2006)
//  https://doi.org/10.1256/qj.04.147
// https://www2.mmm.ucar.edu/gcss-wg4/gcss/case4.html
#pragma once
#include <random>
#include "Anelastic.hpp"
#include "detail/formulas.hpp"
#include "detail/LBA_sounding/input_sounding.hpp"
#include "detail/LBA_sounding/radiative_cooling.hpp"

namespace cases
{
  namespace daytime_convection
  {
    namespace hydrostatic = libcloudphxx::common::hydrostatic;
    namespace theta_std   = libcloudphxx::common::theta_std;
    namespace theta_dry   = libcloudphxx::common::theta_dry;
    namespace lognormal   = libcloudphxx::common::lognormal;
    namespace const_cp    = libcloudphxx::common::const_cp;

     // RH T and p to rv assuming RH = r_v / r_vs
    inline quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return RH * const_cp::r_vs<real_t>(T, p);
    }
  
    const quantity<si::pressure, real_t> 
      p_0 = real_t(99130) * si::pascals;
    const quantity<si::temperature, real_t> 
      T_SST = real_t(296.85) * si::kelvins;
    const quantity<si::length, real_t> 
      Z_def    = 20000 * si::metres,
      X_def    = 64000 * si::metres,
      Y_def    = 64000 * si::metres;
    const real_t z_abs = 16000;
    const quantity<si::length, real_t> z_rlx = 100 * si::metres;


    // returned units: [Celsius], [%], [m/s], [hPa]
    inline real_t interpolate_LBA_sounding(const std::string valname, real_t pos)
    {
      assert(pos>=0);
      assert(valname == "T" || valname == "RH" || valname == "u" || valname == "v" || valname == "p");

      const auto &heights(LBA_sounding_z);

      const auto pos_up = std::upper_bound(heights.begin(), heights.end(), pos); // pos_up is the height index above pos
      																			 // pos_up-1 will be the index below

      if(pos_up == heights.end())
        throw std::runtime_error("UWLCM: The initial sounding is not high enough");
      if(pos_up == heights.begin())
        throw std::runtime_error("UWLCM: The initial profile has incorrect first element");

      const auto &sounding(
        valname == "T" ? LBA_sounding_T :
        valname == "RH" ? LBA_sounding_RH :
        valname == "u" ? LBA_sounding_u :
        valname == "v" ? LBA_sounding_v :
        LBA_sounding_p);

      const auto sounding_upper = sounding.begin() + std::distance(heights.begin(), pos_up);
      return real_t(*(sounding_upper-1) + (pos - *(pos_up-1)) / (*pos_up - *(pos_up-1)) * (*sounding_upper - *(sounding_upper-1)));
    }


    inline quantity<si::velocity, real_t> u_lba(const real_t &z)
    {
      return interpolate_LBA_sounding("u", z) * si::meters / si::seconds;
    }

    inline quantity<si::velocity, real_t> v_lba(const real_t &z)
    {
      return interpolate_LBA_sounding("v", z) * si::meters / si::seconds;
    }

    inline quantity<si::temperature, real_t> T_lba(const real_t &z)
    {
      return (interpolate_LBA_sounding("T", z) + real_t(273.15)) * si::kelvins; //converting celsius to kelvins
    }

    inline quantity<si::dimensionless, real_t> RH_lba(const real_t &z)
    {
      return interpolate_LBA_sounding("RH", z) / real_t(100); //converting % to dimensionless
    }

    inline quantity<si::pressure, real_t> p_lba(const real_t &z)
    {
      return interpolate_LBA_sounding("p", z) * real_t(100) * si::pascals; //converting hPa to Pa
    }


    // returned units: [K/day]
    inline real_t interpolate_rad_cooling(real_t pos, real_t t) {
      assert(pos>=0); //height in m
      assert(t>=0);  //time in s

      t = t / real_t(60); // converting to minutes to match the rad_cooling_time

      const auto &heights = rad_cooling_z;
      const auto &times = rad_cooling_time;
      const auto &sounding = rad_cooling;

      auto pos_up = std::upper_bound(heights.begin(), heights.end(), pos); // pos_up is the height index above pos
      auto pos_low = pos_up - 1; // pos_up-1 will be the index below

      auto t_up = std::upper_bound(times.begin(), times.end(), t); // t_up is the time index after t
      auto t_low = t_up - 1; // t_up-1 will be the index before

      if(t_up == times.end())
        throw std::runtime_error("UWLCM: Time too long for radiative cooling interpolation");

      if(t_up == times.begin()) //use radiative cooling at 10 min between model start and t=10 min
      {
        const auto &time_data = *sounding.begin();
        if(pos_up == heights.end()) //radiative cooling above the last level provided should be set to zero
        {
          return real_t(0);
        }
        if(pos_up == heights.begin())
        {
          return real_t(*time_data.begin()); //radiative cooling at the first level can be used as the rate below this level
        }
        const auto sounding_upper = time_data.begin() + std::distance(heights.begin(), pos_up);
        return real_t(*(sounding_upper-1) + (pos - *(pos_up-1)) / (*pos_up - *(pos_up-1)) * (*sounding_upper - *(sounding_upper-1)));
      }

      const auto next_time_data_iter = sounding.begin() + std::distance(times.begin(), t_up);
      const auto previous_time_data_iter = sounding.begin() + std::distance(times.begin(), t_low);
      const auto &next_time_data = *next_time_data_iter;
      const auto &previous_time_data = *previous_time_data_iter;

      if(pos_up == heights.end()) //radiative heating above the last level provided should be set to zero
        return real_t(0);
      if(pos_up == heights.begin()) //radiative heating at the first level can be used as the rate below this level
      {
        const auto sounding_next_time = real_t(*next_time_data.begin());
        const auto sounding_previous_time = real_t(*previous_time_data.begin());
        return real_t(sounding_previous_time + (t - *(t_low)) / (*t_up - *(t_low)) * (sounding_next_time - sounding_previous_time));
      }

      const auto pos_up_next_time = next_time_data.begin() + std::distance(heights.begin(), pos_up);
      const auto pos_up_prevoius_time = previous_time_data.begin() + std::distance(heights.begin(), pos_up);
      const auto interp_next_time = real_t(*(pos_up_next_time-1) + (pos - *(pos_low)) / (*pos_up - *(pos_low)) * (*pos_up_next_time - *(pos_up_next_time-1)));
      const auto interp_previous_time = real_t(*(pos_up_prevoius_time-1) + (pos - *(pos_low)) / (*pos_up - *(pos_low)) * (*pos_up_prevoius_time - *(pos_up_prevoius_time-1)));
      return real_t(interp_previous_time + (t - *(t_low)) / (*t_up - *(t_low)) * (interp_next_time - interp_previous_time));
    }


    template<class case_ct_params_t, int n_dims>
    class DaytimeConvection_common : public Anelastic<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return T_lba(z) / theta_std::exner(p_lba(z));
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return RH_T_p_to_rv(RH_lba(z), T_lba(z), p_lba(z));
      }

      // water mixing ratio at height z
      struct r_t_fctr
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return RH_T_p_to_rv(RH_lba(z), T_lba(z), p_lba(z));
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };

      // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return T_lba(z) / theta_std::exner(p_lba(z)) / si::kelvins;
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

        u_t() : hori_vel_t(&u_lba) {}

        BZ_DECLARE_FUNCTOR(u_t);
      };

      u_t u;


      // time-varying radiative cooling [K/s]
      struct th_LS_var_fctr
      {
        real_t t;

        th_LS_var_fctr(const real_t &t): t(t) {}

        real_t operator()(const real_t &z) const // height in [m], time in [s]
        {
          return interpolate_rad_cooling(z, t) / real_t(86400); //converting from K/day to K/s
        }
        BZ_DECLARE_FUNCTOR(th_LS_var_fctr);
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
        params.cdrag = 0.;
      }
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
        params.buoyancy_wet = true;
        params.subsidence = subs_t::none;
        params.vel_subsidence = false;
        params.friction = false;
        params.coriolis = false;
        params.radiation = false;

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
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it
          prtrb = where(index * dz >= 1000., 0., prtrb); // no perturbation above 1km
          th_global += prtrb;
          this->make_cyclic(th_global);
          concurr.advectee_global_set(th_global, ix::th);
        }
        {
          std::mt19937 gen(rng_seed+1); // different seed than in th. NOTE: if the same instance of gen is used in th and rv, for some reason it gives the same sequence in rv as in th despite being advanced in th prtrb
          std::uniform_real_distribution<> dis(-1e-4, 1e-4);
          auto rand = std::bind(dis, gen);
  
          auto rv_global = concurr.advectee_global(ix::rv);
          decltype(concurr.advectee(ix::rv)) prtrb(rv_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it
          prtrb = where(index * dz >= 1000., 0., prtrb); // no perturbation above 1km
          rv_global += prtrb;
          this->make_cyclic(rv_global);
          concurr.advectee_global_set(rv_global, ix::rv);
        }
      }
  
  
  
      // calculate the initial environmental theta and rv profiles
      // like in Wojtek's BabyEulag
      // alse set w_LS and hgt_fctrs
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
  
//        // subsidence rate
//        profs.w_LS = w_LS_fctr()(k * dz);
//        // large-scale horizontal advection
        profs.th_LS = th_LS_var_fctr(real_t(0))(k * dz);
//        profs.rv_LS = rv_LS_fctr()(k * dz);
      }

      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                 blitz::Array<real_t, n_dims> th_ground,    // value of th on the ground
                                 blitz::Array<real_t, n_dims> U_ground,     // magnitude of horizontal ground wind
                                 const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        using libcloudphxx::common::moist_air::c_pd;
        real_t flux = formulas::surf_flux_function(timestep * dt);
        surf_flux_sens = -  real_t(270) * std::pow(flux, real_t(1.5)) / (c_pd<real_t>() / si::joules * si::kilogram * si::kelvin ); // [K kg / (m^2 s)]; *= -1 because gradient is taken later and negative gradient of upward flux means inflow
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        using libcloudphxx::common::const_cp::l_tri;
        real_t flux = formulas::surf_flux_function(timestep * dt);
        surf_flux_lat =  - real_t(554) * std::pow(flux, real_t(1.3)) / (l_tri<real_t>() / si::joules * si::kilograms); // [kg / (m^2 s)]
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground (without mean)
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind (total, including mean)
                               const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy, const real_t &uv_mean) override
      {
        surf_flux_uv = real_t(0); // [ kg m/s / (m^2 s) ]
      }


      void update_th_LS(blitz::Array<real_t, 1> th_LS,
                        int timestep, const real_t dt, real_t dz)
      {
        blitz::firstIndex k;
        th_LS = th_LS_var_fctr(timestep * dt)(k * dz);
      };


      void init()
      {
        this->p_0 = p_0;
        this->mean_rd1 = real_t(.03e-6) * si::metres,
        this->mean_rd2 = real_t(.14e-6) * si::metres;
        this->sdev_rd1 = real_t(1.28),
        this->sdev_rd2 = real_t(1.75);
        this->n1_stp = real_t(90e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(15e6) / si::cubic_metres;  // 65 || 16
        this->ForceParameters.coriolis_parameter = 0.; // [1/s]
        this->z_rlx = z_rlx;
      }

      public:

        // TODO: we can delete the "window" option?
      // ctor
      DaytimeConvection_common(const real_t _X, const real_t _Y, const real_t _Z, const bool window)
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
    class DaytimeConvection;


    template<class case_ct_params_t>
    class DaytimeConvection<case_ct_params_t, 2> : public DaytimeConvection_common<case_ct_params_t, 2>
    {
      using parent_t = DaytimeConvection_common<case_ct_params_t, 2>;
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
    class DaytimeConvection<case_ct_params_t, 3> : public DaytimeConvection_common<case_ct_params_t, 3>
    {
      using parent_t = DaytimeConvection_common<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // southerly wind
      struct v_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        v_t() : hori_vel_t(&v_lba) {}

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
        // geostrophic wind equal to the initial velocity profile, doesnt include mean, because its only used in coriolis = u-geostr
        blitz::firstIndex k;
        real_t dz = (this->Z / si::metres) / (nz-1);
        //profs.geostr[0] = this->u(k * dz);
        //profs.geostr[1] = v(k * dz);
      }

      public:
      DaytimeConvection(const real_t _X, const real_t _Y, const real_t _Z, const bool window):
        parent_t(_X, _Y, _Z, window)
        {
          v.init(window, this->Z);
          this->ForceParameters.uv_mean[1] = v.mean_vel;
        }
    };
  };
};

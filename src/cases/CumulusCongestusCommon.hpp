#pragma once
#include <random>
#include <fstream>
#include "Anelastic.hpp"

namespace cases 
{
  namespace CumulusCongestus
  {
    template<class case_ct_params_t, int n_dims>
    class CumulusCongestusCommon : public Anelastic<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

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
        params.fricvelsq = 0.0784;
      }
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
        params.buoyancy_wet = true;
        params.subsidence = false;
        params.vel_subsidence = false;
        params.friction = true;
        params.coriolis = false;
        params.radiation = false;

        this->setopts_sgs(params);
      }
  
      // RH T and p to rv
      quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
      {
        return RH * libcloudphxx::common::const_cp::r_vs<real_t>(T, p);
      }
  
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, int rng_seed, index_t index, 
                        const real_t &th_prtrb, // [K]
                        const real_t &rv_prtrb, // [kg/kg]
                        const real_t &z_abs)    // [m]
      {
        // we assume here that set_profs was called already, so that *_env profiles are initialized
        int nz = concurr.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (this->Z / si::metres) / (nz-1); 
        // copy the env profiles into 2D/3D arrays
        //concurr.advectee(ix::rv) = rv_env(index); 
        //concurr.advectee(ix::th) = th_std_env(index); 
  
        concurr.advectee(ix::w) = 0;  
       
        // absorbers
        concurr.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (this->Z / si::metres - z_abs)), 2), 0);
        concurr.vab_relaxed_state(0) = concurr.advectee(ix::u);
        concurr.vab_relaxed_state(ix::w) = concurr.advectee(ix::w);
  
        // density profile
        concurr.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
  
        // randomly prtrb tht and rv in the lowest 1km
        // NOTE: all processes do this, but ultimately only perturbation calculated by MPI rank 0 is used
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-1*th_prtrb, th_prtrb);
          auto rand = std::bind(dis, gen);
  
          auto th_global = concurr.advectee_global(ix::th);
          decltype(concurr.advectee(ix::th)) prtrb(th_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?

          prtrb = where(index * dz >= 1000., 0., prtrb); // no perturbation above 1km
          th_global += prtrb;
          this->make_cyclic(th_global);
          concurr.advectee_global_set(th_global, ix::th);
        }
        {
          std::mt19937 gen(rng_seed+1); // different seed than in th. NOTE: if the same instance of gen is used in th and rv, for some reason it gives the same sequence in rv as in th despite being advanced in th prtrb
          std::uniform_real_distribution<> dis(-1*rv_prtrb, rv_prtrb);
          auto rand = std::bind(dis, gen);
  
          auto rv_global = concurr.advectee_global(ix::rv);
          decltype(concurr.advectee(ix::rv)) prtrb(rv_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?

          prtrb = where(index * dz >= 1000., 0., prtrb); // no perturbation above 1km
          rv_global += prtrb;
          this->make_cyclic(rv_global);
          concurr.advectee_global_set(rv_global, ix::rv);
        }
      }
  

      // functions that set surface fluxes per timestep
      void update_surf_flux_sens_hlpr(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy,
                                       const real_t &spinup_flux,      // [K m/s] 
                                       const real_t &gauss_half_width, // [m]    
                                       const real_t &gauss_max_flux    // [K m/s]
                                     ) 
      {
        if(timestep == 0) 
          surf_flux_sens = spinup_flux * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters) * theta_std::exner(this->p_0); // [K kg / (m^2 s)]; -1 because negative gradient of upward flux means inflow
        else if(int((3600. / dt) + 0.5) == timestep)
        {
          if(surf_flux_sens.rank() == 3) // TODO: make it a compile-time decision
            surf_flux_sens = gauss_max_flux * exp( - ( pow(blitz::tensor::i * dx - real_t(0.5) * this->X / si::metres, 2) +  pow(blitz::tensor::j * dy - real_t(0.5) * this->Y / si::metres, 2) ) / (gauss_half_width * gauss_half_width) ) * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters) * theta_std::exner(this->p_0);
          else if(surf_flux_sens.rank() == 2)
            surf_flux_sens = gauss_max_flux * exp( - ( pow(blitz::tensor::i * dx - real_t(0.5) * this->X / si::metres, 2)  ) / (gauss_half_width * gauss_half_width) ) * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters) * theta_std::exner(this->p_0);
        }
      }
      

      void update_surf_flux_lat_hlpr(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy,
                                       const real_t &spinup_flux,      // [kg/kg m/s] 
                                       const real_t &gauss_half_width, // [m]         
                                       const real_t &gauss_max_flux)   // [kg/kg m/s] 
      {
        if(timestep == 0)
          surf_flux_lat = spinup_flux * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters); // [m/s]
        else if(int((3600. / dt) + 0.5) == timestep)
        {
          if(surf_flux_lat.rank() == 3) // TODO: make it a compile-time decision
            surf_flux_lat = gauss_max_flux * exp( - ( pow(blitz::tensor::i * dx - real_t(0.5) * this->X / si::metres, 2) +  pow(blitz::tensor::j * dy - real_t(0.5) * this->Y / si::metres, 2) ) / (gauss_half_width * gauss_half_width) ) * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters);
          else if(surf_flux_lat.rank() == 2)
            surf_flux_lat = gauss_max_flux * exp( - ( pow(blitz::tensor::i * dx - real_t(0.5) * this->X / si::metres, 2)  ) / (gauss_half_width * gauss_half_width) ) * -1 * (this->rhod_0 / si::kilograms * si::cubic_meters);
        }
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind
                               const real_t &U_ground_z,
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy, const real_t &uv_mean) override
      {
        surf_flux_uv = where(U_ground < 1e-4, 
            - 0.0784 * (uv_ground + uv_mean) / real_t(1e-4) * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters), // 0.0784 m^2 / s^2 is the square of friction velocity = 0.28 m / s
            - 0.0784 * (uv_ground + uv_mean) / U_ground * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters)
          );
      }
    };
  };
};

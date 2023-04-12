#pragma once
#include <random>
#include <type_traits>
#include "Anelastic.hpp"

template< class, class = std::void_t<> >
struct has_SO2g : std::false_type { };

template< class T >
struct has_SO2g<T, std::void_t<typename T::SO2g>> : std::true_type { };

namespace cases 
{
  namespace dycoms
  {
    namespace hydrostatic = libcloudphxx::common::hydrostatic;
    namespace theta_std = libcloudphxx::common::theta_std;
    namespace theta_dry = libcloudphxx::common::theta_dry;
    namespace lognormal = libcloudphxx::common::lognormal;

    const quantity<si::pressure, real_t> p_0 = 101780 * si::pascals;
    const quantity<si::length, real_t> 
      z_0   = 0    * si::metres,
      Z_def = 1500 * si::metres;
    const quantity<si::length, real_t> X_def[] = {/*RF1*/3360 * si::metres, /*RF2*/6400 * si::metres};
    const quantity<si::length, real_t> Y_def[] = {/*RF1*/3360 * si::metres, /*RF2*/6400 * si::metres};
    const real_t z_abs = 1250;
    const real_t z_i[] = {/*RF1*/840, /*RF2*/795}; //initial inversion height
    const quantity<si::length, real_t> z_rlx = 25 * si::metres;
    const quantity<si::length, real_t> gccn_max_height = 450 * si::metres; // below cloud
    const quantity<si::frequency, real_t> D = real_t(3.75e-6) / si::seconds; // large-scale wind horizontal divergence

    template <int RF>
    quantity<si::velocity, real_t> u_dycoms(const real_t &z)
    {
      return RF == 1 ? 
        real_t(7) * si::meters / si::seconds :                    // RF01 
        real_t(3. + 4.3 * z / 1000.) * si::meters / si::seconds;  // RF02
    }
    template <int RF>
    quantity<si::velocity, real_t> v_dycoms(const real_t &z)
    {
      return RF == 1 ? 
        real_t(-5.5) * si::meters / si::seconds :                 // RF01 
        real_t(-9 + 5.6 * z / 1000.) * si::meters / si::seconds;  // RF02
    }

    // liquid water potential temperature at height z
    template <int RF>
    quantity<si::temperature, real_t> th_l_dycoms(const real_t &z)
    {
      const quantity<si::temperature, real_t>
        th_below = real_t(RF == 1 ? 289 : 288.3) * si::kelvins;
      const real_t th_above_hlpr = (RF == 1 ? 297.5 : 295) + pow(z - z_i[RF - 1], real_t(1./3));
      const quantity<si::temperature, real_t> th_above = th_above_hlpr * si::kelvins; 
      return z < z_i[RF - 1] ? th_below : th_above;
    }

    template <int RF>
    quantity<si::dimensionless, real_t> r_t_dycoms(const real_t &z)
    {
      const quantity<si::dimensionless, real_t>
        rt_below = RF == 1 ? 9.0e-3 : 9.45e-3;
      const quantity<si::dimensionless, real_t>
        rt_above = RF == 1 ? 1.5e-3 : (5. - 3. * (1. - exp((z_i[RF - 1] - z)/500.))) * 1e-3;
      return z < z_i[RF - 1] ? rt_below : rt_above;
    }

    inline quantity<si::dimensionless, real_t> SO2g_dycoms(const real_t &z)
    {
      return quantity<si::dimensionless, real_t>(5e-10);
    }
  
    template<class case_ct_params_t, int RF, int n_dims>
    class DycomsCommon : public Anelastic<case_ct_params_t, n_dims>
    {
      static_assert(RF == 1 || RF == 2,
                    "only setups based on the first and the second DYCOMS research flights are available");

      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return th_l_dycoms<RF>(z);
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return r_t_dycoms<RF>(z);
      }
    
      // water mixing ratio at height z
      struct r_t_fctr
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return r_t_dycoms<RF>(z);
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };
    
      // initial standard potential temp at height z, assuming theta_std = theta_l (spinup needed)
      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return th_l_dycoms<RF>(z) / si::kelvins;
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

        u_t() : hori_vel_t(&u_dycoms<RF>) {}

        BZ_DECLARE_FUNCTOR(u_t);
      };

      u_t u;
    
      // large-scale vertical wind
      struct w_LS_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return - (D * si::seconds) * z; 
        }
        BZ_DECLARE_FUNCTOR(w_LS_fctr);
      };

      // initial ambient SO2 profile (mixing ratio?)
      struct SO2g_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return SO2g_dycoms(z);
        }
        BZ_DECLARE_FUNCTOR(SO2g_fctr);
      };
    
      // density profile as a function of altitude
      // hydrostatic and assuming constant theta (not used now)
      //struct rhod_fctr
      //{
      //  real_t operator()(real_t z) const
      //  {
      //    quantity<si::pressure, real_t> p = hydrostatic::p(
      //    z * si::metres, th_dry_fctr()(0.) * si::kelvins, r_t()(0.), z_0, p_0
      //    );
      //    
      //    quantity<si::mass_density, real_t> rhod = theta_std::rhod(
      //    p, th_dry_fctr()(0.) * si::kelvins, r_t()(0.)
      //    );
    
      //    return rhod / si::kilograms * si::cubic_metres;
      //  }
    
      //  // to make the rhod() functor accept Blitz arrays as arguments
      //  BZ_DECLARE_FUNCTOR(rhod_fctr);
      //};
      //

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
        params.fricvelsq = 0.0625;
      }
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
//        params.outdir = user_params.outdir;
//        params.outfreq = user_params.outfreq;
//        params.spinup = user_params.spinup;
//        params.w_src = user_params.w_src;
//        params.uv_src = user_params.uv_src;
//        params.th_src = user_params.th_src;
//        params.rv_src = user_params.rv_src;
//        params.rc_src = user_params.rc_src;
//        params.rr_src = user_params.rr_src;
//        params.nc_src = user_params.nc_src;
//        params.nr_src = user_params.nr_src;
//        params.dt = user_params.dt;
//        params.nt = user_params.nt;
//        params.relax_th_rv = user_params.relax_th_rv;
        params.buoyancy_wet = true;
        params.subsidence = true;
        params.vel_subsidence = true;
        params.friction = true;
        params.coriolis = true;
        params.radiation = true;

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
        concurr.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (this->Z / si::metres - z_abs)), 2), 0);
        concurr.vab_relaxed_state(0) = concurr.advectee(ix::u);
        concurr.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
  
        // density profile
        concurr.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
  
        // initial potential temperature
        concurr.advectee(ix::th) = th_std_fctr()(index * dz); 

        // randomly prtrb tht
        // NOTE: all processes do this, but ultimately only perturbation calculated by MPI rank 0 is used
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-0.1, 0.1);
          auto rand = std::bind(dis, gen);
  
          auto th_global = concurr.advectee_global(ix::th);
          decltype(concurr.advectee(ix::th)) prtrb(th_global.shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          th_global += prtrb;
          this->make_cyclic(th_global);
          concurr.advectee_global_set(th_global, ix::th);
        }

        // chemical composition
        if constexpr(has_SO2g<ix>{}())
          concurr.advectee(ix::SO2g) = SO2g_fctr()(index * dz);

    //      config::mixr_helper(this->setup)(index * dz)
    //      * (SO2_g_0 * molar_mass::M_SO2<real_t>()  * si::moles / si::kilograms);

      }
  
      // calculate the initial environmental theta and rv profiles
      // like in Wojtek's BabyEulag
      // alse set w_LS and hgt_fctrs
      // TODO: move hgt_fctrs from cases to main code
      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params) override
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
        profs.th_LS = 0.; // no large-scale horizontal advection
        profs.rv_LS = 0.; 

        //nudging, Zhou et al. 2018
        profs.relax_th_rv_coeff = where(k * dz >= 800, 
          1. / 180.,  // 180s time scale at and above 800m
          1. / 7200. * pow(sin(3.1419 / 2. * (k * dz) / 800.), 2) // 7200s time scale below 800m + sinusoidal factor = 0 at ground
          );
      }

      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
        {
          auto flux_value = RF == 1 ? 15. : 16.; // [W/m^2]
          auto conv_fctr_sens = (libcloudphxx::common::moist_air::c_pd<real_t>() * si::kilograms * si::kelvins / si::joules);
          surf_flux_sens = -flux_value / conv_fctr_sens; // [K * kg / (m^2 * s)]
        }
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
        {
          auto flux_value = RF == 1 ? 115. : 93.; // [W/m^2]
          auto conv_fctr_lat = (libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules);
          surf_flux_lat = -flux_value / conv_fctr_lat; // [kg / (m^2 * s)]
        }
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind
                               const real_t &U_ground_z,
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy, const real_t &uv_mean)
      {
        surf_flux_uv = where(U_ground == 0., 0.,
            - 0.0625 * (uv_ground + uv_mean) / U_ground * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters)// 0.0625 m^2 / s^2 is the square of friction velocity = 0.25 m / s; * -1 because negative gradient of upward flux means inflow
          );
      }

      void init() 
      {
        //aerosol bimodal lognormal dist. - DYCOMS
        this->p_0 = p_0;
        this->mean_rd1 = real_t(.011e-6) * si::metres,
        this->mean_rd2 = real_t(.06e-6) * si::metres;
        this->sdev_rd1 = real_t(1.2),
        this->sdev_rd2 = real_t(1.7);
        this->n1_stp = real_t(125e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(65e6) / si::cubic_metres;  // 65 || 16
        this->ForceParameters.coriolis_parameter = 0.76e-4; // [1/s] @ 31.5 deg N
        this->ForceParameters.D = D * si::seconds; 
        this->z_rlx = z_rlx;
        this->gccn_max_height = gccn_max_height;
      }


      public:
      // ctor
      DycomsCommon(const real_t _X, const real_t _Y, const real_t _Z, const bool window)
      {
        init();

        this->X = _X < 0 ? X_def[RF-1] : _X * si::meters;
        if(n_dims == 3)
          this->Y = _Y < 0 ? Y_def[RF-1] : _Y * si::meters;
        this->Z = _Z < 0 ? Z_def : _Z * si::meters;
        u.init(window, this->Z);

        this->ForceParameters.uv_mean[0] = u.mean_vel;
      }
    };
    
    template<class case_ct_params_t, int RF, int n_dims>
    class Dycoms;

    template<class case_ct_params_t, int RF>
    class Dycoms<case_ct_params_t, RF, 2> : public DycomsCommon<case_ct_params_t, RF, 2>
    {
      using parent_t = DycomsCommon<case_ct_params_t, RF, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);
      };

      // ctor
      using parent_t::parent_t;
    };

    template<class case_ct_params_t, int RF>
    class Dycoms<case_ct_params_t, RF, 3> : public DycomsCommon<case_ct_params_t, RF, 3>
    {
      using parent_t = DycomsCommon<case_ct_params_t, RF, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // southerly wind
      struct v_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        v_t() : hori_vel_t(&v_dycoms<RF>) {}

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

      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
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
        // geostrophic wind equal to the initial velocity profile
        blitz::firstIndex k;
        real_t dz = (this->Z / si::metres) / (nz-1);
        profs.geostr[0] = this->u(k * dz); 
        profs.geostr[1] = v(k * dz); 
      }

      public:
      Dycoms(const real_t _X, const real_t _Y, const real_t _Z, const bool window):
        parent_t(_X, _Y, _Z, window)
        {
          v.init(window, this->Z);
          this->ForceParameters.uv_mean[1] = v.mean_vel;
        }
    };
  };
};

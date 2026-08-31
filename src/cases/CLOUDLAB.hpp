#pragma once
#include <random>
#include <fstream>
#include <libcloudph++/lgrngn/opts.hpp>
#include <libcloudph++/lgrngn/ccn_source.hpp>
#include <libcloudph++/common/ice_nucleation.hpp>
#include "detail/CLOUDLAB_sounding/CLOUDLAB_init_data.hpp"
#include "Anelastic.hpp"

namespace cases 
{
  namespace CLOUDLAB
  {
    const quantity<si::pressure, real_t> p_0 = real_t(CLOUDLAB_pres_sfc) * si::pascals;
    const quantity<si::length, real_t>   Z_def = 1200 * si::metres;
    const quantity<si::length, real_t>   X_def = 3000 * si::metres;
    const quantity<si::length, real_t>   Y_def = 3000 * si::metres;
    const real_t z_abs = 1000;
    constexpr real_t seeding_start_time = 2500;
    constexpr real_t seeding_end_time = seeding_start_time + 6 * 60;
    
    constexpr real_t fricvelsq = 0.25; // forest roughness length of 0.5 m, u* = 0.5 m/s
    
    // returned units: [K], [kg/kg], [m/s]
    inline real_t interpolate_CLOUDLAB_sounding(const std::string valname, real_t pos)
    {
      assert(pos>=0);
      assert(valname == "th_l" || valname == "qt" || valname == "u" || valname == "v");

      const auto &pos_s(CLOUDLAB_profile_height);

      const auto pos_up = 
        // valname == "aerosol_conc_factor" ? std::upper_bound(pos_s.begin(), pos_s.end(), pos, std::greater<double>()) : 
                                           std::upper_bound(pos_s.begin(), pos_s.end(), pos);

      if(pos_up == pos_s.end())
        throw std::runtime_error("UWLCM: The initial sounding is not high enough");
      if(pos_up == pos_s.begin())
        throw std::runtime_error("UWLCM: The initial sounding has an incorrect first element (?)");

      const auto &sndg(
        valname == "th_l" ? CLOUDLAB_profile_th_l :
        valname == "qt" ? CLOUDLAB_profile_qt :
        valname == "u" ? CLOUDLAB_profile_u :
        CLOUDLAB_profile_v);

      const auto s_up = sndg.begin() + std::distance(pos_s.begin(), pos_up);
      return real_t(*(s_up-1) + (pos - *(pos_up-1)) / (*pos_up - *(pos_up-1)) * (*s_up - *(s_up-1)));
    }

    // we rotate coordinate system to have v=0 and u>0 at 350 m (seeding height), to make it easier to seed perpendicular to the wind
    // note: this breaks the geostrophic wind assumption, but we don't use it anyway, so it's fine
    // note2: setting window=1 ruins v=0 at 350 m, but seeding is still perpendicular to the wind, so it's fine (?)
    const real_t coord_rotation_angle = atan2(interpolate_CLOUDLAB_sounding("v", 350.), interpolate_CLOUDLAB_sounding("u", 350.)); 

    inline quantity<si::velocity, real_t> u_clab(const real_t &z)
    {
      return real_t(interpolate_CLOUDLAB_sounding("u", z) * cos(coord_rotation_angle) + interpolate_CLOUDLAB_sounding("v", z) * sin(coord_rotation_angle)) * si::meters / si::seconds;
    }
    
    inline quantity<si::velocity, real_t> v_clab(const real_t &z)
    {
      return real_t(-interpolate_CLOUDLAB_sounding("u", z) * sin(coord_rotation_angle) + interpolate_CLOUDLAB_sounding("v", z) * cos(coord_rotation_angle)) * si::meters / si::seconds;
    }

    inline quantity<si::temperature, real_t> th_l_clab(const real_t &z)
    {
      return interpolate_CLOUDLAB_sounding("th_l", z) * si::kelvins;
    }

    inline quantity<si::dimensionless, real_t> r_t_clab(const real_t &z)
    {
      return interpolate_CLOUDLAB_sounding("qt", z);
    }

    template<class case_ct_params_t, int n_dims>
    class CloudlabCommon : public Anelastic<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return th_l_clab(z);
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return r_t_clab(z);
      }

      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return th_l_clab(z) / si::kelvins;
        }
        BZ_DECLARE_FUNCTOR(th_std_fctr);
      };

      struct r_t_fctr 
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return r_t_clab(z);
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };

      struct u_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        u_t() : hori_vel_t(&u_clab) {}

        BZ_DECLARE_FUNCTOR(u_t);
      };

      u_t u;


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
        params.fricvelsq = fricvelsq;
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
        params.subsidence = subs_t::none;
        params.vel_subsidence = true;
        params.friction = true;
        params.coriolis = false;
        params.radiation = false;

        this->setopts_sgs(params);
        setopts_lgrngn_hlpr(params);
      }

      template <class T>
      void setopts_lgrngn_hlpr(
        T &params,
        typename std::enable_if<std::is_same<
          decltype(T::cloudph_opts),
          libcloudphxx::lgrngn::opts_t<real_t>
        >::value>::type* = 0
      )
      {
        params.cloudph_opts_init.inp_type = libcloudphxx::common::ice_nucleation::INP_t::AgI;
        params.cloudph_opts_init.src_type = libcloudphxx::lgrngn::src_t::simple;
        params.cloudph_opts_init.src_x0 = 1400;
        params.cloudph_opts_init.src_x1 = 1500;
        params.cloudph_opts_init.src_y0 = 1300;
        params.cloudph_opts_init.src_y1 = 1700;
        params.cloudph_opts_init.src_z0 = 340;
        params.cloudph_opts_init.src_z1 = 360;

        params.cloudph_opts.src_dry_sizes.emplace(
          libcloudphxx::lgrngn::kappa_soluble_fraction_t<real_t>(real_t(1.2), real_t(0.5)),
          std::map<real_t, std::tuple<real_t, int, int>>{
            {real_t(0.252e-6), {real_t(1e5), 1, 1}}
          }
        );
      }

      template <class T>
      void setopts_lgrngn_hlpr(
        T &,
        typename std::enable_if<!std::is_same<
          decltype(T::cloudph_opts),
          libcloudphxx::lgrngn::opts_t<real_t>
        >::value>::type* = 0
      )
      {}

      template <class T>
      void setopts_ante_step_hlpr(
        T &params,
        const int timestep,
        typename std::enable_if<std::is_same<
          decltype(T::cloudph_opts),
          libcloudphxx::lgrngn::opts_t<real_t>
        >::value>::type* = 0
      )
      {
        const real_t time = timestep * params.dt;
        params.cloudph_opts.src = time >= seeding_start_time && time < seeding_end_time;
      }

      template <class T>
      void setopts_ante_step_hlpr(
        T &,
        const int,
        typename std::enable_if<!std::is_same<
          decltype(T::cloudph_opts),
          libcloudphxx::lgrngn::opts_t<real_t>
        >::value>::type* = 0
      )
      {}

      void setopts(rt_params_t &params, const int timestep) override
      {
        setopts_ante_step_hlpr(params, timestep);
      }
  

      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        int nz = rhod.extent(0) - 1;
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
      }  

      // calculate the initial environmental theta and rv profiles
      // alse set w_LS and hgt_fctrs
      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
      {
        parent_t::set_profs(profs, nz, user_params);

        this->env_prof(profs, nz);
        this->ref_prof(profs, nz);

        profs.w_LS = 0.; // no subsidence
        profs.th_LS = 0.; // no large-scale horizontal advection
        profs.rv_LS = 0.; 

        //nudging, Zhou et al. 2018
        // blitz::firstIndex k;
        // real_t dz = (this->Z / si::metres) / (nz-1);
        // profs.relax_th_rv_coeff = where(k * dz >= 800, 
        //   1. / 180.,  // 180s time scale at and above 800m
        //   1. / 7200. * pow(sin(3.1419 / 2. * (k * dz) / 800.), 2) // 7200s time scale below 800m + sinusoidal factor = 0 at ground
        //   );
      }

      // functions that set surface fluxes per timestep
      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
        {
          auto flux_value = 0; // [W/m^2]
          auto conv_fctr_sens = (libcloudphxx::common::moist_air::c_pd<real_t>() * si::kilograms * si::kelvins / si::joules);
          surf_flux_sens = -flux_value / conv_fctr_sens; // [K * kg / (m^2 * s)]
        }
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
        {
          auto flux_value = 0; // [W/m^2]
          auto conv_fctr_lat = (libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules);
          surf_flux_sens = -flux_value / conv_fctr_lat; // [kg / (m^2 * s)]
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
        surf_flux_uv = where(U_ground == 0., 0.,
            - fricvelsq * (uv_ground + uv_mean) / U_ground * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters)
          );
      }

      void init()
      {
        this->p_0 = p_0;
        //aerosol bimodal lognormal dist. - as in RICO, but 4x conc
        this->mean_rd1 = real_t(.03e-6) * si::metres,
        this->mean_rd2 = real_t(.14e-6) * si::metres;
        this->sdev_rd1 = real_t(1.28),
        this->sdev_rd2 = real_t(1.75);
        this->n1_stp = real_t(3 * 90e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(3 * 15e6) / si::cubic_metres;  // 65 || 16
        this->z_rlx = real_t(1e2) * si::metres;
      }

      public:
      // ctor
      CloudlabCommon(const real_t _X, const real_t _Y, const real_t _Z, const bool window)
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
    class Cloudlab;

    template<class case_ct_params_t>
    class Cloudlab<case_ct_params_t, 2> : public CloudlabCommon<case_ct_params_t, 2>
    {
      using parent_t = CloudlabCommon<case_ct_params_t, 2>;
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
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed, const int nps[2]) override
      {
        blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);
      };

      // ctor
      using parent_t::parent_t;
    };

    template<class case_ct_params_t>
    class Cloudlab<case_ct_params_t, 3> : public CloudlabCommon<case_ct_params_t, 3>
    {
      using parent_t = CloudlabCommon<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // southerly wind
      struct v_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        v_t() : hori_vel_t(&v_clab) {}

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
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed, const int nps[3]) override
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);

        int nz = nps[2];
        real_t dz = (this->Z / si::metres) / (nz-1);

        concurr.advectee(ix::v)= v(k * dz);
        concurr.vab_relaxed_state(1) = concurr.advectee(ix::v);
      }

      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
      {
        parent_t::set_profs(profs, nz, user_params);
      }

      public:
      Cloudlab(const real_t _X, const real_t _Y, const real_t _Z, const bool window):
        parent_t(_X, _Y, _Z, window)
        {
          v.init(window, this->Z);
          this->ForceParameters.uv_mean[1] = v.mean_vel;
        }
    };
  };
};

// Cumulus Congestus case from International Cloud Modeling Workshop 2020
// However, some changes were made (?) Actual setup as in DOI: 10.5194/gmd-17-759-2024 

#pragma once
#include <random>
#include <fstream>
#include "CumulusCongestusCommon.hpp"
#include "detail/CAMP2EX_sounding/input_sounding_camp2ex_000907-pbl4.hpp"
#include "detail/CAMP2EX_sounding/aerosol_profile_factor.hpp"

namespace cases 
{
  namespace CumulusCongestus
  {
    // returned units: [K], [g/kg], [m/s], [1]
    inline real_t interpolate_CAMP2EX_sounding(const std::string valname, real_t pos)
    {
      assert(pos>=0);
      assert(valname == "theta" || valname == "rv" || valname == "u" || valname == "v" || valname == "aerosol_conc_factor");

      const auto &pos_s(
        valname == "aerosol_conc_factor" ? CAMP2EX_aerosol_profile_density : CAMP2EX_sounding_z);

      const auto pos_up = 
        valname == "aerosol_conc_factor" ? std::upper_bound(pos_s.begin(), pos_s.end(), pos, std::greater<double>()) : 
                                           std::upper_bound(pos_s.begin(), pos_s.end(), pos);

      if(pos_up == pos_s.end())
        throw std::runtime_error("UWLCM: The initial sounding is not high enough");
      if(pos_up == pos_s.begin())
        throw std::runtime_error("UWLCM: The initial hsa incorrect first element (?)");

      const auto &sndg(
        valname == "theta" ? CAMP2EX_sounding_theta :
        valname == "rv" ? CAMP2EX_sounding_rv :
        valname == "u" ? CAMP2EX_sounding_u :
        valname == "v" ? CAMP2EX_sounding_v :
        CAMP2EX_aerosol_profile_factor);

      const auto s_up = sndg.begin() + std::distance(pos_s.begin(), pos_up);
      return real_t(*(s_up-1) + (pos - *(pos_up-1)) / (*pos_up - *(pos_up-1)) * (*s_up - *(s_up-1)));
    }

    inline quantity<si::velocity, real_t> u_cc_icmw24(const real_t &z)
    {
      return interpolate_CAMP2EX_sounding("u", z) * si::meters / si::seconds;
      //return 0 * si::meters / si::seconds;
    }

    inline quantity<si::velocity, real_t> v_cc_icmw24(const real_t &z)
    {
      return interpolate_CAMP2EX_sounding("v", z) * si::meters / si::seconds;
      //return 0 * si::meters / si::seconds;
    }

    inline quantity<si::temperature, real_t> th_l_cc_icmw24(const real_t &z)
    {
      return interpolate_CAMP2EX_sounding("theta", z) * si::kelvins;
    }

    inline quantity<si::dimensionless, real_t> r_t_cc_icmw24(const real_t &z)
    {
      return interpolate_CAMP2EX_sounding("rv", z) * 1e-3; // to [kg/kg]
    }

    inline quantity<si::dimensionless, real_t> acf_cc_icmw24(const real_t &rhod) // aerosol concentration factor
    {
      return interpolate_CAMP2EX_sounding("aerosol_conc_factor", rhod);
    }

    template<class case_ct_params_t, int n_dims>
    class CumulusCongestusCommon_icmw24 : public CumulusCongestusCommon<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = CumulusCongestusCommon<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      quantity<si::temperature, real_t> th_l(const real_t &z) override
      {
        return th_l_cc_icmw24(z);
      }

      quantity<si::dimensionless, real_t> r_t(const real_t &z) override
      {
        return r_t_cc_icmw24(z);
      }

      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return th_l_cc_icmw24(z) / si::kelvins;
        }
        BZ_DECLARE_FUNCTOR(th_std_fctr);
      };

      struct r_t_fctr
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return r_t_cc_icmw24(z);
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };

      struct u_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        u_t() : hori_vel_t(&u_cc_icmw24) {}

        BZ_DECLARE_FUNCTOR(u_t);
      };

      u_t u;

      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, int rng_seed, index_t index) 
      {
        int nz = concurr.advectee_global().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (this->Z / si::metres) / (nz-1);

        concurr.advectee(ix::th) = th_std_fctr{}(index * dz);
        concurr.advectee(ix::rv) = r_t_fctr{}(index * dz);
        concurr.advectee(ix::u) = u(index * dz);

        parent_t::intcond_hlpr(concurr, rhod, rng_seed, index, 0.1, 0.025e-3, (this->Z / si::metres) - 1000);
      }

      template <class T, class U>
      void setopts_hlpr(T &params, const int nz, const U &user_params)
      {
        params.aerosol_independent_of_rhod=true;
        for(int i=0; i<nz; ++i)
        {
          params.aerosol_conc_factor.push_back(acf_cc_icmw24((*params.rhod)(i)));
        }
        parent_t::setopts_hlpr(params, user_params);
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
      }

      // functions that set surface fluxes per timestep
      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        parent_t::update_surf_flux_sens_hlpr(surf_flux_sens, th_ground, U_ground, U_ground_z, timestep, dt, dx, dy, 0.01, 3000, .3, true);
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        parent_t::update_surf_flux_lat_hlpr(surf_flux_sens, th_ground, U_ground, U_ground_z, timestep, dt, dx, dy, 4e-5, 3000, 2.4e-4, true);
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

      void init()
      {
        this->p_0 = real_t(101800) * si::pascals;
        //aerosol bimodal lognormal dist. - as in RICO with 11x conc following the ICMW2020 setup
        this->mean_rd1 = real_t(.09e-6) * si::metres,
        this->mean_rd2 = real_t(.5e-6) * si::metres;
        this->sdev_rd1 = real_t(1.65),
        this->sdev_rd2 = real_t(1.65);
        this->n1_stp = real_t(680e6) / si::cubic_metres, 
        this->n2_stp = real_t(2.24e6) / si::cubic_metres;
        this->z_rlx = real_t(1e2) * si::metres;
      }

      // ctor
      CumulusCongestusCommon_icmw24(const real_t _X, const real_t _Y, const real_t _Z, const bool window)
      {
        init();

        this->Z = _Z < 0 ? 11e3 * si::meters : _Z * si::meters;

        u.init(window, this->Z);
        this->ForceParameters.uv_mean[0] = u.mean_vel;
      }
    };
    
    template<class case_ct_params_t, int n_dims>
    class CumulusCongestus_icmw24;

    template<class case_ct_params_t>
    class CumulusCongestus_icmw24<case_ct_params_t, 2> : public CumulusCongestusCommon_icmw24<case_ct_params_t, 2>
    {
      using parent_t = CumulusCongestusCommon_icmw24<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, nps[1], user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);
      }

      public:
      CumulusCongestus_icmw24(const real_t _X, const real_t _Y, const real_t _Z, const bool window):
        parent_t(_X, _Y, _Z, window)
      {
        this->X = _X < 0 ? 12e3 * si::meters : _X * si::meters;
      }
    };

    template<class case_ct_params_t>
    class CumulusCongestus_icmw24<case_ct_params_t, 3> : public CumulusCongestusCommon_icmw24<case_ct_params_t, 3>
    {
      using parent_t = CumulusCongestusCommon_icmw24<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      struct v_t : hori_vel_t
      {
        real_t operator()(const real_t &z) const
        {
          return hori_vel_t::operator()(z);
        }

        v_t() : hori_vel_t(&v_cc_icmw24) {}

        BZ_DECLARE_FUNCTOR(v_t);
      };

      v_t v;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, nps[2], user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Y / si::metres) / (nps[1]-1);
        params.dk = (this->Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        int nz = concurr.advectee_global().extent(ix::w);
        real_t dz = (this->Z / si::metres) / (nz-1); 
  
        concurr.advectee(ix::v)= v(k * dz);
        concurr.vab_relaxed_state(1) = concurr.advectee(ix::v);

        this->intcond_hlpr(concurr, rhod, rng_seed, k);
      }

      public:
      CumulusCongestus_icmw24(const real_t _X, const real_t _Y, const real_t _Z, const bool window):
        parent_t(_X, _Y, _Z, window)
      {
        this->X = _X < 0 ? 12e3 * si::meters : _X * si::meters;
        this->Y = _Y < 0 ? 12e3 * si::meters : _Y * si::meters;
        v.init(window, this->Z);
        this->ForceParameters.uv_mean[1] = v.mean_vel;
      }
    };
  };
};

  


#pragma once

#include <iostream>

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>

#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

#include "../detail/user_params.hpp"
#include "../detail/concurr_types.hpp"

namespace setup 
{
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;

  // container for constants that appear in forcings, some are not needed in all cases, etc...
  // TODO: make forcing functions part of case class
  struct ForceParameters_t
  {
    real_t q_i, heating_kappa, F_0, F_1, rho_i, D, coriolis_parameter;
    bool surf_latent_flux_in_watts_per_square_meter;
    bool surf_sensible_flux_in_watts_per_square_meter;
  };
 
  // CAUTION: new profiles have to be added to both structs and in copy_profiles below
  // TODO: try a different design where it is not necessary ?
  struct profiles_t
  {
    arr_1D_t th_e, p_e, rv_e, rl_e, th_ref, rhod, w_LS, hgt_fctr, th_LS, rv_LS, mix_len;
    std::array<arr_1D_t, 2> geostr;

    profiles_t(int nz) :
    // rhod needs to be bigger, cause it divides vertical courant number
    // TODO: should have a halo both up and down, not only up like now; then it should be interpolated in courant calculation
      th_e(nz), p_e(nz), rv_e(nz), rl_e(nz), th_ref(nz), rhod(nz+1), w_LS(nz), hgt_fctr(nz), th_LS(nz), rv_LS(nz), mix_len(nz)
    {
      geostr[0].resize(nz);
      geostr[1].resize(nz);

      // set to zero just to have predicatble output in cases without Coriolis
      geostr[0] = 0.;
      geostr[1] = 0.;
    }
  };
  struct profile_ptrs_t
  {
    arr_1D_t *th_e, *p_e, *rv_e, *rl_e, *th_ref, *rhod, *w_LS, *hgt_fctr, *geostr[2], *th_LS, *rv_LS, *mix_len;
  };

  // copy external profiles into rt_parameters
  // TODO: more elegant way
  template<class params_t>
  inline void copy_profiles(profiles_t &profs, params_t &p)
  {
    std::vector<std::pair<std::reference_wrapper<setup::arr_1D_t*>, std::reference_wrapper<setup::arr_1D_t>>> tobecopied = {
      {p.hgt_fctr     , profs.hgt_fctr     },
      {p.th_e         , profs.th_e         },
      {p.p_e          , profs.p_e          },
      {p.rv_e         , profs.rv_e         },
      {p.rl_e         , profs.rl_e         },
      {p.th_ref       , profs.th_ref       },
      {p.rhod         , profs.rhod         },
      {p.w_LS         , profs.w_LS         },
      {p.th_LS        , profs.th_LS        },
      {p.rv_LS        , profs.rv_LS        },
      {p.geostr[0]    , profs.geostr[0]    },
      {p.geostr[1]    , profs.geostr[1]    },
      {p.mix_len      , profs.mix_len      }
    };

    for (auto dst_src : tobecopied)
    {
      dst_src.first.get() = new setup::arr_1D_t(dst_src.second.get().dataFirst(), dst_src.second.get().shape(), blitz::neverDeleteData);
    }
  }

  template<class case_ct_params_t, int n_dims>
  class CasesCommon
  {
    public:

    using ix = typename case_ct_params_t::ix;
    using rt_params_t = typename case_ct_params_t::rt_params_t;

    using concurr_any_t = libmpdataxx::concurr::any<
      real_t, 
      n_dims
    >;

    //th, rv and surface fluxes relaxation time and height
    const quantity<si::time, real_t> tau_rlx = 300 * si::seconds;

    // domain size
    quantity<si::length, real_t> X,Y,Z;

    //aerosol bimodal lognormal dist. - VOCALS by default
    quantity<si::length, real_t>
      mean_rd1 = real_t(.02e-6) * si::metres,
      mean_rd2 = real_t(.075e-6) * si::metres;
    quantity<si::dimensionless, real_t>
      sdev_rd1 = real_t(1.4),
      sdev_rd2 = real_t(1.6);
    quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
      n1_stp = real_t(70.47e6) / si::cubic_metres, // gives 60e6 at surface of moist thermal
      n2_stp = real_t(46.98e6) / si::cubic_metres;  // gives 40e6 at surface of moist thermal

    // hygroscopicity kappa of the aerosol 
    quantity<si::dimensionless, real_t> kappa = .61; // defaults to ammonium sulphate; CCN-derived value from Table 1 in Petters and Kreidenweis 2007

    real_t div_LS = 0.; // large-scale wind divergence (same as ForceParameters::D), 0. to turn off large-scale subsidence of SDs, TODO: add a process switch in libcloudph++ like for coal/cond/etc

    ForceParameters_t ForceParameters;

    template<bool enable_sgs = case_ct_params_t::enable_sgs>
    void setopts_sgs(rt_params_t &params,
                     typename std::enable_if<!enable_sgs>::type* = 0)
    {}

    template<bool enable_sgs = case_ct_params_t::enable_sgs>
    void setopts_sgs(rt_params_t &params,
                     typename std::enable_if<enable_sgs>::type* = 0)
    {
      params.c_m = 0.0856;
      params.smg_c = 0.165;
      params.prandtl_num = 0.42;
      params.cdrag = 0.;
      params.friction = 0; // disable explicit momentum surface fluxes when using sgs scheme
    }

    virtual void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params) {assert(false);};
    virtual void intcond(concurr_any_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed) =0;
    virtual void set_profs(profiles_t &profs, int nz, const user_params_t &user_params)
    {
      // set SGS mixing length
      {
        blitz::firstIndex k;
        real_t dz = (Z / si::metres) / (nz-1);

        real_t sgs_delta;
        if (user_params.sgs_delta > 0)
        {
          sgs_delta = user_params.sgs_delta;
        }
        else
        {
          sgs_delta = dz;
        }
        profs.mix_len = min(max(k, 1) * dz * 0.845, sgs_delta);
      }
    }

    virtual void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens, 
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy = 0)
    {if(timestep==0) surf_flux_sens = 0.;};

    virtual void update_surf_flux_lat(blitz::Array<real_t, n_dims>  surf_flux_lat, 
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy = 0)
    {if(timestep==0) surf_flux_lat = 0.;};

    virtual void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv,
                               blitz::Array<real_t, n_dims>  uv_ground,   
                               blitz::Array<real_t, n_dims>  U_ground,   
                               const real_t &U_ground_z,
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy = 0)
    {if(timestep==0) surf_flux_uv = 0.;};

    // ctor
    // TODO: these are DYCOMS definitions, move them there
    CasesCommon()
    {
      ForceParameters.heating_kappa = 85; // m^2/kg
      ForceParameters.F_0 = 70; // w/m^2
      ForceParameters.F_1 = 22; // w/m^2
      ForceParameters.q_i = 8e-3; // kg/kg
      ForceParameters.rho_i = 1.12; // kg/m^3
      ForceParameters.surf_latent_flux_in_watts_per_square_meter = true; // otherwise it's considered to be in [m/s]
      ForceParameters.surf_sensible_flux_in_watts_per_square_meter = true; // otherwise it's considered to be in [K m/s]
      ForceParameters.coriolis_parameter = 0.;
      ForceParameters.D = 0.; // large-scale wind horizontal divergence [1/s], needed in the radiation procedure of DYCOMS
      X = 0 * si::metres;
      Y = 0 * si::metres;
      Z = 0 * si::metres;
    }

    virtual ~CasesCommon() = default;

    protected:
  
    // function enforcing cyclic values in horizontal directions
    // 2D version
    template<class arr_t>
    void make_cyclic(arr_t arr,
      typename std::enable_if<arr_t::rank_ == 2>::type* = 0)
    { arr(arr.extent(0) - 1, blitz::Range::all()) = arr(0, blitz::Range::all()); }
  
    // 3D version
    template<class arr_t>
    void make_cyclic(arr_t arr,
      typename std::enable_if<arr_t::rank_ == 3>::type* = 0)
    { 
      arr(arr.extent(0) - 1, blitz::Range::all(), blitz::Range::all()) = 
        arr(0, blitz::Range::all(), blitz::Range::all()); 
      arr(blitz::Range::all(), arr.extent(1) - 1, blitz::Range::all()) = 
        arr(blitz::Range::all(), 0, blitz::Range::all());
    }
  };
};

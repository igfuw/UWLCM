#pragma once

#include <iostream>

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>


// simulation parameters container
// TODO: write them to rt_params directly in main()
struct user_params_t
{
  int nt, outfreq, spinup, rng_seed;
  setup::real_t dt, z_rlx_sclr;
  std::string outdir, model_case;
  bool th_src, rv_src, rc_src, rr_src, uv_src, w_src;
};


namespace setup 
{
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;

  const real_t D =  3.75e-6; // large-scale wind horizontal divergence [1/s]

  // container for constants that appear in forcings, some are not needed in all cases, etc...
  // TODO: make forcing functions part of case class
  struct ForceParameters_t
  {
    real_t q_i, heating_kappa, F_0, F_1, rho_i, D, u_fric;
    bool surf_latent_flux_in_watts_per_square_meter;
    bool surf_sensible_flux_in_watts_per_square_meter;
  };
 
  // CAUTION: new profiles have to be added to both structs and in copy_profiles below
  // TODO: try a different design where it is not necessary ?
  struct profiles_t
  {
    arr_1D_t th_e, p_e, rv_e, rl_e, th_ref, rhod, w_LS, hgt_fctr_sclr, hgt_fctr_vctr;

    profiles_t(int nz) :
    // rhod needs to be bigger, cause it divides vertical courant number
    // TODO: should have a halo both up and down, not only up like now; then it should be interpolated in courant calculation
      th_e(nz), p_e(nz), rv_e(nz), rl_e(nz), th_ref(nz), rhod(nz+1), w_LS(nz), hgt_fctr_vctr(nz), hgt_fctr_sclr(nz)
    {}
  };
  struct profile_ptrs_t
  {
    arr_1D_t *th_e, *p_e, *rv_e, *rl_e, *th_ref, *rhod, *w_LS, *hgt_fctr_sclr, *hgt_fctr_vctr;
  };
  // copy external profiles into rt_parameters
  // TODO: more elegant way
  template<class params_t>
  void copy_profiles(profiles_t &profs, params_t &p)
  {
    std::vector<std::pair<std::reference_wrapper<setup::arr_1D_t*>, std::reference_wrapper<setup::arr_1D_t>>> tobecopied = {
      {p.hgt_fctr_sclr, profs.hgt_fctr_sclr},
      {p.hgt_fctr_vctr, profs.hgt_fctr_vctr},
      {p.th_e         , profs.th_e         },
      {p.p_e          , profs.p_e          },
      {p.rv_e         , profs.rv_e         },
      {p.rl_e         , profs.rl_e         },
      {p.th_ref       , profs.th_ref       },
      {p.rhod         , profs.rhod         },
      {p.w_LS         , profs.w_LS         }
    };

    for (auto dst_src : tobecopied)
    {
      dst_src.first.get() = new setup::arr_1D_t(dst_src.second.get().dataFirst(), dst_src.second.get().shape(), blitz::neverDeleteData);
    }
  }

  template<class concurr_t>
  class CasesCommon
  {
    public:

    static constexpr int n_dims = concurr_t::solver_t::ct_params_t_::n_dims;

    using concurr_any_t = libmpdataxx::concurr::any<
      real_t, 
      n_dims
    >;

    ForceParameters_t ForceParameters;

    //th, rv and surface fluxes relaxation time and height
    const quantity<si::time, real_t> tau_rlx = 300 * si::seconds;

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
    real_t div_LS = 0.; // large-scale wind divergence (same as ForceParameters::D), 0. to turn off large-scale subsidence of SDs, TODO: add a process switch in libcloudph++ like for coal/cond/etc

    // hygroscopicity kappa of the aerosol 
    quantity<si::dimensionless, real_t> kappa = .61; // defaults to ammonium sulphate; CCN-derived value from Table 1 in Petters and Kreidenweis 2007

    virtual void setopts(typename concurr_t::solver_t::rt_params_t &params, const int nps[], const user_params_t &user_params) {assert(false);};
    virtual void intcond(concurr_any_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed) =0;
    virtual void env_prof(profiles_t &profs, int nz, const user_params_t &user_params) = 0;
    virtual void update_surf_flux_sens(typename concurr_t::solver_t::arr_sub_t &surf_flux_sens, int timestep, real_t dt) {if(timestep==0) surf_flux_sens = 0.;}; 
    virtual void update_surf_flux_lat(typename concurr_t::solver_t::arr_sub_t &surf_flux_lat, int timestep, real_t dt) {if(timestep==0) surf_flux_lat = 0.;};

    // ctor
    // TODO: these are DYCOMS definitions, move them there
    CasesCommon()
    {
      ForceParameters.heating_kappa = 85; // m^2/kg
      ForceParameters.F_0 = 70; // w/m^2
      ForceParameters.F_1 = 22; // w/m^2
      ForceParameters.q_i = 8e-3; // kg/kg
      ForceParameters.D = D; // large-scale wind horizontal divergence [1/s]
      ForceParameters.rho_i = 1.12; // kg/m^3
      ForceParameters.u_fric = 0.25; // m/s; friction velocity
      ForceParameters.surf_latent_flux_in_watts_per_square_meter = true; // otherwise it's considered to be in [m/s]
      ForceParameters.surf_sensible_flux_in_watts_per_square_meter = true; // otherwise it's considered to be in [K m/s]
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

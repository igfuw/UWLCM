#pragma once

#include <iostream>

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>

//#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

#include "../detail/user_params.hpp"
#include "../detail/concurr_types.hpp"
#include "../detail/ForceParameters.hpp"
#include "../detail/profiles.hpp"
#include "../detail/subs_t.hpp"

namespace cases
{
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;

  using real_t = setup::real_t;
  using arr_1D_t = setup::arr_1D_t;

  // helper, could be moved somewhere else
  struct hori_vel_t
  {
    bool initialized;
    real_t mean_vel; 
    const std::function<quantity<si::velocity, real_t>(real_t)> f_vel_prof;

    real_t operator()(const real_t &z) const
    {
      assert(initialized && "called uninitialized hori_vel_t");
      return f_vel_prof(z) * si::seconds / si::meters - mean_vel;
    }

    void init(bool window, quantity<si::length, real_t> Z) 
    {
      initialized=true;
      mean_vel = 0;

      if(window) // calculate mean of the velocity profile
      {
        for(int i=0; i < setup::mean_horvel_npts; ++i)
          mean_vel += f_vel_prof(i * (Z / si::meters) / (setup::mean_horvel_npts-1)) * si::seconds / si::meters;
        mean_vel /= setup::mean_horvel_npts;
      }
    }

    hori_vel_t(std::function<quantity<si::velocity, real_t>(real_t)> f): f_vel_prof(f), initialized(false) {}
  };

  template<class case_ct_params_t, int n_dims>
  class CasesCommon
  {
    static constexpr real_t karman_c = 0.41; // von Karman constant
    static constexpr real_t smg_c = 0.165; // Smagorinsky constant

    public: 

    using ix = typename case_ct_params_t::ix;
    using rt_params_t = typename case_ct_params_t::rt_params_t;

    using concurr_any_t = libmpdataxx::concurr::any<
      real_t, 
      n_dims
    >;

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

    // e-folding height for surface fluxes
    quantity<si::length, real_t> z_rlx = 0 * si::metres; // 0 indicates that there are no surf fluxes

    // hygroscopicity kappa of the aerosol 
    quantity<si::dimensionless, real_t> kappa = .61; // defaults to ammonium sulphate; CCN-derived value from Table 1 in Petters and Kreidenweis 2007

    real_t div_LS = 0.; // large-scale wind divergence (same as ForceParameters::D), 0. to turn off large-scale subsidence of SDs, TODO: add a process switch in libcloudph++ like for coal/cond/etc

    quantity<si::length, real_t> gccn_max_height; // GCCN added (at init and via relaxation) only up to this level

    template<bool enable_sgs = case_ct_params_t::enable_sgs>
    void setopts_sgs(rt_params_t &params,
                     typename std::enable_if<!enable_sgs>::type* = 0)
    {}

    template<bool enable_sgs = case_ct_params_t::enable_sgs>
    void setopts_sgs(rt_params_t &params,
                     typename std::enable_if<enable_sgs>::type* = 0)
    {
      params.c_m = 0.0856;
      params.smg_c = smg_c;
      params.prandtl_num = 0.42;
      params.karman_c = karman_c;
      params.cdrag = 0;
      params.fricvelsq = 0;
    }


    detail::ForceParameters_t ForceParameters;

    virtual void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params) {assert(false);};
    virtual void intcond(concurr_any_t &concurr, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed) =0;

    // virtual void set_profs(detail::profiles_t &profs, int nx, int nz, const user_params_t &user_params)
    virtual void set_profs(detail::profiles_t &profs, const int nps[n_dims], const user_params_t &user_params)
    {
      static_assert(n_dims == 2 || n_dims == 3, "set_profs: n_dims must be 2 or 3");
      // TODO: get these constants from user_params (or rt_params_t)
      // const real_t smg_c = 0.165; // Smagornsky constant
      // const real_t karman_c = 0.41; // von Karman constant

      const int nz = nps[n_dims-1];

      // grid spacing
      real_t dx, dy, dz;
      if constexpr (n_dims == 2)
      {
        dx = (X / si::metres) / (nps[0]-1);
        dy = dx; // used in horizontal mixing length...
        dz = (Z / si::metres) / (nps[1]-1);
      }
      else if (n_dims == 3)
      {
        dx = (X / si::metres) / (nps[0]-1);
        dy = (Y / si::metres) / (nps[1]-1);
        dz = (Z / si::metres) / (nps[2]-1);    
      }             

      // set SGS mixing length
      {
        using blitz::pow2;
        blitz::firstIndex k;

        real_t sgs_delta; // size of a grid cell used as characteristic length scale in SGS models
        if (user_params.sgs_delta > 0)
        {
          sgs_delta = user_params.sgs_delta;
        }
        else
        {
          sgs_delta = dz; // default
        }

        // size of SGS eddies, used in SD SGS models (note that they are isotropic!)
        profs.sgs_eddy_size = min(max(k, 0.5) * dz * 0.845, sgs_delta);

        // isotropic mixing length - used in isotropic SMG and in SD SGS models
        // like in Smagorinsky 1963, but with default characteristic length scale (sgs_delta) equal to dz (not dx*dy*dz)^(1/3))
        // similar approach used e.g. in the SAM model
        // also there's a reduction near ground 
        // NOTE: it's not multiplied by the smagorinsky constant here (that is done in slvr_sgs), because mix_len_iso is used also in SD SGS (opts_lgrngn.hpp) where no such constant is used
        // profs.mix_len_iso = min(max(k, 0.5) * dz * 0.845, sgs_delta);
        // profs.mix_len_iso_sq = pow2(smg_c * min(max(k, 0.5) * dz * 0.845, sgs_delta));
        profs.mix_len_iso_sq = pow2(smg_c * profs.sgs_eddy_size);

        // anisotropic mixing length - used in anisotropic SMG
        // based on Mason & Thomson 1992, but anisotropic, like in Simon & Chow 2021 (https://doi.org/10.1007/s10546-021-00642-0)
        // but we do not use grid size in the normal direction, only along (i.e. dx, dy for horizontal and dz for vertical)
        // TODO: to be consistent, setting user_params.sgs_delta or turning on SD SGS models should not be allowed with anisotropic SMG (?)
        // NOTE: it's not really mixing length, but grid size used to calculate mixing length; mixing length calculation is in slvr_sgs.hpp; 
        //       here we just set the characteristic grid size along each direction to simplify potential future option to override it (as with isotropic sgs_delta)
        real_t dxy = pow(dx * dy, real_t(1./2));
        profs.mix_len_hori_sq = pow2(real_t(1) / (real_t(1) / pow2(karman_c * max(k, 1) * dz) + real_t(1) / pow2(smg_c * dxy)));
        profs.mix_len_hori_sq = pow2(real_t(1) / (real_t(1) / pow2(karman_c * max(k, 1) * dz) + real_t(1) / pow2(smg_c * dz)));
        // real_t(1) / (real_t(1) / pow(this->smg_c * params.di, 2) + (*this->params.mix_len)(this->vert_idx))
        // profs.mix_len_vert_sq = dz; // real_t(1) / pow2(karman_c * max(k, 1) * dz);
        // add dx dy terms; also add sgs_delta, from user params...
        // TODO: save these, use them in slvr_sgs, ...
      }

      // set profile of vertical distribution of surface fluxes
      const real_t z_0 = z_rlx / si::metres;
      if(z_0 > 0.)
      {
        blitz::firstIndex k;
        // ------ version with flux added starting from 0th cell -----
        // (fraction of surface flux going out through upper cell boundary - fraction of surface flux going in through lower cell boundary) / cell height
        profs.hgt_fctr = (exp(- (k + 1 - 0.5) * dz / z_0) - exp(- (k - 0.5) * dz / z_0)) / dz;

        // uppermost and lowermost cells are lower
        profs.hgt_fctr(0) = (exp(- 0.5 * dz / z_0) - 1.) / (0.5 * dz);
        profs.hgt_fctr(nz-1) = (exp(- (nz - 1.) * dz / z_0) - exp(- (nz - 1.5) * dz / z_0) ) / (0.5 * dz);


        // ------ version with flux added starting from 1st cell -----
        /*
        // (fraction of surface flux going out through upper cell boundary - fraction of surface flux going in through lower cell boundary) / cell height
        profs.hgt_fctr = (exp(- (k + 1 - 1.0) * dz / z_0) - exp(- (k - 1.0) * dz / z_0)) / dz;

        // uppermost and lowermost cells are lower
        profs.hgt_fctr(0) = 0;
        profs.hgt_fctr(nz-1) = (exp(- (nz - 1.5) * dz / z_0) - exp(- (nz - 2.0) * dz / z_0) ) / (0.5 * dz);
*/
      }
    }

    virtual void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy = 0)
    {if(timestep==0) surf_flux_sens = 0.;};

    virtual void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy = 0)
    {if(timestep==0) surf_flux_lat = 0.;};

    virtual void update_surf_flux_uv(blitz::Array<real_t, n_dims> surf_flux_uv,
                                     blitz::Array<real_t, n_dims> uv_ground,   
                                     blitz::Array<real_t, n_dims> U_ground,   
                                     const real_t &U_ground_z,
                                     const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy = 0, const real_t &uv_mean = 0)
    {if(timestep==0) surf_flux_uv = 0.;};

    virtual void update_rv_LS(blitz::Array<real_t, 1> rv_LS,
                              int timestep, real_t dt, real_t dz)
    {};

    virtual void update_th_LS(blitz::Array<real_t, 1> th_LS,
                              int timestep, real_t dt, real_t dz)
    {};

    // ctor
    // TODO: these are DYCOMS definitions, move them there
    CasesCommon()
    {
      ForceParameters.heating_kappa = 85; // m^2/kg
      ForceParameters.F_0 = 70; // w/m^2
      ForceParameters.F_1 = 22; // w/m^2
      ForceParameters.q_i = 8e-3; // kg/kg
      ForceParameters.rho_i = 1.12; // kg/m^3
      ForceParameters.coriolis_parameter = 0.;
      ForceParameters.D = 0.; // large-scale wind horizontal divergence [1/s], needed in the radiation procedure of DYCOMS
      ForceParameters.uv_mean[0] = 0; // mean horizontal wind speed (for 'moving-window' simulations)
      ForceParameters.uv_mean[1] = 0; // mean horizontal wind speed (for 'moving-window' simulations)
      X = 0 * si::metres;
      Y = 0 * si::metres;
      Z = 0 * si::metres;
      gccn_max_height = 0 * si::metres;
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

#pragma once
#include "../cases/CasesCommon.hpp"
#include "slvr_dim.hpp"
#include <chrono>
#include <libmpdata++/git_revision.hpp>
#include <libcloudph++/git_revision.h>
#include "../../git_revision.h"

struct smg_tag  {};
struct iles_tag {};

template <class ct_params_t>
class slvr_common : public slvr_dim<ct_params_t>
{
  using parent_t = slvr_dim<ct_params_t>;

  public:
  using real_t = typename ct_params_t::real_t;
  using ix = typename ct_params_t::ix;
  using sgs_tag = typename std::conditional<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles, iles_tag, smg_tag>::type;

  protected:

  using clock = std::chrono::high_resolution_clock; // TODO: option to disable timing, as it may affect performance a little?
  // timing fields
  // TODO: timing slows down simulations
  //       either remove it and use profiling tools (e.g. vtune)
  //       or add some compile-time flag to turn it off
  clock::time_point tbeg, tend, tbeg_loop;
  std::chrono::milliseconds tdiag, tupdate, tsync, tsync_wait, tasync, tasync_wait, tloop, tvip_rhs, tnondelayed_step;

  int spinup; // number of timesteps

  // array with index of inversion
  blitz::Array<real_t, parent_t::n_dims-1> k_i;

  // array with sensible and latent heat surface flux
  blitz::Array<real_t, parent_t::n_dims> &surf_flux_sens;
  blitz::Array<real_t, parent_t::n_dims> &surf_flux_lat;
  // surface flux array filled with zeros ... TODO: add a way to set zero flux directly in libmpdata
  blitz::Array<real_t, parent_t::n_dims> &surf_flux_zero;

  // global arrays, shared among threads, TODO: in fact no need to share them?
  typename parent_t::arr_t &tmp1,
                           &r_l,
                           &F,       // forcings helper
                           &alpha,   // 'explicit' rhs part - does not depend on the value at n+1
                           &beta,    // 'implicit' rhs part - coefficient of the value at n+1
                           &radiative_flux,
                           &diss_rate; // TODO: move to slvr_sgs to save memory in iles simulations !;

  // surface precip stuff
  std::ofstream f_puddle; // output precipitation file
  
  // spinup stuff
  virtual bool get_rain() = 0;
  virtual void set_rain(bool) = 0;
  
  virtual void sgs_scalar_forces(const std::vector<int>&) {}
  virtual typename parent_t::arr_t get_rc(typename parent_t::arr_t&) = 0;

  // get shape from a rng_t or an idx_t
  inline int shape(const rng_t &rng) { return rng.length();}
  template<int n_dims>
  blitz::TinyVector<int, n_dims> shape(const idx_t<n_dims> &rng) { return rng.ubound() - rng.lbound() + 1;}

  void buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv);
  void radiation(typename parent_t::arr_t &rv);
  void rv_src();
  void th_src(typename parent_t::arr_t &rv);
  void w_src(typename parent_t::arr_t &th, typename parent_t::arr_t &rv, const int at);

  void surf_sens_impl(smg_tag);
  void surf_sens_impl(iles_tag);

  void surf_latent_impl(smg_tag);
  void surf_latent_impl(iles_tag);

  void surf_sens();
  void surf_latent();

  void subsidence(const int&);
  void coriolis(const int&);

  void hook_ante_loop(int nt); 
  void hook_ante_step(); 
  void hook_post_step(); 
  void vip_rhs_expl_calc(); 
  virtual void diag(); 
  void record_all(); 
  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  );

  public:

  // note dual inheritance to get profile pointers
  struct rt_params_t : parent_t::rt_params_t, setup::profile_ptrs_t
  { 
    int spinup = 0, // number of timesteps during which autoconversion is to be turned off
        nt;         // total number of timesteps
    bool rv_src, th_src, uv_src, w_src, subsidence, coriolis, friction, buoyancy_wet, radiation;
    bool rc_src, rr_src; // these two are only relevant for blk_1m, but need to be here so that Cases can have access to it
    typename ct_params_t::real_t dz; // vertical grid size
    setup::ForceParameters_t ForceParameters;
    user_params_t user_params; // copy od user_params needed only for output to const.h5, since the output has to be done at the end of hook_ante_loop

    // functions for updating surface fluxes per timestep
    std::function<void(typename parent_t::arr_t, int, const real_t&, const real_t&, const real_t&)> update_surf_flux_sens;
    std::function<void(typename parent_t::arr_t, int, const real_t&, const real_t&, const real_t&)> update_surf_flux_lat;
  };

  // per-thread copy of params
  rt_params_t params;

  // ctor
  slvr_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p),
    spinup(p.spinup),
    tmp1(args.mem->tmp[__FILE__][0][0]),
    r_l(args.mem->tmp[__FILE__][0][1]),
    F(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    radiative_flux(args.mem->tmp[__FILE__][0][5]),
    diss_rate(args.mem->tmp[__FILE__][0][6]),
    surf_flux_sens(args.mem->tmp[__FILE__][1][0]),
    surf_flux_lat(args.mem->tmp[__FILE__][1][1]),
    surf_flux_zero(args.mem->tmp[__FILE__][1][2])
  {
    k_i.resize(this->shape(this->hrzntl_domain)); // TODO: resize to hrzntl_subdomain
    r_l = 0.;
    surf_flux_zero = 0.;
  }

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 7); // tmp1, tmp2, r_l, alpha, beta, F, diss_rate, radiative_flux
    parent_t::alloc_tmp_sclr(mem, __FILE__, 3, "", true); // surf_flux_sens, surf_flux_lat, surf_flux_zero
  }
};

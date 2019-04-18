#pragma once
#include "cases/CasesCommon.hpp"
#include "slvr_dim.hpp"
#include <chrono>
#include <libmpdata++/git_revision.hpp>
#include <libcloudph++/git_revision.h>
#include "../git_revision.h"

template <class ct_params_t>
class slvr_common : public slvr_dim<ct_params_t>
{
  using parent_t = slvr_dim<ct_params_t>;

  public:
  using real_t = typename ct_params_t::real_t;
  using ix = typename ct_params_t::ix;

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
  blitz::Array<real_t, parent_t::n_dims-1> surf_flux_sens;
  blitz::Array<real_t, parent_t::n_dims-1> surf_flux_lat;

  // global arrays, shared among threads, TODO: in fact no need to share them?
  typename parent_t::arr_t &tmp1,
                           &r_l,
                           &F,       // forcings helper
                           &alpha,   // 'explicit' rhs part - does not depend on the value at n+1
                           &beta;    // 'implicit' rhs part - coefficient of the value at n+1

  // surface precip stuff
  std::ofstream f_puddle; // output precipitation file
  
  // spinup stuff
  virtual bool get_rain() = 0;
  virtual void set_rain(bool) = 0;

  void hook_ante_loop(int nt) 
  {
    if (spinup > 0)
    {
      set_rain(false);
    }
    else
      set_rain(true);

    parent_t::hook_ante_loop(nt); 

    // open file for output of precitpitation volume
    if(this->rank==0)
      f_puddle.open(this->outdir+"/prec_vol.dat");

    // record user_params and profiles
    if(this->rank==0)
    {
      this->record_aux_const(std::string("UWLCM git_revision : ") + UWLCM_GIT_REVISION, -44);  
#ifdef LIBMPDATAXX_GIT_REVISION
      this->record_aux_const(std::string("LIBMPDATAXX git_revision : ") + LIBMPDATAXX_GIT_REVISION, -44);  
#else
      throw std::runtime_error("LIBMPDATAXX_GIT_REVISION is not defined, update your libmpdata++ library");
#endif
#ifdef LIBCLOUDPHXX_GIT_REVISION
      this->record_aux_const(std::string("LIBCLOUDPHXX git_revision : ") + LIBCLOUDPHXX_GIT_REVISION, -44);  
#else
      throw std::runtime_error("LIBCLOUDPHXX_GIT_REVISION is not defined, update your libcloudph++ library");
#endif
      this->record_aux_const("omp_max_threads (on MPI rank 0)", omp_get_max_threads());  
      this->record_aux_const(std::string("user_params case : ") + params.user_params.model_case, -44);  
      this->record_aux_const("user_params nt", params.user_params.nt);  
      this->record_aux_const("user_params dt", params.user_params.dt);  
      this->record_aux_const("user_params outfreq", params.user_params.outfreq);  
      this->record_aux_const(std::string("user_params outdir : ") +  params.user_params.outdir, -44);  
      this->record_aux_const("user_params spinup", params.user_params.spinup);  
      this->record_aux_const("user_params rng_seed", params.user_params.rng_seed);  
      this->record_aux_const("user_params z_rlx_sclr", params.user_params.z_rlx_sclr);  
      this->record_aux_const("user_params th_src", params.user_params.th_src);  
      this->record_aux_const("user_params rv_src", params.user_params.rv_src);  
      this->record_aux_const("user_params uv_src", params.user_params.uv_src);  
      this->record_aux_const("user_params w_src", params.user_params.w_src);  

      this->record_aux_const("rt_params th_src", params.th_src);  
      this->record_aux_const("rt_params rv_src", params.rv_src);  
      this->record_aux_const("rt_params uv_src", params.uv_src);  
      this->record_aux_const("rt_params w_src", params.w_src);  
      this->record_aux_const("rt_params spinup", params.spinup);  
      this->record_aux_const("rt_params subsidence", params.subsidence);  
      this->record_aux_const("rt_params coriolis", params.coriolis);  
      this->record_aux_const("rt_params friction", params.friction);  
      this->record_aux_const("rt_params buoyancy_wet", params.buoyancy_wet);  

      this->record_aux_const("ForceParameters q_i", params.ForceParameters.q_i);  
      this->record_aux_const("ForceParameters heating_kappa", params.ForceParameters.heating_kappa);  
      this->record_aux_const("ForceParameters F_0", params.ForceParameters.F_0);  
      this->record_aux_const("ForceParameters F_1", params.ForceParameters.F_1);  
      this->record_aux_const("ForceParameters rho_i", params.ForceParameters.rho_i);  
      this->record_aux_const("ForceParameters D", params.ForceParameters.D);  
      this->record_aux_const("ForceParameters u_fric", params.ForceParameters.u_fric);  
      this->record_aux_const("ForceParameters coriolis_parameter", params.ForceParameters.coriolis_parameter);  
      this->record_aux_const("ForceParameters surf_latent_flux_in_watts_per_square_meter", params.ForceParameters.surf_latent_flux_in_watts_per_square_meter);  
      this->record_aux_const("ForceParameters surf_sensible_flux_in_watts_per_square_meter", params.ForceParameters.surf_sensible_flux_in_watts_per_square_meter);  
     
      // recording profiles
      this->record_prof_const("th_e", params.th_e->data()); 
      this->record_prof_const("p_e", params.p_e->data()); 
      this->record_prof_const("rv_e", params.rv_e->data()); 
      this->record_prof_const("th_ref", params.th_ref->data()); 
      this->record_prof_const("rhod", params.rhod->data()); 
      this->record_prof_const("u_geostr", params.geostr[0]->data()); 
      this->record_prof_const("v_geostr", params.geostr[1]->data()); 
    }
 
    // initialize surf fluxes with timestep==0
    params.update_surf_flux_sens(surf_flux_sens, 0, this->dt);
    params.update_surf_flux_lat(surf_flux_lat, 0, this->dt);
  }

  void hook_ante_step()
  {
    if (spinup != 0 && spinup == this->timestep)
    {
      // turn autoconversion on only after spinup (if spinup was specified)
      set_rain(true);
    }
    parent_t::hook_ante_step(); 
  }


  // get shape from a rng_t or an idx_t
  inline int shape(const rng_t &rng) { return rng.length();}
  template<int n_dims>
  blitz::TinyVector<int, n_dims> shape(const idx_t<n_dims> &rng) { return rng.ubound() - rng.lbound() + 1;}

  void buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv);
  void radiation(typename parent_t::arr_t &rv);
  void rv_src();
  void th_src(typename parent_t::arr_t &rv);
  void w_src(typename parent_t::arr_t &th, typename parent_t::arr_t &rv, const int at);
  void surf_sens();
  void surf_latent();
  void subsidence(const int&);
  void coriolis(const int&);

  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  )
  {
    parent_t::update_rhs(rhs, dt, at); // zero-out rhs
    this->mem->barrier();
    if(this->rank == 0)
      tbeg = clock::now();

    using ix = typename ct_params_t::ix;

    const auto &ijk = this->ijk;

    // forcing
    switch (at)
    {
      // for eulerian integration or used to init trapezoidal integration
      case (0):
      {
        // ---- water vapor sources ----
        rv_src();
        rhs.at(ix::rv)(ijk) += alpha(ijk) + beta(ijk) * this->state(ix::rv)(ijk);

        // ---- potential temp sources ----
        th_src(this->state(ix::rv));
        rhs.at(ix::th)(ijk) += alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk);

        // vertical velocity sources
        if(params.w_src && (!ct_params_t::piggy))
        {
          w_src(this->state(ix::th), this->state(ix::rv), 0);
          rhs.at(ix::w)(ijk) += alpha(ijk);
        }

        // horizontal velocity sources 
        if(params.uv_src)
        {
          for(auto type : this->hori_vel)
          {
            // subsidence
            subsidence(type);
            rhs.at(type)(ijk) += F(ijk);

            // Coriolis
            coriolis((type+1) % this->hori_vel.size());
            if(type == ix::u)
              rhs.at(type)(ijk) += F(ijk);
            else
              rhs.at(type)(ijk) -= F(ijk);
          }
        }
        
        break;
      }
      case (1):
      // trapezoidal rhs^n+1
      {
        // ---- water vapor sources ----
//        rv_src();
//        rhs.at(ix::rv)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::rv)(ijk)) / (1. - 0.5 * this->dt * beta(ijk));
        // TODO: alpha should also take (possibly impolicit) estimate of rv^n+1 too
        //       becomes important when nudging is introduced?


        // ---- potential temp sources ----
//        tmp2(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
        // todo: once rv_src beta!=0 (e.g. nudging), rv^n+1 estimate should be implicit here
//        th_src(tmp2);
//        rhs.at(ix::th)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk)) / (1. - 0.5 * this->dt * beta(ijk));
        // TODO: alpha should also take (possibly impolicit) estimate of th^n+1 too
        //       becomes important when nudging is introduced?

        // vertical velocity sources
        if(params.w_src && (!ct_params_t::piggy))
        {
          // temporarily use beta to store the th^n+1 estimate
//          beta(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
          // todo: oncethv_src beta!=0 (e.g. nudging), th^n+1 estimate should be implicit here

          // temporarily use F to store the rv^n+1 estimate
  //        F(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
          // todo: once rv_src beta!=0 (e.g. nudging), rv^n+1 estimate should be implicit here

    //      w_src(beta, F, 1);
          w_src(this->state(ix::th), this->state(ix::rv), 1);
          rhs.at(ix::w)(ijk) += alpha(ijk);
        }

        // horizontal velocity sources 
        // large-scale vertical wind
/*
        if(params.uv_src)
        {
          for(auto type : this->hori_vel)
          {
            subsidence(type); 
            rhs.at(type)(ijk) += F(ijk);
          }
        }
*/
        break;
      }
    }
    nancheck(rhs.at(ix::th)(this->ijk), "RHS of th after rhs_update");
    nancheck(rhs.at(ix::rv)(this->ijk), "RHS of rv after rhs_update");
    nancheck(rhs.at(ix::w)(this->ijk), "RHS of w after rhs_update");
    for(auto type : this->hori_vel)
      {nancheck(rhs.at(type)(this->ijk), (std::string("RHS of horizontal velocity after rhs_update, type: ") + std::to_string(type)).c_str());}
    this->mem->barrier();
    if(this->rank == 0)
    {
      tend = clock::now();
      tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
    }
  }


  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();

    if(!params.friction) return;
  
    this->mem->barrier();
    if(this->rank == 0)
      tbeg = clock::now();
    // kinematic momentum flux  = -u_fric^2 * u_i / |U| * exponential decay
    typename parent_t::arr_sub_t U_ground(this->shape(this->hrzntl_subdomain));
    U_ground = this->calc_U_ground();

    // loop over horizontal dimensions
    for(int it = 0; it < parent_t::n_dims-1; ++it)
    {
      this->vip_rhs[it](this->ijk).reindex(this->zero) += 
        where(U_ground(blitz::tensor::i, blitz::tensor::j) == 0., 0., 
          -2 * pow(params.ForceParameters.u_fric,2) *  // 2, because it is multiplied by 0.5 in vip_rhs_apply
          this->vip_ground[it](blitz::tensor::i, blitz::tensor::j) /              // u_i at z=0
          U_ground(blitz::tensor::i, blitz::tensor::j) *  // |U| at z=0
          (*params.hgt_fctr_vctr)(this->vert_idx)                                       // hgt_fctr 
        );
    }

    for(int it = 0; it < parent_t::n_dims-1; ++it)
      {nancheck(this->vip_rhs[it](this->ijk), (std::string("vip_rhs after vip_rhs_expl_calc type: ") + std::to_string(it)).c_str());} 
    this->mem->barrier();
    if(this->rank == 0)
    {
      tend = clock::now();
      tvip_rhs += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
    }
  }

  void hook_post_step()
  {
    negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv at start of slvr_common::hook_post_step");
    parent_t::hook_post_step(); // includes output
    this->mem->barrier();

    if (this->rank == 0) 
    {
      // there's no hook_post_loop, so we imitate it here to write out computation times, TODO: move to destructor?
      if(this->timestep == params.nt) // timestep incremented before post_step
      {
        tend = clock::now();
        tloop = std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg_loop );
        std::cout <<  "wall time in milliseconds: " << std::endl
          << "loop: " << tloop.count() << std::endl
          << "custom rhs update: " << tupdate.count() << " ("<< setup::real_t(tupdate.count())/tloop.count()*100 <<"%)" << std::endl
          << "custom vip_rhs: " << tvip_rhs.count() << " ("<< setup::real_t(tvip_rhs.count())/tloop.count()*100 <<"%)" << std::endl
          << "diag: " << tdiag.count() << " ("<< setup::real_t(tdiag.count())/tloop.count()*100 <<"%)" << std::endl
          << "sync: " << tsync.count() << " ("<< setup::real_t(tsync.count())/tloop.count()*100 <<"%)" << std::endl
          << "nondelayed step: " << tnondelayed_step.count() << " ("<< setup::real_t(tnondelayed_step.count())/tloop.count()*100 <<"%)" << std::endl
          << "async: " << tasync.count() << " ("<< setup::real_t(tasync.count())/tloop.count()*100 <<"%)" << std::endl
          << "async_wait: " << tasync_wait.count() << " ("<< setup::real_t(tasync_wait.count())/tloop.count()*100 <<"%)" << std::endl
          << "sync_wait: " << tsync_wait.count() << " ("<< setup::real_t(tsync_wait.count())/tloop.count()*100 <<"%)" << std::endl;
      }
    }
    negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv at end of slvr_common::hook_post_step");
  }

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
    std::function<void(typename parent_t::arr_sub_t&, int, real_t)> update_surf_flux_sens;
    std::function<void(typename parent_t::arr_sub_t&, int, real_t)> update_surf_flux_lat;
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
    r_l(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    F(args.mem->tmp[__FILE__][0][1])
  {
    k_i.resize(this->shape(this->hrzntl_domain)); // TODO: resize to hrzntl_subdomain
    surf_flux_sens.resize(this->shape(this->hrzntl_domain)); // TODO: resize to hrzntl_subdomain
    surf_flux_lat.resize(this->shape(this->hrzntl_domain)); // TODO: resize to hrzntl_subdomain
    r_l = 0.;
  }

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 5);
  }
};

#pragma once
#include "../cases/CasesCommon.hpp"
#include "slvr_dim.hpp"
#include <chrono>
#include <libmpdata++/git_revision.hpp>
#include <libcloudph++/git_revision.h>
#include <libcloudph++/common/output.hpp>
#include "../detail/get_uwlcm_git_revision.hpp"
#include "../detail/ForceParameters.hpp"

struct smg_tag  {};
struct iles_tag {};

namespace cmn = libcloudphxx::common;

template <class ct_params_t>
class slvr_common : public slvr_dim<ct_params_t>
{
  using parent_t = slvr_dim<ct_params_t>;

  public:
  using real_t = typename ct_params_t::real_t;
  using ix = typename ct_params_t::ix;
  using sgs_tag = typename std::conditional<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles, iles_tag, smg_tag>::type;

  protected:

#if defined(UWLCM_TIMING)
  using clock = std::chrono::system_clock; 
  using timer = std::chrono::milliseconds; 
  timer tsync, tsync_gpu, tsync_wait, tasync, tasync_gpu, tasync_wait, tasync_wait_in_record_all; // timings used in lgrngn solver TODO: move them to slvr_lgrngn
#else
  using timer = void; 
  using clock = void; 
#endif

  static constexpr int n_flxs = ct_params_t::n_dims + 1; // number of surface fluxes = number of hori velocities + th + rv

  // array with index of inversion
  blitz::Array<real_t, parent_t::n_dims-1> k_i; // TODO: allocate k_i with alloc surf + in MPI calc average k_i over all processes

/*
  TODO: an array (map?) of surf fluxes, something like:
  // array with sensible and latent heat surface flux
  // one of them is a surface flux array filled with zeros ... TODO: add a way to set zero flux directly in libmpdata
  std::array<blitz::Array<real_t, parent_t::n_dims>, n_flxs> surf_fluxes;
  */

  blitz::Array<real_t, parent_t::n_dims> surf_flux_sens,
                                         surf_flux_lat,
                                         surf_flux_u,
                                         surf_flux_v,
                                         surf_flux_tmp,
                                         surf_flux_zero, // zero-filled array, find a way to avoid this
                                         U_ground;

  // global arrays, shared among threads, TODO: in fact no need to share them?
  typename parent_t::arr_t &tmp1,
                           &r_l,
                           &F,       // forcings helper
                           &alpha,   // 'explicit' rhs part - does not depend on the value at n+1
                           &beta,    // 'implicit' rhs part - coefficient of the value at n+1
                           &radiative_flux,
                           &diss_rate; // TODO: move to slvr_sgs to save memory in iles simulations !;

  setup::arr_1D_t th_mean_prof, rv_mean_prof; // profiles with mean th/rv at each level
                                              // NOTE: these profiles hold the same values for all threads (and MPI processes),
                                              //       but each thread has it's own copy. We could have one per MPI process to save memory

  // precip output
  std::map<cmn::output_t, real_t> puddle;
  const int n_puddle_scalars = cmn::output_names.size();
  virtual void get_puddle() = 0;

  // spinup stuff
  virtual bool get_rain() = 0;
  virtual void set_rain(bool) = 0;
  
  virtual void sgs_scalar_forces(const std::vector<int>&) {}
  virtual typename parent_t::arr_t get_rc(typename parent_t::arr_t&) = 0;

  //void common_water_src(int, int);

  void hook_ante_loop(int nt)
  {
    if (params.user_params.spinup > 0)
    {
      set_rain(false);
    }
    else
      set_rain(true);

    parent_t::hook_ante_loop(nt);

    // record user_params and profiles
    if(this->rank==0)
    {
      this->record_aux_const("UWLCM git_revision", "git_revisions", get_uwlcm_git_revision());  
#ifdef LIBMPDATAXX_GIT_REVISION
      this->record_aux_const("LIBMPDATAXX git_revision", "git_revisions", LIBMPDATAXX_GIT_REVISION);  
#else
      static_assert(false, "LIBMPDATAXX_GIT_REVISION is not defined, update your libmpdata++ library");
#endif
#ifdef LIBCLOUDPHXX_GIT_REVISION
      this->record_aux_const("LIBCLOUDPHXX git_revision", "git_revisions", LIBCLOUDPHXX_GIT_REVISION);  
#else
      static_assert(false, "LIBCLOUDPHXX_GIT_REVISION is not defined, update your libcloudph++ library");
#endif
      std::string run_command;
      for(int i=0; i<ac; ++i)
        run_command += std::string(av[i]) + std::string(" ");
      this->record_aux_const("run command", "misc", run_command);  
      this->record_aux_const("omp_max_threads (on MPI rank 0)", omp_get_max_threads());  
      this->record_aux_const("MPI compiler (true/false)", "MPI details", 
#ifdef USE_MPI
        1
#else
        0
#endif
        );
      this->record_aux_const("MPI size", "MPI details", this->mem->distmem.size());

      this->record_aux_const("case", "user_params", params.user_params.model_case);
      this->record_aux_const("nt", "user_params", params.user_params.nt);
      this->record_aux_const("dt", "user_params", params.user_params.dt);  
      this->record_aux_const("X", "user_params", params.user_params.X);  
      this->record_aux_const("Y", "user_params", params.user_params.Y);  
      this->record_aux_const("Z", "user_params", params.user_params.Z);  
      this->record_aux_const("outfreq", "user_params", params.user_params.outfreq);  
      this->record_aux_const("outdir", "user_params", params.user_params.outdir);
      this->record_aux_const("spinup", "user_params", params.user_params.spinup);  
      this->record_aux_const("rng_seed", "user_params", params.user_params.rng_seed);  
      this->record_aux_const("rng_seed_init", "user_params", params.user_params.rng_seed_init);  
      this->record_aux_const("sgs_delta", "user_params", params.user_params.sgs_delta);  
      this->record_aux_const("relax_th_rv", "user_params", params.user_params.relax_th_rv);  
      this->record_aux_const("case_n_stp_multiplier", "user_params", params.user_params.case_n_stp_multiplier);  
      this->record_aux_const("n_fra_iter", "user_params", params.user_params.n_fra_iter);  

      this->record_aux_const("th_src", "rt_params", params.th_src);  
      this->record_aux_const("rv_src", "rt_params", params.rv_src);  
      this->record_aux_const("uv_src", "rt_params", params.uv_src);  
      this->record_aux_const("w_src",  "rt_params", params.w_src);  
      this->record_aux_const("outfreq", "rt_params", params.outfreq);  
      this->record_aux_const("outdir", "rt_params", params.outdir);

      this->record_aux_const("subsidence", "rt_params", params.subsidence);  
      this->record_aux_const("vel_subsidence", "rt_params", params.vel_subsidence);  
      this->record_aux_const("coriolis", "rt_params", params.coriolis);  
      this->record_aux_const("friction", "rt_params", params.friction);  
      this->record_aux_const("buoyancy_wet", "rt_params", params.buoyancy_wet);  
      this->record_aux_const("radiation", "rt_params", params.radiation);  
      this->record_aux_const("dz", "rt_params", params.dz);  

      this->record_aux_const("q_i", "ForceParameters", params.ForceParameters.q_i);  
      this->record_aux_const("heating_kappa", "ForceParameters", params.ForceParameters.heating_kappa);  
      this->record_aux_const("F_0", "ForceParameters", params.ForceParameters.F_0);  
      this->record_aux_const("F_1", "ForceParameters", params.ForceParameters.F_1);  
      this->record_aux_const("rho_i", "ForceParameters", params.ForceParameters.rho_i);  
      this->record_aux_const("D", "ForceParameters", params.ForceParameters.D);  
      this->record_aux_const("coriolis_parameter", "ForceParameters", params.ForceParameters.coriolis_parameter);  

      // TODO: need to update ref files in tests to include this in output
      this->record_aux_const("mean_rd1", "user_params", params.user_params.mean_rd1 / si::metres);  
      this->record_aux_const("sdev_rd1", "user_params", params.user_params.sdev_rd1);
      this->record_aux_const("n1_stp", "user_params", params.user_params.n1_stp * si::cubic_metres);
      this->record_aux_const("kappa1", "user_params", params.user_params.kappa1);
      this->record_aux_const("mean_rd2", "user_params", params.user_params.mean_rd2 / si::metres);  
      this->record_aux_const("sdev_rd2", "user_params", params.user_params.sdev_rd2);
      this->record_aux_const("n2_stp", "user_params", params.user_params.n2_stp * si::cubic_metres);
      this->record_aux_const("kappa2", "user_params", params.user_params.kappa2);

      // recording profiles
      this->record_prof_const("th_e", params.profs.th_e.data()); 
      this->record_prof_const("p_e", params.profs.p_e.data()); 
      this->record_prof_const("rv_e", params.profs.rv_e.data()); 
      this->record_prof_const("rl_e", params.profs.rl_e.data()); 
      this->record_prof_const("th_reference", params.profs.th_reference.data()); 
      this->record_prof_const("rhod", params.profs.rhod.data()); 
      this->record_prof_const("w_LS", params.profs.w_LS.data()); 
      this->record_prof_const("th_LS", params.profs.th_LS.data()); 
      this->record_prof_const("rv_LS", params.profs.rv_LS.data()); 
      this->record_prof_const("hgt_fctr", params.profs.hgt_fctr.data()); 
      this->record_prof_const("mix_len", params.profs.mix_len.data());
      this->record_prof_const("relax_th_rv_coeff", params.profs.relax_th_rv_coeff.data()); 
      if(parent_t::n_dims==3)
      {
        this->record_prof_const("u_geostr", params.profs.geostr[0].data()); 
        this->record_prof_const("v_geostr", params.profs.geostr[1].data()); 
      }

      this->record_prof_const("refined th_e", params.profs_ref.th_e.data(), true); 
      this->record_prof_const("refined p_e", params.profs_ref.p_e.data(), true); 
      this->record_prof_const("refined rv_e", params.profs_ref.rv_e.data(), true); 
      this->record_prof_const("refined rl_e", params.profs_ref.rl_e.data(), true); 
      this->record_prof_const("refined th_reference", params.profs_ref.th_reference.data(), true); 
      this->record_prof_const("refined rhod", params.profs_ref.rhod.data(), true); 
      this->record_prof_const("refined w_LS", params.profs_ref.w_LS.data(), true); 
      this->record_prof_const("refined th_LS", params.profs_ref.th_LS.data(), true); 
      this->record_prof_const("refined rv_LS", params.profs_ref.rv_LS.data(), true); 
      this->record_prof_const("refined hgt_fctr", params.profs_ref.hgt_fctr.data(), true); 
      this->record_prof_const("refined mix_len", params.profs_ref.mix_len.data(), true);
      this->record_prof_const("refined relax_th_rv_coeff", params.profs_ref.relax_th_rv_coeff.data(), true); 
      if(parent_t::n_dims==3)
      {
        this->record_prof_const("refined u_geostr", params.profs_ref.geostr[0].data(), true); 
        this->record_prof_const("refined v_geostr", params.profs_ref.geostr[1].data(), true); 
      }
    }
    
    this->mem->barrier();
    
    // initialize surf fluxes with timestep==0
    U_ground(this->hrzntl_slice(0)) = this->calc_U_ground();

    params.update_surf_flux_sens(
      surf_flux_sens(this->hrzntl_slice(0)).reindex(this->origin),
      this->state(ix::th)(this->hrzntl_slice(0)).reindex(this->origin),
      U_ground(this->hrzntl_slice(0)).reindex(this->origin),
      params.dz / 2, 0, this->dt, this->di, this->dj
    );
    params.update_surf_flux_lat(
      surf_flux_lat(this->hrzntl_slice(0)).reindex(this->origin),
      this->state(ix::rv)(this->hrzntl_slice(0)).reindex(this->origin),
      U_ground(this->hrzntl_slice(0)).reindex(this->origin), 
      params.dz / 2, 0, this->dt, this->di, this->dj
    );
    params.update_surf_flux_uv(
      surf_flux_u(this->hrzntl_slice(0)).reindex(this->origin),
      this->state(ix::vip_i)(this->hrzntl_slice(0)).reindex(this->origin),
      U_ground(this->hrzntl_slice(0)).reindex(this->origin), 
      params.dz / 2, 0, this->dt, this->di, this->dj
    );
    if(parent_t::n_dims==3)
    {
      params.update_surf_flux_uv(
        surf_flux_v(this->hrzntl_slice(0)).reindex(this->origin),
        this->state(ix::vip_j)(this->hrzntl_slice(0)).reindex(this->origin),
        U_ground(this->hrzntl_slice(0)).reindex(this->origin),
        params.dz / 2, 0, this->dt, this->di, this->dj
      );
    }
  }

  void hook_ante_step()
  {
    if (params.user_params.spinup != 0 && params.user_params.spinup == this->timestep)
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

  // get base from a rng_t or an idx_t
  inline blitz::TinyVector<int, 1> base(const rng_t &rng) { return blitz::TinyVector<int, 1>(rng.first());}
  template<int n_dims>
  blitz::TinyVector<int, n_dims> base(const idx_t<n_dims> &rng) { return rng.lbound();}

  void buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv);
  void radiation(typename parent_t::arr_t &rv);
  void rv_src();
  void th_src(typename parent_t::arr_t &rv);
  void w_src(typename parent_t::arr_t &th, typename parent_t::arr_t &rv, const int at);

  void surf_sens_impl(smg_tag);
  void surf_sens_impl(iles_tag);

  void surf_latent_impl(smg_tag);
  void surf_latent_impl(iles_tag);

  void surf_u_impl(smg_tag);
  void surf_u_impl(iles_tag);

  void surf_v_impl(smg_tag);
  void surf_v_impl(iles_tag);

  void surf_sens();
  void surf_latent();
  void surf_u();
  void surf_v();

  void subsidence(const int&);
  void coriolis(const int&);
  void relax_th_rv(const int&);

  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  )
  {
    parent_t::update_rhs(rhs, dt, at); // zero-out rhs
    this->mem->barrier();

    using ix = typename ct_params_t::ix;

    const auto &ijk = this->ijk;

    // forcing
    switch (at)
    {
      // for eulerian integration or used to init trapezoidal integration
      case (0):
      {
        // calculate surface wind magnitude, TODO: not needed if there are no surface fluxes
        U_ground(this->hrzntl_slice(0)) = this->calc_U_ground();

        // calculate mean th and rv at each level, needed for nudging of the horizontal mean
        if(params.user_params.relax_th_rv)
        {
          this->hrzntl_mean(this->state(ix::rv), rv_mean_prof);
          this->hrzntl_mean(this->state(ix::th), th_mean_prof);
        }

        // ---- water vapor sources ----
        rv_src();
        rhs.at(ix::rv)(ijk) += alpha(ijk);// + beta(ijk) * this->state(ix::rv)(ijk);

        // ---- potential temp sources ----
        th_src(this->state(ix::rv));
        rhs.at(ix::th)(ijk) += alpha(ijk);// + beta(ijk) * this->state(ix::th)(ijk);

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
            // TODO: move these to uv_src
            // subsidence
            if(params.vel_subsidence)
            {
              subsidence(type);
              rhs.at(type)(ijk) += F(ijk);
            }

            // Coriolis
            coriolis((type+1) % this->hori_vel.size());
            if(type == ix::vip_i)
              rhs.at(type)(ijk) += F(ijk);
            else if(type == ix::vip_j)
              rhs.at(type)(ijk) -= F(ijk);
          }

          // surface flux with exp folding in vertical (only in ILES)
          if(params.friction && ct_params_t::sgs_scheme == libmpdataxx::solvers::iles)
          {
            for(auto type : this->hori_vel)
            {
              if(type == ix::vip_i)
                surf_u();
              else if(type == ix::vip_j)
                surf_v();
              rhs.at(type)(ijk) += F(ijk);
            }
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
  }


/*
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
          (params.profs.hgt_fctr_vctr)(this->vert_idx)                                       // hgt_fctr
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
*/

  void hook_post_step()
  {
    negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv at start of slvr_common::hook_post_step");
    parent_t::hook_post_step(); // includes output
    this->mem->barrier();
    negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv at end of slvr_common::hook_post_step");
  }

  void hook_mixed_rhs_ante_step()
  {
    this->update_rhs(this->rhs, this->dt, 0); // TODO: update_rhs called twice per step causes halo filling twice (done by parent_t::update_rhs), probably not needed - we just need to set rhs to zero (?)
    this->apply_rhs(this->dt);

    // rv might be negative due to large negative RHS from SD fluctuations + large-scale subsidence?
    // turn all negative rv into rv = 0... CHEATING
    negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv after mixed_rhs_ante_step apply rhs");
    nancheck(this->mem->advectee(ix::th)(this->ijk), "th after mixed_rhs_ante_step apply rhs");
    nancheck(this->mem->advectee(ix::rv)(this->ijk), "rv after mixed_rhs_ante_step apply rhs");
    negcheck(this->mem->advectee(ix::th)(this->ijk), "th after mixed_rhs_ante_step apply rhs");
    negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after mixed_rhs_ante_step apply rhs");
  }

  void hook_mixed_rhs_post_step()
  {
    negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after advection");

    this->update_rhs(this->rhs, this->dt, 1);
    negcheck(this->rhs.at(ix::rv)(this->ijk), "RHS rv after update_rhs in mixed_rhs_post_step");
    this->apply_rhs(this->dt);

    // no negtozero after apply, because only w changed here
    // TODO: add these nanchecks/negchecks to apply_rhs, since they are repeated twice now
    nancheck(this->mem->advectee(ix::th)(this->ijk), "th after mixed_rhs_post_step apply rhs");
    nancheck(this->mem->advectee(ix::rv)(this->ijk), "rv after mixed_rhs_post_step apply rhs");
    negcheck2(this->mem->advectee(ix::th)(this->ijk), this->rhs.at(ix::th)(this->ijk), "th after mixed_rhs_post_step apply rhs (+ output of th rhs)");
    negcheck2(this->mem->advectee(ix::rv)(this->ijk), this->rhs.at(ix::rv)(this->ijk), "rv after mixed_rhs_post_step apply rhs (+ output of rv rhs)");
  }

  // called before record_all() but by all openmp ranks
  void hook_ante_record_all() override
  {
    parent_t::hook_ante_record_all();
    this->reconstruct_refinee(ix::th);
    this->reconstruct_refinee(ix::rv);
  }

  virtual void diag()
  {
    assert(this->rank == 0);
    this->record_aux_dsc("radiative_flux", radiative_flux); 

    auto conv_fctr_sens = (cmn::moist_air::c_pd<real_t>() * si::kilograms * si::kelvins / si::joules);
    surf_flux_tmp = - surf_flux_sens * conv_fctr_sens;
    this->record_aux_dsc("sensible surface flux", surf_flux_tmp, true); 

    auto conv_fctr_lat = (cmn::const_cp::l_tri<real_t>() * si::kilograms / si::joules);
    surf_flux_tmp = - surf_flux_lat * conv_fctr_lat;
    this->record_aux_dsc("latent surface flux", surf_flux_tmp, true); 

    get_puddle();
    for(int i=0; i < n_puddle_scalars; ++i)
    {
      real_t sum = this->mem->distmem.sum(puddle.at(static_cast<cmn::output_t>(i)));
      this->record_aux_scalar(cmn::output_names.at(static_cast<cmn::output_t>(i)), "puddle", sum);
    }

    this->record_aux_dsc_refined("refined th", this->mem->refinee(this->ix_r2r.at(ix::th)));
    this->record_aux_dsc_refined("refined rv", this->mem->refinee(this->ix_r2r.at(ix::rv)));
  } 

  void record_all()
  {
    assert(this->rank == 0);

    // plain (no xdmf) hdf5 output
    parent_t::parent_t::parent_t::parent_t::record_all();
    this->diag();
    // xmf markup
    this->write_xmfs();
  }

  public:
  // case-specific parameters
  // note dual inheritance to get profile pointers
  struct rt_params_t : parent_t::rt_params_t//, detail::profile_ptrs_t
  {
    bool rv_src = true, th_src = true, uv_src = true, w_src = true;
    bool subsidence = false,
         coriolis = false, 
         friction = false, 
         buoyancy_wet = false, 
         radiation = false,
         vel_subsidence = false; // should subsidence be also applied to velocitiy fields - False eg. in RICO
    bool rc_src = true, rr_src = true; // these two are only relevant for blk schemes, but need to be here so that Cases can have access to it
    bool nc_src = true, nr_src = true; // these two are only relevant for blk_2m, but need to be here so that Cases can have access to them
    typename ct_params_t::real_t dz; // vertical grid size
    detail::ForceParameters_t ForceParameters;
    user_params_t user_params; // copy od user_params
    detail::profiles_t profs, profs_ref;

    // functions for updating surface fluxes per timestep
    std::function<void(typename parent_t::arr_t, typename parent_t::arr_t, typename parent_t::arr_t, const real_t&, int, const real_t&, const real_t&, const real_t&)> update_surf_flux_uv, update_surf_flux_sens, update_surf_flux_lat;
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
    tmp1(args.mem->tmp[__FILE__][0][0]),
    r_l(args.mem->tmp[__FILE__][0][1]),
    F(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    radiative_flux(args.mem->tmp[__FILE__][0][5]),
    diss_rate(args.mem->tmp[__FILE__][0][6]),
    surf_flux_sens(args.mem->tmp[__FILE__][1][0]),
    surf_flux_lat(args.mem->tmp[__FILE__][1][1]),
    surf_flux_zero(args.mem->tmp[__FILE__][1][2]),
    U_ground(args.mem->tmp[__FILE__][1][3]),
    surf_flux_tmp(args.mem->tmp[__FILE__][1][4]),
    surf_flux_u(args.mem->tmp[__FILE__][1][5]),
    surf_flux_v(args.mem->tmp[__FILE__][1][6]) // flux_v needs to be last
  {
    k_i.resize(this->shape(this->hrzntl_subdomain)); 
    k_i.reindexSelf(this->base(this->hrzntl_subdomain));
    r_l = 0.;
    surf_flux_zero = 0.;
    th_mean_prof.resize(this->vert_rng.length());
    rv_mean_prof.resize(this->vert_rng.length());
  }

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 7); // tmp1, tmp2, r_l, alpha, beta, F, diss_rate, radiative_flux
    parent_t::alloc_tmp_sclr(mem, __FILE__, n_flxs+3, "", true); // surf_flux sens/lat/hori_vel/zero/tmp, U_ground
  }
};

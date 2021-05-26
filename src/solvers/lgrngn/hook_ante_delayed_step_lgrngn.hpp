#pragma once
#include "../slvr_lgrngn.hpp"
#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif
#include "../../detail/func_time.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_ante_delayed_step()
{
  parent_t::hook_ante_delayed_step();
  if (this->rank == 0) 
  {
    // assuring previous sync step finished ...
#if defined(STD_FUTURE_WORKS)
    if (
      params.async
    ) {
      assert(ftr.valid());
#if defined(UWLCM_TIMING)
      tbeg = parent_t::clock::now();
#endif
#if defined(UWLCM_TIMING)
      parent_t::tsync_gpu += ftr.get();
#else
      ftr.get();
#endif
#if defined(UWLCM_TIMING)
      tend = parent_t::clock::now();
      parent_t::tsync_wait += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
#endif
    } else assert(!ftr.valid()); 
#endif
  }
  this->mem->barrier();

  // add microphysics contribution to th and rv
  if(params.cloudph_opts.cond)
  {
    // with cyclic bcond, th and rv in corresponding edge cells needs to change by the same amount
    this->avg_edge_sclr(rv_post_cond, this->ijk);
    this->avg_edge_sclr(th_post_cond, this->ijk);

    this->state(ix::rv)(this->ijk) += rv_post_cond(this->ijk) - rv_pre_cond(this->ijk); 
    this->state(ix::th)(this->ijk) += th_post_cond(this->ijk) - th_pre_cond(this->ijk); 
    // microphysics could have led to rv < 0 ?
    negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv after condensation");
    nancheck(this->mem->advectee(ix::th)(this->ijk), "th after condensation");
    nancheck(this->mem->advectee(ix::rv)(this->ijk), "rv after condensation");
    negcheck(this->mem->advectee(ix::th)(this->ijk), "th after condensation");
    negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after condensation");
  }

  // store liquid water content (post-cond, pre-adve and pre-subsidence)
  diag_rl();
  if(ct_params_t::sgs_scheme == libmpdataxx::solvers::smg)
    diag_rc();
    
  if (this->rank == 0) 
  {
    // running asynchronous stuff
    {
      using libcloudphxx::lgrngn::particles_t;
      using libcloudphxx::lgrngn::CUDA;
      using libcloudphxx::lgrngn::multi_CUDA;
#if defined(UWLCM_TIMING)
      tbeg = parent_t::clock::now();
#endif
#if defined(STD_FUTURE_WORKS)
      if (params.async)
      {
        assert(!ftr.valid());
        if(params.backend == CUDA)
          ftr = async_timing_launcher<typename parent_t::clock, typename parent_t::timer>(
            &particles_t<real_t, CUDA>::step_async, 
            dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
            params.cloudph_opts
          );
        else if(params.backend == multi_CUDA)
          ftr = async_timing_launcher<typename parent_t::clock, typename parent_t::timer>(
            &particles_t<real_t, multi_CUDA>::step_async, 
            dynamic_cast<particles_t<real_t, multi_CUDA>*>(prtcls.get()),
            params.cloudph_opts
          );
        assert(ftr.valid());
      } else 
#endif
        prtcls->step_async(params.cloudph_opts);

#if defined(UWLCM_TIMING)
      tend = parent_t::clock::now();
      parent_t::tasync += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
#endif
    }
  }
  this->mem->barrier();

  // subsidence of rl
  // TODO: very similar code to subsidence function in forcings.hppp
  if(params.subsidence)
  {
    auto& tmp1(this->tmp1);
    auto& r_l(this->r_l);
    const auto& ijk(this->ijk);
    const auto& params(this->params);
    auto& F(this->F);

    tmp1(ijk) = r_l(ijk);
    // fill halos for gradient calculation
    // TODO: no need to xchng in horizontal, which potentially causes MPI communication
    this->xchng_sclr(tmp1, this->ijk, this->halo);
    this->vert_grad_cnt(tmp1, F, params.dz);
    F(ijk).reindex(this->zero) *= - (*params.w_LS)(this->vert_idx);
    r_l(ijk) += F(ijk) * this->dt;

    tmp1(ijk) = r_c(ijk);
    // fill halos for gradient calculation
    // TODO: no need to xchng in horizontal, which potentially causes MPI communication
    this->xchng_sclr(tmp1, this->ijk, this->halo);
    this->vert_grad_cnt(tmp1, F, params.dz);
    F(ijk).reindex(this->zero) *= - (*params.w_LS)(this->vert_idx);
    r_c(ijk) += F(ijk) * this->dt;
  }

  // advect r_l and r_c (1st-order)
  this->self_advec_donorcell(this->r_l);
  this->self_advec_donorcell(this->r_c);
  negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv at the end of ante delayed step");
}

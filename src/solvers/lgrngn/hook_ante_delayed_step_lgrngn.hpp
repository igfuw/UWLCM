#pragma once
#include "../slvr_lgrngn.hpp"
#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif
#include "../../detail/func_time.hpp"
#include <libmpdata++/formulae/refined_grid.hpp>

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
      tbeg = setup::clock::now();
#endif
#if defined(UWLCM_TIMING)
      parent_t::tsync_gpu += ftr.get();
#else
      ftr.get();
#endif
#if defined(UWLCM_TIMING)
      tend = setup::clock::now();
      parent_t::tsync_wait += std::chrono::duration_cast<setup::timer>( tend - tbeg );
#endif
    } else assert(!ftr.valid()); 
#endif
  }
  this->mem->barrier();

  // add microphysics contribution to th and rv
  if(params.cloudph_opts.cond)
  {
    // calculate th and rv changes as averages from the refined grid

//    std::cerr << "rv post cond: " << rv_post_cond(this->ijk_ref) << std::endl;
//    std::cerr << "th post cond: " << th_post_cond(this->ijk_ref) << std::endl;

    // refined drv and dth stored in _post_cond
    rv_post_cond(this->ijk_ref) = rv_post_cond(this->ijk_ref) - rv_pre_cond(this->ijk_ref); 
    th_post_cond(this->ijk_ref) = th_post_cond(this->ijk_ref) - th_pre_cond(this->ijk_ref); 

//    std::cerr << "rv diff cond: " << rv_post_cond(this->ijk_ref) << std::endl;
//    std::cerr << "th diff cond: " << th_post_cond(this->ijk_ref) << std::endl;

    this->mem->barrier();
    this->spatial_average_ref2reg(th_post_cond, dth);
    this->spatial_average_ref2reg(rv_post_cond, drv);

 //   std::cerr << "dth: " << dth(this->ijk) << std::endl;
 //   std::cerr << "drv: " << drv(this->ijk) << std::endl;

    // with cyclic bcond, th and rv in corresponding edge cells needs to change by the same amount
//    this->avg_edge_sclr(dth, this->ijk);
//    this->avg_edge_sclr(drv, this->ijk);

    this->state(ix::th)(this->ijk) += dth(this->ijk); 
    this->state(ix::rv)(this->ijk) += drv(this->ijk); 

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
      tbeg = setup::clock::now();
#endif
#if defined(STD_FUTURE_WORKS)
      if (params.async)
      {
        assert(!ftr.valid());
        if(params.backend == CUDA)
          ftr = async_launcher(
            &particles_t<real_t, CUDA>::step_async, 
            dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
            params.cloudph_opts
          );
        else if(params.backend == multi_CUDA)
          ftr = async_launcher(
            &particles_t<real_t, multi_CUDA>::step_async, 
            dynamic_cast<particles_t<real_t, multi_CUDA>*>(prtcls.get()),
            params.cloudph_opts
          );
        assert(ftr.valid());
      } else 
#endif
        prtcls->step_async(params.cloudph_opts);

#if defined(UWLCM_TIMING)
      tend = setup::clock::now();
      parent_t::tasync += std::chrono::duration_cast<setup::timer>( tend - tbeg );
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
    F(ijk).reindex(this->zero) *= - (params.profs.w_LS)(this->vert_idx);
    r_l(ijk) += F(ijk) * this->dt;

    tmp1(ijk) = r_c(ijk);
    // fill halos for gradient calculation
    // TODO: no need to xchng in horizontal, which potentially causes MPI communication
    this->xchng_sclr(tmp1, this->ijk, this->halo);
    this->vert_grad_cnt(tmp1, F, params.dz);
    F(ijk).reindex(this->zero) *= - (params.profs.w_LS)(this->vert_idx);
    r_c(ijk) += F(ijk) * this->dt;
  }

  // advect r_l and r_c (1st-order)
  this->self_advec_donorcell(this->r_l);
  this->self_advec_donorcell(this->r_c);
  negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv at the end of ante delayed step");
}

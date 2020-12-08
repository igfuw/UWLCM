#pragma once
#include "../slvr_lgrngn.hpp"
#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_mixed_rhs_ante_step()
{
  negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv at start of mixed_rhs_ante_step");

  rv_pre_cond(this->ijk) = this->state(ix::rv)(this->ijk); 
  th_pre_cond(this->ijk) = this->state(ix::th)(this->ijk); 

  this->mem->barrier();

  // pass Eulerian fields to microphysics 
  if (this->rank == 0) 
  {
    // temporarily Cx & Cz are multiplied by this->rhod ...
    auto 
      Cx = this->mem->GC[0](this->Cx_domain).copy(),
      Cy = this->mem->GC[1](this->Cy_domain).copy(),
      Cz = this->mem->GC[ix::w](this->Cz_domain).copy(); 
    nancheck(Cx, "Cx after copying from mpdata");
    nancheck(Cy, "Cy after copying from mpdata");
    nancheck(Cz, "Cz after copying from mpdata");

    // ... and now dividing them by this->rhod (TODO: z=0 is located at k=1/2)
    {
      Cx.reindex(this->zero) /= (*params.rhod)(this->vert_idx);
      Cy.reindex(this->zero) /= (*params.rhod)(this->vert_idx);
      Cz.reindex(this->zero) /= (*params.rhod)(this->vert_idx); // TODO: should be interpolated, since theres a shift between positions of rhod and Cz
    }

    // start synchronous stuff timer
#if defined(UWLCM_TIMING)
    tbeg = parent_t::clock::now();
#endif

    using libcloudphxx::lgrngn::particles_t;
    using libcloudphxx::lgrngn::CUDA;
    using libcloudphxx::lgrngn::multi_CUDA;

    prtcls->sync_in(
      make_arrinfo(this->mem->advectee(ix::th)),
      make_arrinfo(this->mem->advectee(ix::rv)),
      libcloudphxx::lgrngn::arrinfo_t<real_t>(),
      make_arrinfo(Cx),
      this->n_dims == 2 ? libcloudphxx::lgrngn::arrinfo_t<real_t>() : make_arrinfo(Cy),
      make_arrinfo(Cz),
      (ct_params_t::sgs_scheme == libmpdataxx::solvers::iles) || (!params.cloudph_opts.turb_cond && !params.cloudph_opts.turb_adve && !params.cloudph_opts.turb_coal) ?
                                  libcloudphxx::lgrngn::arrinfo_t<real_t>() :
                                  make_arrinfo(this->diss_rate(this->domain).reindex(this->zero))
    );

    // start sync/async run of step_cond
    // step_cond takes th and rv only for sync_out purposes - the values of th and rv before condensation come from sync_in, i.e. before apply_rhs

#if defined(STD_FUTURE_WORKS)
    if (params.async)
    {
      assert(!ftr.valid());
      if(params.backend == CUDA)
        ftr = std::async(
          std::launch::async, 
          &particles_t<real_t, CUDA>::step_cond, 
          dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
          params.cloudph_opts,
          make_arrinfo(th_post_cond(this->domain).reindex(this->zero)),
          make_arrinfo(rv_post_cond(this->domain).reindex(this->zero)),
          std::map<enum libcloudphxx::common::chem::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> >()
        );
      else if(params.backend == multi_CUDA)
        ftr = std::async(
          std::launch::async, 
          &particles_t<real_t, multi_CUDA>::step_cond, 
          dynamic_cast<particles_t<real_t, multi_CUDA>*>(prtcls.get()),
          params.cloudph_opts,
          make_arrinfo(th_post_cond(this->domain).reindex(this->zero)),
          make_arrinfo(rv_post_cond(this->domain).reindex(this->zero)),
          std::map<enum libcloudphxx::common::chem::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> >()
        );
      assert(ftr.valid());
    } else 
#endif
    {
      prtcls->step_cond(
        params.cloudph_opts,
        make_arrinfo(th_post_cond(this->domain).reindex(this->zero)),
        make_arrinfo(rv_post_cond(this->domain).reindex(this->zero))
      );
    }

#if defined(UWLCM_TIMING)
    tend = parent_t::clock::now();
    parent_t::tsync += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
#endif
  }
  this->mem->barrier();

  parent_t::hook_mixed_rhs_ante_step();
}

#pragma once
#include "../slvr_lgrngn.hpp"
#include "../../detail/func_time.hpp"
#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_mixed_rhs_ante_step()
{
/*
  courants[0] = NAN;
  courants[1] = NAN;
  courants[2] = NAN;

  courants[0] = 0;
  courants[1] = 0;
  courants[2] = 0;

  this->mem->barrier();

  courants[0](this->ijk_ref) = 0;
  courants[1](this->ijk_ref) = 0;
  courants[2](this->ijk_ref) = 0;

  this->mem->barrier();
  */

  this->interpolate_refinee(ix::u);
  this->interpolate_refinee(ix::v);
  this->interpolate_refinee(ix::w);

  this->interpolate_refined_courants(courants, this->uvw_ref);

  //this->xchng_ref(ix::u);
  //this->xchng_ref(ix::v);
  //this->xchng_ref(ix::w);

//  courants[0](this->ijkm_ref) = this->mem->psi_ref.at(this->ix_r2r.at(ix::u))(this->ijkm_ref);
//  courants[0](this->ijkm_ref) = this->mem->refinee(this->ix_r2r.at(ix::u))(this->ijkm_ref) * this->dt / this->di;
//  courants[1](this->ijkm_ref) = this->mem->refinee(this->ix_r2r.at(ix::v))(this->ijkm_ref) * this->dt / this->dj;
//  courants[2](this->ijkm_ref) = this->mem->refinee(this->ix_r2r.at(ix::w))(this->ijkm_ref) * this->dt / this->dk;



  negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv at start of mixed_rhs_ante_step");

  this->generate_stretching_parameters(std::random_device{}());
  this->reconstruct_refinee(ix::th);
  this->reconstruct_refinee(ix::rv);

  negtozero(this->mem->refinee(this->ix_r2r.at(ix::rv))(this->ijk_ref), "refined rv at start of mixed_rhs_ante_step");

  rv_pre_cond(this->ijk_ref) = this->mem->refinee(this->ix_r2r.at(ix::rv))(this->ijk_ref); 
  th_pre_cond(this->ijk_ref) = this->mem->refinee(this->ix_r2r.at(ix::th))(this->ijk_ref); 

  this->mem->barrier();

  // pass Eulerian fields to microphysics 
  if (this->rank == 0) 
  {

    nancheck(courants[0], "courants[0] after interpolation from refined velocities");
    nancheck(courants[1], "courants[1] after interpolation from refined velocities");
    nancheck(courants[2], "courants[2] after interpolation from refined velocities");

    // temporarily Cx & Cz are multiplied by this->rhod ...
    auto 
      Cx = this->mem->GC[0](this->Cx_domain).copy(),
      Cy = this->mem->GC[1](this->Cy_domain).copy(), // TODO: no need to copy in 2D
      Cz = this->mem->GC[ix::w](this->Cz_domain).copy(); 
    nancheck(Cx, "Cx after copying from mpdata");
    nancheck(Cy, "Cy after copying from mpdata");
    nancheck(Cz, "Cz after copying from mpdata");

    // ... and now dividing them by this->rhod (TODO: z=0 is located at k=1/2)
    {
      Cx.reindex(this->zero) /= (params.profs.rhod)(this->vert_idx);
      Cy.reindex(this->zero) /= (params.profs.rhod)(this->vert_idx);
      Cz.reindex(this->zero) /= (params.profs.rhod)(this->vert_idx); // TODO: should be interpolated, since theres a shift between positions of rhod and Cz
    }

    std::cerr << "courants[0]: " << courants[0] << std::endl;
    std::cerr << "Cx: " << Cx << std::endl;

    std::cerr << "this->mem->refinee(this->ix_r2r.at(ix::u)): " << this->mem->refinee(this->ix_r2r.at(ix::u)) << std::endl;
    std::cerr << "this->mem->advectee(ix::u): " << this->mem->advectee(ix::u) << std::endl;

    // assuring previous async step finished ...
#if defined(STD_FUTURE_WORKS)
    if (
      params.async && 
      this->timestep != 0 && // ... but not in first timestep ...
      ((this->timestep ) % this->outfreq != 0) // ... and not after diag call, note: timestep is updated after ante_step
    ) {
      assert(ftr.valid());
#if defined(UWLCM_TIMING)
      tbeg = parent_t::clock::now();
#endif
#if defined(UWLCM_TIMING)
      parent_t::tasync_gpu += ftr.get();
#else
      ftr.get();
#endif
#if defined(UWLCM_TIMING)
      tend = parent_t::clock::now();
      parent_t::tasync_wait += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
#endif
    } else assert(!ftr.valid()); 
#endif

    // change src and rlx flags after the first step. needs to be done after async finished, because async uses opts reference
    if(this->timestep == 1)
    {
      // turn off aerosol src, because it was only used to initialize gccn below some height
      params.cloudph_opts.src = false;
      // if relaxation is to be done, turn it on after gccn were created by src
      if(params.user_params.relax_ccn)
        params.cloudph_opts.rlx = true;
    }

    // start synchronous stuff timer
#if defined(UWLCM_TIMING)
    tbeg = parent_t::clock::now();
#endif

    using libcloudphxx::lgrngn::particles_t;
    using libcloudphxx::lgrngn::CUDA;
    using libcloudphxx::lgrngn::multi_CUDA;

/*
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
    */

/*
    std::cerr << "Cx: " << Cx << std::endl;
    std::cerr << "Cy: " << Cy << std::endl;
    std::cerr << "Cz: " << Cz << std::endl;

    std::cerr << "curants[0]: " << courants[0] << std::endl;
    std::cerr << "curants[1]: " << courants[1] << std::endl;
    std::cerr << "curants[2]: " << courants[2] << std::endl;

    // TODO: make courants have the desired sizes (?)
    auto 
      Cx_ref = courants[0].copy(),
      Cy_ref = courants[1].copy(),
      Cz_ref = courants[2].copy();


    std::cerr << "Cx_ref: " << Cx_ref << std::endl;
    std::cerr << "Cy_ref: " << Cy_ref << std::endl;
    std::cerr << "Cz_ref: " << Cz_ref << std::endl;
    */

    prtcls->sync_in(
      make_arrinfo(this->mem->refinee(this->ix_r2r.at(ix::th))),
      make_arrinfo(this->mem->refinee(this->ix_r2r.at(ix::rv))),
      libcloudphxx::lgrngn::arrinfo_t<real_t>(),
      make_arrinfo(courants[0]),
      this->n_dims == 2 ? libcloudphxx::lgrngn::arrinfo_t<real_t>() : make_arrinfo(courants[1]),
      make_arrinfo(courants[2])
      /*,
      (ct_params_t::sgs_scheme == libmpdataxx::solvers::iles) || (!params.cloudph_opts.turb_cond && !params.cloudph_opts.turb_adve && !params.cloudph_opts.turb_coal) ?
                                  libcloudphxx::lgrngn::arrinfo_t<real_t>() :
                                  make_arrinfo(this->diss_rate(this->domain).reindex(this->zero))
                                  */
    );

    // start sync/async run of step_cond
    // step_cond takes th and rv only for sync_out purposes - the values of th and rv before condensation come from sync_in, i.e. before apply_rhs

#if defined(STD_FUTURE_WORKS)
    if (params.async)
    {
      assert(!ftr.valid());
      if(params.backend == CUDA)
  #if defined(UWLCM_TIMING)
        ftr = async_timing_launcher<typename parent_t::clock, typename parent_t::timer>(
  #else
        ftr = std::async(std::launch::async,
  #endif
          &particles_t<real_t, CUDA>::step_cond, 
          dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
          params.cloudph_opts,
          make_arrinfo(th_post_cond(this->domain_ref).reindex(this->zero)),
          make_arrinfo(rv_post_cond(this->domain_ref).reindex(this->zero)),
          std::map<enum libcloudphxx::common::chem::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> >()
        );
      else if(params.backend == multi_CUDA)
  #if defined(UWLCM_TIMING)
        ftr = async_timing_launcher<typename parent_t::clock, typename parent_t::timer>(
  #else
        ftr = std::async(std::launch::async,
  #endif
          &particles_t<real_t, multi_CUDA>::step_cond, 
          dynamic_cast<particles_t<real_t, multi_CUDA>*>(prtcls.get()),
          params.cloudph_opts,
          make_arrinfo(th_post_cond(this->domain_ref).reindex(this->zero)),
          make_arrinfo(rv_post_cond(this->domain_ref).reindex(this->zero)),
          std::map<enum libcloudphxx::common::chem::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> >()
        );
      assert(ftr.valid());
    } else 
#endif
    {
      prtcls->step_cond(
        params.cloudph_opts,
        make_arrinfo(th_post_cond(this->domain_ref).reindex(this->zero)),
        make_arrinfo(rv_post_cond(this->domain_ref).reindex(this->zero))
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

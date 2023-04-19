#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::step_cond()
{
  assert(this->rank == 0);

  using libcloudphxx::lgrngn::particles_t;
  using libcloudphxx::lgrngn::CUDA;
  using libcloudphxx::lgrngn::multi_CUDA;

  // start sync/async run of step_cond
  // step_cond takes th and rv only for sync_out purposes - the values of th and rv before condensation come from sync_in, i.e. before apply_rhs

#if defined(STD_FUTURE_WORKS)
  if (params.async)
  {
    assert(!ftr.valid());
    if(params.backend == CUDA)
      ftr = async_launcher(
        &particles_t<real_t, CUDA>::step_cond, 
        dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
        params.cloudph_opts,
        make_arrinfo(th_post_cond(this->domain).reindex(this->zero)),
        make_arrinfo(rv_post_cond(this->domain).reindex(this->zero)),
        std::map<enum libcloudphxx::common::chem::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> >()
      );
    else if(params.backend == multi_CUDA)
      ftr = async_launcher(
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
}

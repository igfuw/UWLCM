#pragma once
#include "../slvr_lgrngn_chem.hpp"

template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::step_cond()
{
  assert(this->rank == 0);

  using libcloudphxx::lgrngn::particles_t;
  using libcloudphxx::lgrngn::CUDA;
  using libcloudphxx::lgrngn::multi_CUDA;
  
  boost::assign::insert(this->ambient_chem_post_cond)
    (chem_species_t::SO2,  this->make_arrinfo(SO2_post_cond))
    (chem_species_t::O3,   this->make_arrinfo(O3_post_cond))
    (chem_species_t::H2O2, this->make_arrinfo(H2O2_post_cond))
    (chem_species_t::CO2,  this->make_arrinfo(CO2_post_cond))
    (chem_species_t::NH3,  this->make_arrinfo(NH3_post_cond))
    (chem_species_t::HNO3, this->make_arrinfo(HNO3_post_cond));

  // start sync/async run of step_cond
  // step_cond takes th and rv only for sync_out purposes - the values of th and rv before condensation come from sync_in, i.e. before apply_rhs

#if defined(STD_FUTURE_WORKS)
  if (params.async)
  {
    assert(!this->ftr.valid());
    if(params.backend == CUDA)
      this->ftr = async_launcher(
        &particles_t<real_t, CUDA>::step_cond, 
        dynamic_cast<particles_t<real_t, CUDA>*>(this->prtcls.get()),
        params.cloudph_opts,
        this->make_arrinfo(this->th_post_cond(this->domain).reindex(this->zero)),
        this->make_arrinfo(this->rv_post_cond(this->domain).reindex(this->zero)),
        ambient_chem_post_cond
      );
    else if(params.backend == multi_CUDA)
      this->ftr = async_launcher(
        &particles_t<real_t, multi_CUDA>::step_cond, 
        dynamic_cast<particles_t<real_t, multi_CUDA>*>(this->prtcls.get()),
        params.cloudph_opts,
        this->make_arrinfo(this->th_post_cond(this->domain).reindex(this->zero)),
        this->make_arrinfo(this->rv_post_cond(this->domain).reindex(this->zero)),
        ambient_chem_post_cond
      );
    assert(this->ftr.valid());
  } else 
#endif
  {
    this->prtcls->step_cond(
      params.cloudph_opts,
      this->make_arrinfo(this->th_post_cond(this->domain).reindex(this->zero)),
      this->make_arrinfo(this->rv_post_cond(this->domain).reindex(this->zero)),
      ambient_chem_post_cond
    );
      std::cerr <<  "SO2_post_cond right after step_cond with async==0: " << SO2_post_cond(this->ijk) << std::endl;

  }
}

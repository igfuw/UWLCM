#pragma once
#include "../slvr_lgrngn_chem.hpp"

template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::sync_e2l()
{
  assert(this->rank==0);

  this->prtcls->sync_in(
    this->make_arrinfo(this->mem->advectee(ix::th)),
    this->make_arrinfo(this->mem->advectee(ix::rv)),
    libcloudphxx::lgrngn::arrinfo_t<real_t>(),
    this->make_arrinfo(this->Cx),
    this->n_dims == 2 ? libcloudphxx::lgrngn::arrinfo_t<real_t>() : this->make_arrinfo(this->Cy),
    this->make_arrinfo(this->Cz),
    (ct_params_t::sgs_scheme == libmpdataxx::solvers::iles) || (!params.cloudph_opts.turb_cond && !params.cloudph_opts.turb_adve && !params.cloudph_opts.turb_coal) ?
      libcloudphxx::lgrngn::arrinfo_t<real_t>() :
      this->make_arrinfo(this->diss_rate(this->domain).reindex(this->zero)),
    ambient_chem
  );
}

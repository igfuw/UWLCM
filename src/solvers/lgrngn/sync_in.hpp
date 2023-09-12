#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::sync_in()
{
  assert(this->rank == 0);

  // pass euler variables to lagrangian microphysics
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
}

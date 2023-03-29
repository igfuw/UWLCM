#pragma once
#include "../slvr_lgrngn_chem.hpp"

template <class ct_params_t>
virtual void slvr_lgrngn_chem<ct_params_t>::sync_e2l()
{
  assert(this->rank==0);

  std::map<enum chem_species_t, const libcloudphxx::lgrngn::arrinfo_t<real_t> > ambient_chem;
  boost::assign::insert(ambient_chem_init)
    (chem_species_t::SO2,  this->make_arrinfo(this->mem->advectee(ix::SO2g)))
    (chem_species_t::O3,   this->make_arrinfo(this->mem->advectee(ix::O3g)))
    (chem_species_t::H2O2, this->make_arrinfo(this->mem->advectee(ix::H2O2g)))
    (chem_species_t::CO2,  this->make_arrinfo(this->mem->advectee(ix::CO2g)))
    (chem_species_t::NH3,  this->make_arrinfo(this->mem->advectee(ix::NH3g)))
    (chem_species_t::HNO3, this->make_arrinfo(this->mem->advectee(ix::HNO3g)));

  this->prtcls->sync_in(
    make_arrinfo(this->mem->advectee(ix::th)),
    make_arrinfo(this->mem->advectee(ix::rv)),
    libcloudphxx::lgrngn::arrinfo_t<real_t>(),
    make_arrinfo(this->Cx),
    this->n_dims == 2 ? libcloudphxx::lgrngn::arrinfo_t<real_t>() : make_arrinfo(this->Cy),
    make_arrinfo(this->Cz),
    (ct_params_t::sgs_scheme == libmpdataxx::solvers::iles) || (!params.cloudph_opts.turb_cond && !params.cloudph_opts.turb_adve && !params.cloudph_opts.turb_coal) ?
      libcloudphxx::lgrngn::arrinfo_t<real_t>() :
      make_arrinfo(this->diss_rate(this->domain).reindex(this->zero)),
    ambient_chem
  );
}

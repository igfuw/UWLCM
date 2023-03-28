#pragma once
#include "../slvr_lgrngn_chem.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::init_prtcls() override
{
  assert(this->rank == 0);
  assert(params.cloudph_opts_init.chem_switch == true);

  // temporary array of densities - prtcls cant be init'd with 1D profile
  typename parent_t::arr_t rhod(this->mem->advectee(ix::th).shape()); // TODO: replace all rhod arrays with this->mem->G
  rhod = (*params.rhod)(this->vert_idx);

  // temporary array of pressure - prtcls cant be init'd with 1D profile
  typename parent_t::arr_t p_e(this->mem->advectee(ix::th).shape());
  p_e = (*params.p_e)(this->vert_idx);

  std::map<enum chem_species_t, const libcloudphxx::lgrngn::arrinfo_t<real_t> > ambient_chem_init;
  boost::assign::insert(ambient_chem_init)
    (chem_species_t::SO2,  this->make_arrinfo(this->mem->advectee(ix::SO2g)))
    (chem_species_t::O3,   this->make_arrinfo(this->mem->advectee(ix::O3g)))
    (chem_species_t::H2O2, this->make_arrinfo(this->mem->advectee(ix::H2O2g)))
    (chem_species_t::CO2,  this->make_arrinfo(this->mem->advectee(ix::CO2g)))
    (chem_species_t::NH3,  this->make_arrinfo(this->mem->advectee(ix::NH3g)))
    (chem_species_t::HNO3, this->make_arrinfo(this->mem->advectee(ix::HNO3g)));

  this->prtcls->init(
    this->make_arrinfo(this->mem->advectee(ix::th)),
    this->make_arrinfo(this->mem->advectee(ix::rv)),
    this->make_arrinfo(rhod),
    this->make_arrinfo(p_e)
    libcloudphxx::lgrngn::arrinfo_t<real_t>(),
    libcloudphxx::lgrngn::arrinfo_t<real_t>(),
    libcloudphxx::lgrngn::arrinfo_t<real_t>(),
    ambient_chem_init
  );
}

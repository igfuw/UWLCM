#pragma once
#include "../slvr_lgrngn_chem.hpp"

template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::hook_ante_loop(int nt)
{
  params.flag_chem_rct = params.cloudph_opts.chem_rct;

/*
    this->coal = parent_t::params.cloudph_opts.coal;
    this->sedi = parent_t::params.cloudph_opts.sedi;
    {
      blitz::secondIndex j;
      namespace molar_mass  = libcloudphxx::common::molar_mass;

      // initialise trace gases profiles
      this->state(ix::SO2g)  =
        config::mixr_helper(this->setup)(j * parent_t::params.dz)
        * (this->setup.SO2_g_0  * molar_mass::M_SO2<real_t>()  * si::moles / si::kilograms);
      this->state(ix::O3g)   =
        config::mixr_helper(this->setup)(j * parent_t::params.dz)
        * (this->setup.O3_g_0   * molar_mass::M_O3<real_t>()   * si::moles / si::kilograms);
      this->state(ix::H2O2g) =
        config::mixr_helper(this->setup)(j * parent_t::params.dz)
        * (this->setup.H2O2_g_0 * molar_mass::M_H2O2<real_t>() * si::moles / si::kilograms);
      this->state(ix::CO2g)  =
        config::mixr_helper(this->setup)(j * parent_t::params.dz)
        * (this->setup.CO2_g_0  * molar_mass::M_CO2<real_t>()  * si::moles / si::kilograms);
      this->state(ix::NH3g)  =
        config::mixr_helper(this->setup)(j * parent_t::params.dz)
        * (this->setup.NH3_g_0  * molar_mass::M_NH3<real_t>()  * si::moles / si::kilograms);
      this->state(ix::HNO3g) =
        config::mixr_helper(this->setup)(j * parent_t::params.dz)
        * (this->setup.HNO3_g_0 * molar_mass::M_HNO3<real_t>() * si::moles / si::kilograms);
    }
    */

    parent_t::hook_ante_loop(nt);

    if (this->rank == 0)
    {
      this->record_aux_const("flag_chem_rct", "lgrngn_chem", params.flag_chem_rct);
      this->record_aux_const("chem_rct", "lgrngn_chem", params.cloudph_opts.chem_rct);
      this->record_aux_const("chem_dsc", "lgrngn_chem", params.cloudph_opts.chem_dsc);
      this->record_aux_const("chem_dsl", "lgrngn_chem", params.cloudph_opts.chem_dsl);
    }
}

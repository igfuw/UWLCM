#pragma once
#include "../slvr_lgrngn_chem.hpp"

// recording mass of H and S_VI in wet radius bins
template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::diag_pH()
{
  assert(this->rank == 0);
  {
    int rng_num = 0;
    for (auto &rng_moms : params.out_wet_pH)
    {
      //wet
      auto &rng(rng_moms.first);
      parent_t::prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
      for (auto &mom : rng_moms.second)
      {
        parent_t::prtcls->diag_chem(chem_species_t::S_VI);
        this->record_aux(this->aux_name("chem_S_VI_rw", rng_num, mom), parent_t::prtcls->outbuf());

        parent_t::prtcls->diag_chem(chem_species_t::H);
        this->record_aux(this->aux_name("chem_H_rw", rng_num, mom), parent_t::prtcls->outbuf());
      }
      rng_num++;
    }
  }
}

template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::diag_chem()
{
  assert(this->rank == 0);
  {
    // chem
    for (auto &rng_moms : params.out_chem)
    {
      auto &rng(rng_moms.first);
      parent_t::prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);

      //TODO: for (auto &chem_enum : chem_species_t)
      {
        parent_t::prtcls->diag_chem(chem_species_t::SO2);
        this->record_aux("chem_S_IV_aq", parent_t::prtcls->outbuf());
        parent_t::prtcls->diag_chem(chem_species_t::S_VI);
        this->record_aux("chem_S_VI_aq", parent_t::prtcls->outbuf());

        parent_t::prtcls->diag_chem(chem_species_t::O3);
        this->record_aux("chem_O3_aq", parent_t::prtcls->outbuf());
        parent_t::prtcls->diag_chem(chem_species_t::H2O2);
        this->record_aux("chem_H2O2_aq", parent_t::prtcls->outbuf());

        parent_t::prtcls->diag_chem(chem_species_t::H);
        this->record_aux("chem_H_aq", parent_t::prtcls->outbuf());

        parent_t::prtcls->diag_chem(chem_species_t::CO2);
        this->record_aux("chem_C_IV_aq", parent_t::prtcls->outbuf());

        parent_t::prtcls->diag_chem(chem_species_t::NH3);
        this->record_aux("chem_N_III_aq", parent_t::prtcls->outbuf());

        parent_t::prtcls->diag_chem(chem_species_t::HNO3);
        this->record_aux("chem_N_V_aq", parent_t::prtcls->outbuf());
      }
    }
  }
}

template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::diag()
{
  parent_t::diag();
  diag_pH();
  diag_chem();
}

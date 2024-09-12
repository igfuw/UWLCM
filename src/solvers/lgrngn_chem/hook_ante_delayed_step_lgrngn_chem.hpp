#pragma once
#include "../slvr_lgrngn_chem.hpp"

template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::hook_ante_delayed_step()
{
  parent_t::hook_ante_delayed_step();

  // add microphysics contribution to chem species
  //if(params.cloudph_opts.chem)
  {
    std::array chem_syncout = {
      std::tuple{ix::SO2g,  SO2_pre_cond,  SO2_post_cond,  "SO2g"},
      std::tuple{ix::CO2g,  CO2_pre_cond,  CO2_post_cond,  "CO2g"},
      std::tuple{ix::O3g,   O3_pre_cond,   O3_post_cond,   "O3g"},
      std::tuple{ix::H2O2g, H2O2_pre_cond, H2O2_post_cond, "H2O2g"},
      std::tuple{ix::NH3g,  NH3_pre_cond,  NH3_post_cond,  "NH3g"},
      std::tuple{ix::HNO3g, HNO3_pre_cond, HNO3_post_cond, "HNO3g"}
    };

    for(auto &cs: chem_syncout)
    {
      // with cyclic bcond, chem in corresponding edge cells needs to change by the same amount
      this->avg_edge_sclr(std::get<2>(cs), this->ijk);
      this->state(std::get<0>(cs))(this->ijk) += std::get<2>(cs)(this->ijk) - std::get<1>(cs)(this->ijk); 
      nancheck(this->mem->advectee(std::get<0>(cs))(this->ijk), std::string(std::get<3>(cs)) + std::string(" after condensation"));
      negcheck(this->mem->advectee(std::get<0>(cs))(this->ijk), std::string(std::get<3>(cs)) + std::string(" after condensation"));
    }
  }
}

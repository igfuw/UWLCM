#pragma once
#include "../slvr_lgrngn_chem.hpp"

template <class ct_params_t>
void slvr_lgrngn_chem<ct_params_t>::hook_mixed_rhs_ante_step()
{
  SO2_pre_cond(this->ijk)  = this->state(ix::SO2g)(this->ijk);
  O3_pre_cond(this->ijk)   = this->state(ix::O3g)(this->ijk);
  H2O2_pre_cond(this->ijk) = this->state(ix::H2O2g)(this->ijk);
  CO2_pre_cond(this->ijk)  = this->state(ix::CO2g)(this->ijk);
  NH3_pre_cond(this->ijk)  = this->state(ix::NH3g)(this->ijk);
  HNO3_pre_cond(this->ijk) = this->state(ix::HNO3g)(this->ijk);

  parent_t::hook_mixed_rhs_ante_step();
}

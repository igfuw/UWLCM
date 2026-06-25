//TODO: move calc_forces and forcings to case class?
#pragma once
#include "../../forcings/subsidence.hpp"
#include "../common/calc_forces_common.hpp"

// single-moment bulk forcing functions
template <class ct_params_t>
void slvr_blk_1m_ice_common<ct_params_t>::ria_src()
{
  const auto &ijk = this->ijk;
  if(params.ria_src)
  {
    // large-scale vertical wind
    parent_t::subsidence(ix::ria); 

    this->alpha(ijk) = this->F(ijk);
  }
  else
    this->alpha(ijk) = 0.;

  this->beta(ijk) = 0.;
}

template <class ct_params_t>
void slvr_blk_1m_ice_common<ct_params_t>::rib_src()
{
  const auto &ijk = this->ijk;
  if(params.rib_src)
  {
    // large-scale vertical wind
    parent_t::subsidence(ix::rib); 

    this->alpha(ijk) = this->F(ijk);
  }
  else
    this->alpha(ijk) = 0.;

  this->beta(ijk) = 0.;
}

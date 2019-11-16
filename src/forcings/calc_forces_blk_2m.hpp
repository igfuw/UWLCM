//TODO: move calc_forces and forcings to case class?
#pragma once
#include "subsidence.hpp"
#include "calc_forces_common.hpp"
#include "../solvers/slvr_blk_2m.hpp"

// double-moment bulk forcing functions
// rc_src
template <class ct_params_t>
void slvr_blk_2m_common<ct_params_t>::rc_src()
{
  const auto &ijk = this->ijk;
  if(params.rc_src)
  {
    // large-scale vertical wind
    parent_t::subsidence(ix::rc); 
    this->alpha(ijk) = this->F(ijk);
  }
  else
    this->alpha(ijk) = 0.;

  this->beta(ijk) = 0.;
}

// nc_src
template <class ct_params_t>
void slvr_blk_2m_common<ct_params_t>::nc_src()
{
  const auto &ijk = this->ijk;
  if(params.nc_src)
  {
    // large-scale vertical wind
    parent_t::subsidence(ix::nc); 
    this->alpha(ijk) = this->F(ijk);
  }
  else
    this->alpha(ijk) = 0.;

  this->beta(ijk) = 0.;
}

// rr_src
template <class ct_params_t>
void slvr_blk_2m_common<ct_params_t>::rr_src()
{
  const auto &ijk = this->ijk;
  if(params.rr_src)
  {
    // large-scale vertical wind
    parent_t::subsidence(ix::rr); 
    this->alpha(ijk) = this->F(ijk);
  }
  else
    this->alpha(ijk) = 0.;

  this->beta(ijk) = 0.;
}

// nr_src
template <class ct_params_t>
void slvr_blk_2m_common<ct_params_t>::nr_src()
{
  const auto &ijk = this->ijk;
  if(params.nr_src)
  {
    // large-scale vertical wind
    parent_t::subsidence(ix::nr); 
    this->alpha(ijk) = this->F(ijk);
  }
  else
    this->alpha(ijk) = 0.;

  this->beta(ijk) = 0.;
}

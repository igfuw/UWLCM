//TODO: move calc_forces and forcings to case class?
#pragma once
#include "subsidence.hpp"
#include "calc_forces_common.hpp"
#include "../solvers/slvr_blk_2m.hpp"

// single-moment bulk forcing functions
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
  // nudging, todo: use some other coeff than vab_coeff
//  this->alpha(ijk).reindex(this->zero) += (*this->mem->vab_coeff)(ijk).reindex(this->zero) * (*params.rv_e)(this->vert_idx); // TODO: its a constant, cache it
//  this->beta(ijk) = - (*this->mem->vab_coeff)(ijk);
}


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
  // nudging, todo: use some other coeff than vab_coeff
//  this->alpha(ijk).reindex(this->zero) += (*this->mem->vab_coeff)(ijk).reindex(this->zero) * (*params.rv_e)(this->vert_idx); // TODO: its a constant, cache it
//  this->beta(ijk) = - (*this->mem->vab_coeff)(ijk);
}

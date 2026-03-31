//TODO: move calc_forces and forcings to case class?
#pragma once
#include "../../forcings/subsidence.hpp"
#include "../common/calc_forces_common.hpp"
#include "../slvr_blk_1m.hpp"

/**
 * @brief Source term for cloud water mixing ratio (rc) due to large-scale forcings
 *
 * Applies large-scale vertical wind (subsidence) and optional nudging to the
 * cloud water field. Updates the `alpha` and `beta` arrays that are used in
 * the solver for forcing calculations.
 *
 * - If `params.rc_src` is true, applies subsidence from parent class and
 *   sets `alpha` from forcing array `F`.
 * - Otherwise, `alpha` is set to zero.
 */
template <class ct_params_t>
void slvr_blk_1m_common<ct_params_t>::rc_src()
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

/**
 * @brief Source term for rain water mixing ratio (rr) due to large-scale forcings
 *
 * Applies large-scale vertical wind (subsidence) and optional nudging to the
 * rain water field. Updates the `alpha` and `beta` arrays that are used in
 * the solver for forcing calculations.
 *
 * - If `params.rr_src` is true, applies subsidence from parent class and
 *   sets `alpha` from forcing array `F`.
 * - Otherwise, `alpha` is set to zero.
 */
template <class ct_params_t>
void slvr_blk_1m_common<ct_params_t>::rr_src()
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

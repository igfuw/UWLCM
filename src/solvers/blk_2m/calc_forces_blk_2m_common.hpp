//TODO: move calc_forces and forcings to case class?
#pragma once
#include "../../forcings/subsidence.hpp"
#include "../common/calc_forces_common.hpp"
#include "../slvr_blk_2m.hpp"

/**
 * @brief Apply cloud water (rc) source term due to large-scale vertical motion.
 *
 * If `params.rc_src` is true, the subsidence is applied and the forcing coefficient
 * `alpha` is set from the large-scale forcing function `F()`. The `beta` coefficient
 * is set to zero.
 */
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

/**
 * @brief Apply cloud droplet number (nc) source term due to large-scale vertical motion.
 *
 * If `params.nc_src` is true, the subsidence is applied and the forcing coefficient
 * `alpha` is set from the large-scale forcing function `F()`. The `beta` coefficient
 * is set to zero.
 */
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

/**
 * @brief Apply rain water (rr) source term due to large-scale vertical motion.
 *
 * If `params.rr_src` is true, the subsidence is applied and the forcing coefficient
 * `alpha` is set from the large-scale forcing function `F()`. The `beta` coefficient
 * is set to zero.
 */
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

/**
 * @brief Apply rain droplet number (nr) source term due to large-scale vertical motion.
 *
 * If `params.nr_src` is true, the subsidence is applied and the forcing coefficient
 * `alpha` is set from the large-scale forcing function `F()`. The `beta` coefficient
 * is set to zero.
 */
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

/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

// headers common to uwlcm.cpp and run_hlpr.cpp

#pragma once

#include "detail/setup.hpp"

#include "opts/opts_common.hpp"
#include "detail/ct_params.hpp"
#include "solvers/common/calc_forces_common.hpp"

#if !defined(UWLCM_DISABLE_2D_LGRNGN) || !defined(UWLCM_DISABLE_3D_LGRNGN)
  #include "solvers/slvr_lgrngn.hpp"
  #include "solvers/lgrngn/diag_lgrngn.hpp"
  #include "solvers/lgrngn/hook_ante_delayed_step_lgrngn.hpp"
  #include "solvers/lgrngn/hook_ante_loop_lgrngn.hpp"
  #include "solvers/lgrngn/hook_ante_step_lgrngn.hpp"
  #include "solvers/lgrngn/hook_mixed_rhs_ante_step_lgrngn.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_BLK_1M) || !defined(UWLCM_DISABLE_3D_BLK_1M)
  #include "solvers/slvr_blk_1m.hpp"
  #include "solvers/blk_1m/calc_forces_blk_1m_common.hpp"
  #include "solvers/blk_1m/update_rhs_blk_1m_common.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_BLK_2M) || !defined(UWLCM_DISABLE_3D_BLK_2M)
  #include "solvers/slvr_blk_2m.hpp"
  #include "solvers/blk_2m/calc_forces_blk_2m_common.hpp"
  #include "solvers/blk_2m/update_rhs_blk_2m_common.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_NONE) || !defined(UWLCM_DISABLE_3D_NONE)
  #include "solvers/slvr_dry.hpp"
#endif


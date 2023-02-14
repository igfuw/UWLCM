/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "opts_common.hpp"
#include "../solvers/slvr_dry.hpp"

// simulation and output parameters for micro=lgrngn
template <class solver_t, class user_params_t, class case_ptr_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params, 
  const user_params_t &user_params,
  const case_ptr_t &case_ptr,  
  typename std::enable_if<std::is_same<
    decltype(solver_t::solver_family),
    uwlcm_dry_family_tag
  >::value>::type* = 0

)
{
}

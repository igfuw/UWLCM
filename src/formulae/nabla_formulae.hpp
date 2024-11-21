/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/common.hpp>
#include <libmpdata++/formulae/nabla_formulae.hpp>

template <class arr_1D_t, class real_t>
void grad_fwd(arr_1D_t in, arr_1D_t out, real_t dz, rng_t rng)
{     
  // extrapolate upward, top cell is two times lower
  in(rng.last() + 1) = 1.5*in(rng.last()) - .5 * in(rng.last()-1);
  out(rng) = ( in(rng+1) - in(rng)) / dz;
  // top nad bottom cells are two times lower
  out(rng.last()) *= 2;
  out(0) *= 2;
}

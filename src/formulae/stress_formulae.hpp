/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/common.hpp>
#include <libmpdata++/formulae/nabla_formulae.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace stress
    {
      using arakawa_c::h;
      using opts::opts_t;
      
      // Compact formulation
      
      // surface stress
      // 2D version
      template <int nd, opts_t opts, class arr_t, class arrvec_t, class real_t, class ijk_t, class ijkm_t>
      inline void calc_drag_cmpct_fricvel(arrvec_t &tau,
                                          const arrvec_t &v,
                                          const arr_t &rho,
                                          const real_t fricvelsq,
                                          const ijk_t &ijk,
                                          const ijkm_t &ijkm,
                                          typename std::enable_if<nd == 2>::type* = 0)
      {
        auto zro = rng_t(0, 0);

        // tau[0] as tmp storage of |U|
        tau[0](ijkm[0] + h, zro) = abs((v[0](ijkm[0] + 1, zro) + v[0](ijkm[0], zro)));

        tau[0](ijkm[0] + h, zro) = where(tau[0](ijkm[0] + h, zro) == real_t(0), real_t(0),
                                     fricvelsq / 8 / 
                                     tau[0](ijkm[0] + h, zro) *
                                     (v[0](ijkm[0] + 1, zro) + v[0](ijkm[0], zro)) *
                                     (  G<opts, 0>(rho, ijkm[0] + 1, zro)
                                      + G<opts, 0>(rho, ijkm[0]    , zro) )
                                   );
      }
  
      // 3D version
      template <int nd, opts_t opts, class arr_t, class arrvec_t, class real_t, class ijk_t, class ijkm_t>
      inline void calc_drag_cmpct_fricvel(arrvec_t &tau,
                                          const arrvec_t &v,
                                          const arr_t &rho,
                                          const real_t fricvelsq,
                                          const ijk_t &ijk,
                                          const ijkm_t &ijkm,
                                          typename std::enable_if<nd == 3>::type* = 0)
      {
        auto zro = rng_t(0, 0);

        // tau[0] as tmp storage of |U|
        tau[0](ijkm[0] + h, ijk[1], zro) = sqrt(  pow2((v[0](ijkm[0] + 1, ijk[1], zro) + v[0](ijkm[0], ijk[1], zro)))
                                                + pow2((v[1](ijkm[0] + 1, ijk[1], zro) + v[1](ijkm[0], ijk[1], zro))));

        tau[0](ijkm[0] + h, ijk[1], zro) = where(tau[0](ijkm[0] + h, ijk[1], zro) == real_t(0), real_t(0),
                                             fricvelsq / 8 / 
                                             tau[0](ijkm[0] + h, ijk[1], zro) *
                                             (v[0](ijkm[0] + 1, ijk[1], zro) + v[0](ijkm[0], ijk[1], zro)) *
                                             (  G<opts, 0>(rho, ijkm[0] + 1, ijk[1], zro)
                                              + G<opts, 0>(rho, ijkm[0]    , ijk[1], zro) )
                                           );

        // tau[1] as tmp storage of |U|
        tau[1](ijk[0], ijkm[1] + h, zro) = sqrt(  pow2((v[0](ijk[0], ijkm[1] + 1, zro) + v[0](ijk[0], ijkm[1], zro)))
                                                + pow2((v[1](ijk[0], ijkm[1] + 1, zro) + v[1](ijk[0], ijkm[1], zro)))); 

        tau[1](ijk[0], ijkm[1] + h, zro) = where(tau[1](ijk[0], ijkm[1] + h, zro) == real_t(0), real_t(0),
                                             fricvelsq / 8 / 
                                             tau[1](ijk[0], ijkm[1] + h, zro) *
                                             (v[1](ijk[0], ijkm[1] + 1, zro) + v[1](ijk[0], ijkm[1], zro)) *
                                             (  G<opts, 0>(rho, ijk[0], ijkm[1] + 1, zro)
                                              + G<opts, 0>(rho, ijk[0], ijkm[1]    , zro) )
                                           );
      }

      // zeroing drag
      template <int nd, class arrvec_t, class ijk_t, class ijkm_t>
      inline void zero_drag_cmpct(arrvec_t &tau,
                                  const ijk_t &ijk,
                                  const ijkm_t &ijkm,
                                  typename std::enable_if<nd == 2>::type* = 0)
      {
        auto zro = rng_t(0, 0);
        tau[0](ijkm[0] + h, zro) = 0;
      }

      template <int nd, class arrvec_t, class ijk_t, class ijkm_t>
      inline void zero_drag_cmpct(arrvec_t &tau,
                                  const ijk_t &ijk,
                                  const ijkm_t &ijkm,
                                  typename std::enable_if<nd == 3>::type* = 0)
      {
        auto zro = rng_t(0, 0);
        tau[0](ijkm[0] + h, ijk[1], zro) = 0;
        tau[1](ijk[0], ijkm[1] + h, zro) = 0;
      };

    } // namespace stress
  } // namespace formulae
} // namespace libmpdataxx

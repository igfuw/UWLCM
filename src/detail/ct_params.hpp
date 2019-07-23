/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/opts.hpp>
#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>

//
////#include "detail/setup.hpp"
////#include "detail/concurr_types.hpp"
////
////#include "cases/DYCOMS.hpp"
////#include "cases/MoistThermalGrabowskiClark99.hpp"
////#include "cases/DryThermalGMD2015.hpp"
////#include "cases/LasherTrapp2001.hpp"
//
//#include "opts/opts_common.hpp"
////#include "opts/opts_lgrngn.hpp"
////#include "opts/opts_blk_1m.hpp"
//
//#include "solvers/slvr_lgrngn.hpp"
//#include "solvers/slvr_blk_1m.hpp"
//
////#include "detail/panic.hpp"
//#include <map>

// libmpdata++'s compile-time parameters
struct ct_params_common : libmpdataxx::ct_params_default_t
{
  using real_t = setup::real_t;
  enum { rhs_scheme = libmpdataxx::solvers::mixed }; 
  enum { prs_scheme = libmpdataxx::solvers::cr };
  enum { vip_vab = libmpdataxx::solvers::expl };

  enum { opts = libmpdataxx::opts::nug 
#if defined(MPDATA_OPTS_IGA)
    | libmpdataxx::opts::iga 
#endif
#if defined(MPDATA_OPTS_FCT)
    | libmpdataxx::opts::fct 
#endif
#if defined(MPDATA_OPTS_ABS)
    | libmpdataxx::opts::abs
#endif
  };
};

struct ct_params_2D_sd : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 4 };
  struct ix { enum {
    u, w, th, rv, 
    vip_i=u, vip_j=w, vip_den=-1
  }; };
  enum { delayed_step = libmpdataxx::opts::bit(ix::th) | libmpdataxx::opts::bit(ix::rv) };
};

struct ct_params_3D_sd : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 5 };
  struct ix { enum {
    u, v, w, th, rv, 
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
  enum { delayed_step = libmpdataxx::opts::bit(ix::th) | libmpdataxx::opts::bit(ix::rv) };
};

struct ct_params_2D_blk_1m : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 6 };
  struct ix { enum {
    u, w, th, rv, rc, rr,
    vip_i=u, vip_j=w, vip_den=-1
  }; };
};

struct ct_params_3D_blk_1m : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 7 };
  struct ix { enum {
    u, v, w, th, rv, rc, rr, 
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
};

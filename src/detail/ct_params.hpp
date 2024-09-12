/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/opts.hpp>
#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>

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

struct ct_params_2D_lgrngn : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 4 };
  struct ix { enum {
    u, w, th, rv, 
    vip_i=u, vip_j=w, vip_den=-1
  }; };
  enum { delayed_step = libmpdataxx::opts::bit(ix::th) | libmpdataxx::opts::bit(ix::rv) };
};

struct ct_params_3D_lgrngn : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 5 };
  struct ix { enum {
    u, v, w, th, rv, 
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
  enum { delayed_step = libmpdataxx::opts::bit(ix::th) | libmpdataxx::opts::bit(ix::rv) };
};

struct ct_params_2D_lgrngn_chem : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 10 };
  struct ix { enum {
    u, w, th, rv, 
    SO2g, O3g, H2O2g, CO2g, NH3g, HNO3g,
    vip_i=u, vip_j=w, vip_den=-1
  }; };
  enum { delayed_step = libmpdataxx::opts::bit(ix::th) | libmpdataxx::opts::bit(ix::rv) 
                      | libmpdataxx::opts::bit(ix::SO2g) | libmpdataxx::opts::bit(ix::O3g)
                      | libmpdataxx::opts::bit(ix::H2O2g) | libmpdataxx::opts::bit(ix::CO2g)
                      | libmpdataxx::opts::bit(ix::NH3g) | libmpdataxx::opts::bit(ix::HNO3g)
       };
};

struct ct_params_3D_lgrngn_chem : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 11 };
  struct ix { enum {
    u, v, w, th, rv, 
    SO2g, O3g, H2O2g, CO2g, NH3g, HNO3g,
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
  enum { delayed_step = libmpdataxx::opts::bit(ix::th) | libmpdataxx::opts::bit(ix::rv) 
                      | libmpdataxx::opts::bit(ix::SO2g) | libmpdataxx::opts::bit(ix::O3g)
                      | libmpdataxx::opts::bit(ix::H2O2g) | libmpdataxx::opts::bit(ix::CO2g)
                      | libmpdataxx::opts::bit(ix::NH3g) | libmpdataxx::opts::bit(ix::HNO3g)
       };
};

struct ct_params_2D_blk_2m : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 8 };
  struct ix { enum {
    u, w, th, rv, rc, rr, nc, nr,
    vip_i=u, vip_j=w, vip_den=-1
  }; };
};

struct ct_params_3D_blk_2m : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 9 };
  struct ix { enum {
    u, v, w, th, rv, rc, rr, nc, nr, 
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
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

struct ct_params_2D_dry : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 4 };
  struct ix { enum {
    u, w, th, rv,
    vip_i=u, vip_j=w, vip_den=-1
  }; };
};

struct ct_params_3D_dry : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 5 };
  struct ix { enum {
    u, v, w, th, rv,
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
};

// ct_params type traits
template< typename, typename = void >
constexpr bool has_SO2g{};
template< typename T >
constexpr bool has_SO2g<T, std::void_t<decltype(T::SO2g)>> = true;

template< typename, typename = void >
constexpr bool has_HNO3g{};
template< typename T >
constexpr bool has_HNO3g<T, std::void_t<decltype(T::HNO3g)>> = true;

template< typename, typename = void >
constexpr bool has_NH3g{};
template< typename T >
constexpr bool has_NH3g<T, std::void_t<decltype(T::NH3g)>> = true;

template< typename, typename = void >
constexpr bool has_CO2g{};
template< typename T >
constexpr bool has_CO2g<T, std::void_t<decltype(T::CO2g)>> = true;

template< typename, typename = void >
constexpr bool has_H2O2g{};
template< typename T >
constexpr bool has_H2O2g<T, std::void_t<decltype(T::H2O2g)>> = true;

template< typename, typename = void >
constexpr bool has_O3g{};
template< typename T >
constexpr bool has_O3g<T, std::void_t<decltype(T::O3g)>> = true;

#pragma once
#include "setup.hpp"
#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/const_cp.hpp>

struct calc_c_p
{
  setup::real_t operator()(setup::real_t rv) const
  {return libcloudphxx::common::moist_air::c_p<setup::real_t>(rv) * si::kilograms * si::kelvins / si::joules;}
  BZ_DECLARE_FUNCTOR(calc_c_p)
};

/*
struct calc_T
{
  setup::real_t operator()(setup::real_t th, setup::real_t rhod) const
  {return libcloudphxx::common::theta_dry::T<setup::real_t>(th * si::kelvins, rhod * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins;}
  BZ_DECLARE_FUNCTOR2(calc_T)
};
*/

struct calc_exner
{
  setup::real_t operator()(setup::real_t p) const
  {return libcloudphxx::common::theta_std::exner<setup::real_t>(p * si::pascals);}
  BZ_DECLARE_FUNCTOR(calc_exner)
};

struct calc_l_v
{
  setup::real_t operator()(setup::real_t T) const
  {return libcloudphxx::common::const_cp::l_v<setup::real_t>(T * si::kelvins) * si::kilograms / si::joules;}
  BZ_DECLARE_FUNCTOR(calc_l_v)
};

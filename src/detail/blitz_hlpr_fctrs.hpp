#pragma once

// helper functors
struct calc_c_p
{
  setup::real_t operator()(setup::real_t rv) const
  {return libcloudphxx::common::moist_air::c_p<setup::real_t>(rv) * si::kilograms * si::kelvins / si::joules;}
  BZ_DECLARE_FUNCTOR(calc_c_p)
};

struct calc_T
{
  setup::real_t operator()(setup::real_t th, setup::real_t rhod) const
  {return libcloudphxx::common::theta_dry::T<setup::real_t>(th * si::kelvins, rhod * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins;}
  BZ_DECLARE_FUNCTOR2(calc_T)
};

struct calc_exner
{
  setup::real_t operator()(setup::real_t p) const
  {return libcloudphxx::common::theta_std::exner<setup::real_t>(p * si::pascals);}
  BZ_DECLARE_FUNCTOR(calc_exner)
};

#pragma once

#include "setup.hpp"

namespace detail
{
  // container for constants that appear in forcings, some are not needed in all cases, etc...
  // TODO: make forcing functions part of case class ?
  struct ForceParameters_t
  {
    setup::real_t q_i, heating_kappa, F_0, F_1, rho_i, D, coriolis_parameter;
  };
};


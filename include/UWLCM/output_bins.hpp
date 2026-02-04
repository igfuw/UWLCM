#pragma once

#include <vector>
#include "setup.hpp"

namespace
{
  // bin sizes for calc and plot
  std::vector<quantity<si::length, setup::real_t> > bins_dry()
  {
    std::vector<quantity<si::length, setup::real_t> > ret;
    // dry radius bins: .001 ... .01 ... 10 .. 25.119 (22 bins in total)
    for (int i = 0; i < 22; ++i)
      ret.push_back(setup::real_t(1e-6 * pow(10, -3 + i * .2)) * si::metres);
    return ret;
  }
  
  std::vector<quantity<si::length, setup::real_t> > bins_wet()
  {
    std::vector<quantity<si::length, setup::real_t> > ret;
  
    // wet radius bins: .001 ... .01 ... 1 mm (30 bins in total)
    for (int i = 0; i < 30; ++i)
      ret.push_back(setup::real_t(1e-6 * pow(10, -3 + i * .2)) * si::metres);
    return ret;
  }

  std::vector<quantity<si::length, setup::real_t> > bins_ice()
  {
    std::vector<quantity<si::length, setup::real_t> > ret;

    // ice polar and equatorial radius bins: .001 ... .01 ... 1 mm (30 bins in total)
    for (int i = 0; i < 30; ++i)
      ret.push_back(setup::real_t(1e-6 * pow(10, -3 + i * .2)) * si::metres);
    return ret;
  }

};

#pragma once

#include <vector>

namespace
{
  // bin sizes for calc and plot
  std::vector<quantity<si::length> > bins_dry()
  {
    std::vector<quantity<si::length> > ret;
    // dry radius bins: .001 ... .01 ... 10 .. 25.119 (22 bins in total)
    for (int i = 0; i < 22; ++i)
      ret.push_back(1e-6 * pow(10, -3 + i * .2) * si::metres);
    return ret;
  }
  
  std::vector<quantity<si::length> > bins_wet()
  {
    std::vector<quantity<si::length> > ret;
  
    // wet radius bins: .001 ... .01 ... 1 mm (30 bins in total)
    for (int i = 0; i < 30; ++i)
      ret.push_back(1e-6 * pow(10, -3 + i * .2) * si::metres);
    return ret;
  }
};

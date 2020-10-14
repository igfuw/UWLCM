#pragma once

namespace setup 
{
  // special case for api tests - low aerosol concentration to avoid multiplicity overflows in tests with very low nx/ny/nz
  template<class case_t>
  class api_test : public case_t
  {
    public:
    
    api_test()
    {
      this->n1_stp = real_t(2) / si::cubic_metres;
      this->n2_stp = real_t(1) / si::cubic_metres;
    }
  };
};

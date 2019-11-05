#pragma once

#include "../../detail/user_params.hpp"

namespace setup 
{
  // special case for api tests - low aerosol concentration to avoid multiplicity overflows in tests with very low nx/ny/nz
  template<class case_t>
  class api_test : public case_t
  {
    public:
    
    user_params_t user_params;

    api_test()
    {
      //this->user_params.n1_stp = real_t(2) / si::cubic_metres;
      //this->user_params.n2_stp = real_t(1) / si::cubic_metres;
      
      user_params.n1_stp = real_t(2) / si::cubic_metres;
      user_params.n2_stp = real_t(1) / si::cubic_metres;
    }
  };
};

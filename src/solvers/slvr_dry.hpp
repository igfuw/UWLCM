#pragma once
#include "slvr_sgs.hpp"

template <class ct_params_t>
class slvr_dry : public std::conditional_t<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles,
                                              slvr_common<ct_params_t>,
                                              slvr_sgs<ct_params_t>
                                             >
{
  using parent_t = std::conditional_t<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles,
                                    slvr_common<ct_params_t>,
                                    slvr_sgs<ct_params_t>
                                   >;

  bool rain_flag = false;

  public:
  using solver_family = uwlcm_dry_family_tag;


  protected:
  virtual typename parent_t::arr_t get_rc(typename parent_t::arr_t& tmp) final
  {
    return this->r_l; // r_l should be =0 in dry
  }

  void get_puddle() final
  {
    for(int i=0; i < this->n_puddle_scalars; ++i)
    {
      this->puddle[static_cast<cmn::output_t>(i)] = 0;
    }
  }

  void set_rain(bool val) final
  {
    rain_flag = val;
  }

  bool get_rain() final
  {
    return rain_flag;
  }

  void hook_mixed_rhs_ante_loop()
  {}



  public:

  using parent_t::parent_t;
};

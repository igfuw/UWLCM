#pragma once
#include "slvr_lgrngn.hpp"

template <class ct_params_t>
class slvr_lgrngn_chem : public slvr_lgrngn<ct_params_t>
{
  using parent_t = slvr_lgrngn<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  using arr_sub_t = typename parent_t::arr_sub_t;

  private:

  void diag_pH()
  void diag_chem()

  protected:

  void set_chem(bool val)
  {
    parent_t::params.cloudph_opts.chem_rct = val ? params.flag_chem : false;
  };

  void set_rain(bool val)
  {
    parent_t::set_rain(val);
    set_chem(val);
  }

  public:

  struct rt_params_t : parent_t::rt_params_t
  {
    outmom_t<real_t> out_chem, out_wet_pH;
  };

  private:

  // per-thread copy of params
  rt_params_t params;

  public:

  // ctor
  slvr_lgrngn_chem(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    params(p)
  {
  }
};

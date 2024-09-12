#pragma once
#include <boost/assign/ptr_map_inserter.hpp>
#include "slvr_lgrngn.hpp"

template <class ct_params_t>
class slvr_lgrngn_chem : public slvr_lgrngn<ct_params_t>
{
  using parent_t = slvr_lgrngn<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  using arr_sub_t = typename parent_t::arr_sub_t;
  using chem_species_t = libcloudphxx::common::chem::chem_species_t;
  using ambient_chem_t = std::map<enum libcloudphxx::common::chem::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> >;
  using ambient_chem_ct = std::map<enum libcloudphxx::common::chem::chem_species_t, const libcloudphxx::lgrngn::arrinfo_t<real_t> >;

  private:

  typename parent_t::arr_t &SO2_pre_cond,
                           &O3_pre_cond,
                           &H2O2_pre_cond,
                           &CO2_pre_cond,
                           &NH3_pre_cond,
                           &HNO3_pre_cond,

                           &SO2_post_cond,
                           &O3_post_cond,
                           &H2O2_post_cond,
                           &CO2_post_cond,
                           &NH3_post_cond,
                           &HNO3_post_cond;

  void diag_pH();
  void diag_chem();

  protected:

  void hook_ante_loop(int nt) override;
  void hook_ante_delayed_step() override;
  void hook_mixed_rhs_ante_step() override;
  void diag() override;
  void init_prtcls() override;
  void sync_in() override;
  void step_cond() override;


  void set_chem(bool val)
  {
    this->params.cloudph_opts.chem_rct = val ? params.flag_chem_rct : false;
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
    bool flag_chem_rct;
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
    params(p),

    SO2_pre_cond(args.mem->tmp[__FILE__][0][0]),
    O3_pre_cond(args.mem->tmp[__FILE__][0][1]),
    H2O2_pre_cond(args.mem->tmp[__FILE__][0][2]),
    CO2_pre_cond(args.mem->tmp[__FILE__][0][3]),
    NH3_pre_cond(args.mem->tmp[__FILE__][0][4]),
    HNO3_pre_cond(args.mem->tmp[__FILE__][0][5]),

    SO2_post_cond(args.mem->tmp[__FILE__][1][0]),
    O3_post_cond(args.mem->tmp[__FILE__][1][1]),
    H2O2_post_cond(args.mem->tmp[__FILE__][1][2]),
    CO2_post_cond(args.mem->tmp[__FILE__][1][3]),
    NH3_post_cond(args.mem->tmp[__FILE__][1][4]),
    HNO3_post_cond(args.mem->tmp[__FILE__][1][5])

  {
  }

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 6);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 6);
  }

};

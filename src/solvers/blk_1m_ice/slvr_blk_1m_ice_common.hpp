#pragma once
#include "../blk_1m/slvr_blk_1m_common.hpp"

template <class ct_params_t>
class slvr_blk_1m_ice_common : public slvr_blk_1m_common<ct_params_t>
{
  public:

  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using ix = typename ct_params_t::ix; 
  using real_t = typename ct_params_t::real_t;
  using solver_family = uwlcm_blk_1m_ice_family_tag;

  protected:
  typename parent_t::arr_t &iceA_precipitation_rate, &iceB_precipitation_rate; 

  // accumulated water falling out of domain
  real_t ice_puddle;

  void get_puddle() override
  {
    // storing puddle
    for(int i=0; i < this->n_puddle_scalars; ++i)
    {   
      this->puddle[static_cast<cmn::output_t>(i)] = (i == 8 ? parent_t::liquid_puddle : 0);
      this->puddle[static_cast<cmn::output_t>(i)] = (i == 11 ? ice_puddle : 0);
    }
  }

  void diag()
  {
    parent_t::diag();

    // recording precipitation flux
    this->record_aux_dsc("precip_rate_iceA", iceA_precipitation_rate);
    this->record_aux_dsc("precip_rate_iceB", iceB_precipitation_rate);
  } 

  void ria_src();
  void rib_src();
  //void rc_src();
  //void rr_src();
  bool get_rain() { return params.cloudph_opts.hetB; }
  void set_rain(bool val)
  { 
    params.cloudph_opts.hetB = val ? params.flag_hetB : false;
  };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    parent_t::hook_ante_loop(nt); 

    params.flag_hetB = params.cloudph_opts.hetB;

    // recording parameters
    if(this->rank==0)
    {
      this->record_aux_const("microphysics", "single-moment bulk ice");
      this->record_aux_const("cond",   "blk_1m", params.cloudph_opts.cond);  
      this->record_aux_const("cevp",   "blk_1m", params.cloudph_opts.cevp);  
      this->record_aux_const("revp",   "blk_1m", params.cloudph_opts.revp);  
      this->record_aux_const("conv",   "blk_1m", params.flag_conv);  
      this->record_aux_const("accr",   "blk_1m", params.cloudph_opts.accr);  
      this->record_aux_const("sedi",   "blk_1m", params.cloudph_opts.sedi);
      this->record_aux_const("homA1",   "blk_1m", params.cloudph_opts.homA1);
      this->record_aux_const("homA2",   "blk_1m", params.cloudph_opts.homA2);
      this->record_aux_const("hetA",   "blk_1m", params.cloudph_opts.hetA);
      this->record_aux_const("hetB",   "blk_1m", params.cloudph_opts.hetB);
      this->record_aux_const("depA",   "blk_1m", params.cloudph_opts.depA);
      this->record_aux_const("depB",   "blk_1m", params.cloudph_opts.depB);
      this->record_aux_const("rimA",   "blk_1m", params.cloudph_opts.rimA);
      this->record_aux_const("rimB",   "blk_1m", params.cloudph_opts.rimB);
      this->record_aux_const("melA",   "blk_1m", params.cloudph_opts.melA);
      this->record_aux_const("melB",   "blk_1m", params.cloudph_opts.melB);
      this->record_aux_const("r_c0",   "blk_1m", params.cloudph_opts.r_c0);  
      this->record_aux_const("k_acnv", "blk_1m", params.cloudph_opts.k_acnv);  
      this->record_aux_const("r_eps",  "blk_1m", params.cloudph_opts.r_eps);  
      this->record_aux_const("rc_src", "rt_params", params.rc_src);  
      this->record_aux_const("rr_src", "rt_params", params.rr_src);
      // this->record_aux_const("ria_src", "rt_params", params.ria_src);
      // this->record_aux_const("rib_src", "rt_params", params.rib_src);
    }
    this->mem->barrier();
  }

  void hook_ante_step()
  {

    parent_t::hook_ante_step();

    //negtozero(this->mem->advectee(ix::ria)(this->ijk), "ria after first half of rhs");
    //negtozero(this->mem->advectee(ix::rib)(this->ijk), "rib after first half of rhs");

    nancheck(this->mem->advectee(ix::ria)(this->ijk), "ria after first half of rhs");
    negcheck(this->mem->advectee(ix::ria)(this->ijk), "ria after first half of rhs");
    nancheck(this->mem->advectee(ix::rib)(this->ijk), "rib after first half of rhs");
    negcheck(this->mem->advectee(ix::rib)(this->ijk), "rib after first half of rhs");

    this->mem->barrier();
  }

  public:

  struct rt_params_t : parent_t::rt_params_t
  {
    bool flag_hetB; // do we want hetB after spinup
  };

  protected:

  // per-thread copy of params
  // TODO: but slvr_common also has a copy of it's params....
  rt_params_t params;

  public:

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // precipitation_rate for iceA and iceB
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  );

  // ctor
  slvr_blk_1m_ice_common(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    params(p),
    ice_puddle(0),
    iceA_precipitation_rate(args.mem->tmp[__FILE__][0][0]),
    iceB_precipitation_rate(args.mem->tmp[__FILE__][0][1])
  {}  
};

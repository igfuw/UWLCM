#pragma once
#include "../slvr_sgs.hpp"
#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>

template <class ct_params_t>
class slvr_blk_1m_common : public std::conditional_t<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles,
                                                     slvr_common<ct_params_t>,
                                                     slvr_sgs<ct_params_t>
                                                    >
{
  using parent_t = std::conditional_t<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles,
                                    slvr_common<ct_params_t>,
                                    slvr_sgs<ct_params_t>
                                   >;

  public:
  using ix = typename ct_params_t::ix; // TODO: it's now in solver_common - is it needed here?
  using real_t = typename ct_params_t::real_t;
  using arr_sub_t = typename parent_t::arr_sub_t;

  protected:
  typename parent_t::arr_t &precipitation_rate; 

  private:
  // a 2D/3D array with copy of the environmental total pressure of dry air 
  typename parent_t::arr_t &p_e;

  void condevap()
  {
    auto
      th   = this->state(ix::th)(this->ijk), // potential temperature
      rv   = this->state(ix::rv)(this->ijk), // water vapour mixing ratio
      rc   = this->state(ix::rc)(this->ijk), // cloud water mixing ratio
      rr   = this->state(ix::rr)(this->ijk); // rain water mixing ratio
    auto const
      rhod = (*this->mem->G)(this->ijk),
      &p_e_arg = p_e(this->ijk);

/*
    libcloudphxx::blk_1m::adj_cellwise<real_t>( 
      params.cloudph_opts, rhod, th, rv, rc, rr, this->dt
    );
    libcloudphxx::blk_1m::adj_cellwise_constp<real_t>( 
      params.cloudph_opts, rhod, p_e_arg, th, rv, rc, rr, this->dt
    );
*/
    libcloudphxx::blk_1m::adj_cellwise_nwtrph<real_t>( 
      params.cloudph_opts, p_e_arg, th, rv, rc, this->dt
    );
    this->mem->barrier();
  }

  protected:

  // accumulated water falling out of domain
  real_t liquid_puddle;

  void get_puddle() override
  {
    // storing puddle
    for(int i=0; i < this->n_puddle_scalars; ++i)
    {   
      this->puddle[static_cast<cmn::output_t>(i)] = (i == 8 ? liquid_puddle : 0);
    }
  }

  void diag()
  {
    parent_t::diag();

    // recording precipitation flux
    this->record_aux_dsc("precip_rate", precipitation_rate);
  } 

  void rc_src();
  void rr_src();
  bool get_rain() { return params.cloudph_opts.conv; }
  void set_rain(bool val) 
  { 
    params.cloudph_opts.conv = val ? params.flag_conv : false;
  };

  virtual typename parent_t::arr_t get_rc(typename parent_t::arr_t&) final
  {
    return this->state(ix::rc);
  }

  void hook_mixed_rhs_ante_loop()
  {}

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    params.flag_conv = params.cloudph_opts.conv;

    // fill with zeros
    this->state(ix::rc)(this->ijk) = 0;
    this->state(ix::rr)(this->ijk) = 0;

    // init the p_e array
    p_e(this->ijk).reindex(this->zero) = (*params.profs.p_e)(this->vert_idx);

    // deal with initial supersaturation, TODO: don't do it here (vide slvr_lgrngn)
    condevap();

    parent_t::hook_ante_loop(nt); // forcings after adjustments

    // recording parameters
    if(this->rank==0)
    {
      this->record_aux_const("microphysics", "single-moment bulk");  
      this->record_aux_const("cond",   "blk_1m", params.cloudph_opts.cond);  
      this->record_aux_const("cevp",   "blk_1m", params.cloudph_opts.cevp);  
      this->record_aux_const("revp",   "blk_1m", params.cloudph_opts.revp);  
      this->record_aux_const("conv",   "blk_1m", params.flag_conv);  
      this->record_aux_const("accr",   "blk_1m", params.cloudph_opts.accr);  
      this->record_aux_const("sedi",   "blk_1m", params.cloudph_opts.sedi);  
      this->record_aux_const("r_c0",   "blk_1m", params.cloudph_opts.r_c0);  
      this->record_aux_const("k_acnv", "blk_1m", params.cloudph_opts.k_acnv);  
      this->record_aux_const("r_eps",  "blk_1m", params.cloudph_opts.r_eps);  
      this->record_aux_const("rc_src", "rt_params", params.rc_src);  
      this->record_aux_const("rr_src", "rt_params", params.rr_src);  
    }
    this->mem->barrier();
  }

  void hook_ante_step()
  {

    parent_t::hook_ante_step();

    negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv after first half of rhs");
    negtozero(this->mem->advectee(ix::rc)(this->ijk), "rc after first half of rhs");
    negtozero(this->mem->advectee(ix::rr)(this->ijk), "rr after first half of rhs");

    condevap(); 
    nancheck(this->mem->advectee(ix::rv)(this->ijk), "rv after condevap");
    nancheck(this->mem->advectee(ix::rc)(this->ijk), "rc after condevap");
    nancheck(this->mem->advectee(ix::rr)(this->ijk), "rr after condevap");
    negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after condevap");
    negcheck(this->mem->advectee(ix::rc)(this->ijk), "rc after condevap");
    negcheck(this->mem->advectee(ix::rr)(this->ijk), "rr after condevap");
    this->mem->barrier();
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); 
  }

  public:

  struct rt_params_t : parent_t::rt_params_t
  {
    libcloudphxx::blk_1m::opts_t<real_t> cloudph_opts;
    bool flag_conv; // do we want coal after spinup
  };

  protected:

  // per-thread copy of params
  // TODO: but slvr_common also has a copy of it's params....
  rt_params_t params;

  public:

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // p_e, precipitation_rate
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  );

  // ctor
  slvr_blk_1m_common(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    params(p),
    liquid_puddle(0),
    p_e(args.mem->tmp[__FILE__][0][0]),
    precipitation_rate(args.mem->tmp[__FILE__][0][1])
  {}  
};

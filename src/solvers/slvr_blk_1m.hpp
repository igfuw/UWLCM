#pragma once
#include "slvr_sgs.hpp"

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

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
  using clock = typename parent_t::clock;

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
  void hook_mixed_rhs_ante_step()
  {
    update_rhs(this->rhs, this->dt, 0);
    this->apply_rhs(this->dt);
  }
  void hook_mixed_rhs_post_step()
  {
    update_rhs(this->rhs, this->dt, 1);
    this->apply_rhs(this->dt);
  }

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    params.flag_conv = params.cloudph_opts.conv;

    // fill with zeros
    this->state(ix::rc)(this->ijk) = 0;
    this->state(ix::rr)(this->ijk) = 0;

    // init the p_e array
    p_e(this->ijk).reindex(this->zero) = (*params.p_e)(this->vert_idx);

    // deal with initial supersaturation, TODO: don't do it here (vide slvr_lgrngn)
    condevap();

    parent_t::hook_ante_loop(nt); // forcings after adjustments

    // recording parameters
    if(this->rank==0)
    {
      this->record_aux_const("single-moment bulk microphysics", -44);  
      this->record_aux_const("cond", params.cloudph_opts.cond);  
      this->record_aux_const("cevp", params.cloudph_opts.cevp);  
      this->record_aux_const("revp", params.cloudph_opts.revp);  
      this->record_aux_const("conv", params.flag_conv);  
      this->record_aux_const("accr", params.cloudph_opts.accr);  
      this->record_aux_const("sedi", params.cloudph_opts.sedi);  
      this->record_aux_const("r_c0", params.cloudph_opts.r_c0);  
      this->record_aux_const("k_acnv", params.cloudph_opts.k_acnv);  
      this->record_aux_const("r_eps", params.cloudph_opts.r_eps);  
      this->record_aux_const("user_params rc_src", params.user_params.rc_src);  
      this->record_aux_const("user_params rr_src", params.user_params.rr_src);  
      this->record_aux_const("rt_params rc_src", params.rc_src);  
      this->record_aux_const("rt_params rr_src", params.rr_src);  
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

    // store rl for buoyancy
    //this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk);
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    // store rl for buoyancy
    this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk);

    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    // zero-out precipitation rate, will be set in columnwise
    if(at ==0)
      precipitation_rate(this->ijk) = 0;

    this->mem->barrier();
    if(this->rank == 0)
      this->tbeg = clock::now();

    // cell-wise
    // TODO: rozne cell-wise na n i n+1 ?
    if(at == 0)
    {
      auto
	dot_th = rhs.at(ix::th)(this->ijk),
	dot_rv = rhs.at(ix::rv)(this->ijk),
	dot_rc = rhs.at(ix::rc)(this->ijk),
	dot_rr = rhs.at(ix::rr)(this->ijk);
      const auto
	th   = this->state(ix::th)(this->ijk),
	rv   = this->state(ix::rv)(this->ijk),
	rc   = this->state(ix::rc)(this->ijk),
	rr   = this->state(ix::rr)(this->ijk),
        rhod = (*this->mem->G)(this->ijk),
        &p_e_arg = p_e(this->ijk);
      libcloudphxx::blk_1m::rhs_cellwise_nwtrph<real_t>(
          params.cloudph_opts,
          dot_th, dot_rv, dot_rc, dot_rr,
          rhod, p_e_arg, th, rv, rc, rr,
          dt
      );
    }

    // forcing
    switch (at)
    {
      // for eulerian integration or used to init trapezoidal integration
      case (0):
      {
        // ---- cloud water sources ----
        rc_src();
        rhs.at(ix::rc)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rc)(this->ijk);
        nancheck(rhs.at(ix::rc)(this->ijk), "RHS of rc after rc_src");

        // ---- rain water sources ----
        rr_src();
        rhs.at(ix::rr)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rr)(this->ijk);
        nancheck(rhs.at(ix::rr)(this->ijk), "RHS of rr after rr_src");

    
        // when using explicit turbulence model add subgrid forces to rc and rr
        // (th and rv were already applied in slvr_sgs)
        if (ct_params_t::sgs_scheme != libmpdataxx::solvers::iles)
        {
          this->sgs_scalar_forces({ix::rc, ix::rr});
          nancheck(rhs.at(ix::rc)(this->ijk), "RHS of rc after sgs_scalar_forces");
          nancheck(rhs.at(ix::rr)(this->ijk), "RHS of rr after sgs_scalar_forces");
        }
        
        break;
      }

      case (1):
      // trapezoidal rhs^n+1
      {
/*
        // ---- cloud water sources ----
        rc_src();
        rhs.at(ix::rc)(this->ijk) += this->alpha(this->ijk) + this->beta(this->ijk) * this->state(ix::rc)(this->ijk) / (1. - 0.5 * this->dt * this->beta(this->ijk));

        // ---- rain water sources ----
        rr_src();
        rhs.at(ix::rr)(this->ijk) += this->alpha(this->ijk) + this->beta(this->ijk) * this->state(ix::rr)(this->ijk) / (1. - 0.5 * this->dt * this->beta(this->ijk));

*/
        break;
      }
    }
    this->mem->barrier();
    if(this->rank == 0)
    {
      this->tend = clock::now();
      this->tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( this->tend - this->tbeg );
    }
    
  }

  //
  void hook_post_step()
  {
    //condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment
    parent_t::hook_post_step(); // includes the above forcings

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

template <class ct_params_t, class enableif = void>
class slvr_blk_1m
{};

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?

// 2D version
template <class ct_params_t>
class slvr_blk_1m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_1m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;
  using clock = typename parent_t::clock;

  // ctor
  slvr_blk_1m(
    typename parent_t::ctor_args_t args,
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p)
  {}

  protected:
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    this->mem->barrier();
    if(at == 0)
    {
      if(this->rank == 0)
        this->tbeg = clock::now();

      // column-wise
      for (int i = this->i.first(); i <= this->i.last(); ++i)
      {
        auto 
          precipitation_rate_arg = this->precipitation_rate(i, this->j);
        const auto 
          rhod   = (*this->mem->G)(i, this->j),
          rr     = this->state(parent_t::ix::rr)(i, this->j);
        this->liquid_puddle += - libcloudphxx::blk_1m::rhs_columnwise<real_t>(this->params.cloudph_opts, precipitation_rate_arg, rhod, rr, this->params.dz);
      }
      rhs.at(parent_t::ix::rr)(this->ijk) += this->precipitation_rate(this->ijk);


      nancheck(rhs.at(parent_t::ix::rr)(this->ijk), "RHS of rr after rhs_update");
      this->mem->barrier();
      if(this->rank == 0)
      {
        this->tend = clock::now();
        this->tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( this->tend - this->tbeg );
      }
    }
  }
};

// 3D version
template <class ct_params_t>
class slvr_blk_1m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_blk_1m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;
  using clock = typename parent_t::clock;

  // ctor
  slvr_blk_1m(
    typename parent_t::ctor_args_t args,
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p)
  {}

  protected:
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    this->mem->barrier();
    if(at == 0)
    {
      if(this->rank == 0)
        this->tbeg = clock::now();

      // column-wise
      for (int i = this->i.first(); i <= this->i.last(); ++i)
        for (int j = this->j.first(); j <= this->j.last(); ++j)
        {
          auto 
          precipitation_rate_arg = this->precipitation_rate(i, j, this->k);
          const auto 
          rhod   = (*this->mem->G)(i, j, this->k),
          rr     = this->state(parent_t::ix::rr)(i, j, this->k);
          this->liquid_puddle += - libcloudphxx::blk_1m::rhs_columnwise<real_t>(this->params.cloudph_opts, precipitation_rate_arg, rhod, rr, this->params.dz);
        }
      rhs.at(parent_t::ix::rr)(this->ijk) += this->precipitation_rate(this->ijk);

      nancheck(rhs.at(parent_t::ix::rr)(this->ijk), "RHS of rr after rhs_update");
      this->mem->barrier();
      if(this->rank == 0)
      {
        this->tend = clock::now();
        this->tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( this->tend - this->tbeg );
      }
    }
  }
};

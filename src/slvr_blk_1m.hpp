#pragma once
#include "slvr_common.hpp"

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

template <class ct_params_t>
class slvr_blk_1m_common : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix; // TODO: it's now in solver_common - is it needed here?
  using real_t = typename ct_params_t::real_t;
  private:

  void condevap()
  {
    auto 
      th   = this->state(ix::th)(this->ijk), // potential temperature
      rv   = this->state(ix::rv)(this->ijk), // water vapour mixing ratio
      rc   = this->state(ix::rc)(this->ijk), // cloud water mixing ratio
      rr   = this->state(ix::rr)(this->ijk); // rain water mixing ratio
    auto const
      rhod = (*this->mem->G)(this->ijk);
      
    libcloudphxx::blk_1m::adj_cellwise<real_t>( 
      opts, rhod, th, rv, rc, rr, this->dt
    );
    this->mem->barrier(); 
  }

  void zero_if_uninitialised(int e)
  {
    if (!finite(sum(this->state(e)(this->ijk)))) 
      this->state(e)(this->ijk) = 0;
  }

  protected:

  bool get_rain() { return opts.conv; }
  void set_rain(bool val) { opts.conv = val; };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rc);
    zero_if_uninitialised(ix::rr);

    // deal with initial supersaturation
    condevap();

    parent_t::hook_ante_loop(nt); // forcings after adjustments

    // recording parameters
    if(this->rank==0)
    {
      this->record_aux_const("single-moment bulk microphysics", -44);  
      this->record_aux_const("cond", opts.cond);  
      this->record_aux_const("cevp", opts.cevp);  
      this->record_aux_const("revp", opts.revp);  
      this->record_aux_const("conv", opts.conv);  
      this->record_aux_const("accr", opts.accr);  
      this->record_aux_const("sedi", opts.sedi);  
      this->record_aux_const("r_c0", opts.r_c0);  
      this->record_aux_const("k_acnv", opts.k_acnv);  
      this->record_aux_const("r_eps", opts.r_eps);  
    }
  }

  void hook_ante_step()
  {
    parent_t::hook_ante_step();
    // store rl for buoyancy
    this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk);
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    // cell-wise
    {
      auto 
	dot_rc = rhs.at(ix::rc)(this->ijk),
	dot_rr = rhs.at(ix::rr)(this->ijk);
      const auto 
	rc   = this->state(ix::rc)(this->ijk),
	rr   = this->state(ix::rr)(this->ijk);
      libcloudphxx::blk_1m::rhs_cellwise<real_t>(opts, dot_rc, dot_rr, rc, rr);
    }
  }

  // 
  void hook_post_step()
  {
    condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment
    parent_t::hook_post_step(); // includes the above forcings
  }

  libcloudphxx::blk_1m::opts_t<real_t> opts; // local copy of opts from rt_params, why is it needed? use rt_params::cloudph_opts instead?

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    libcloudphxx::blk_1m::opts_t<real_t> cloudph_opts;
  };

  // ctor
  slvr_blk_1m_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
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

    // column-wise
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      auto 
	dot_rr = rhs.at(parent_t::ix::rr)(i, this->j);
      const auto 
        rhod   = (*this->mem->G)(i, this->j),
	rr     = this->state(parent_t::ix::rr)(i, this->j);
      libcloudphxx::blk_1m::rhs_columnwise<real_t>(this->opts, dot_rr, rhod, rr, this->params.dz);
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

    // column-wise
    for (int i = this->i.first(); i <= this->i.last(); ++i)
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        auto 
  	dot_rr = rhs.at(parent_t::ix::rr)(i, j, this->k);
        const auto 
        rhod   = (*this->mem->G)(i, j, this->k),
  	rr     = this->state(parent_t::ix::rr)(i, j, this->k);
        libcloudphxx::blk_1m::rhs_columnwise<real_t>(this->opts, dot_rr, rhod, rr, this->params.dz);
      }
  }
};

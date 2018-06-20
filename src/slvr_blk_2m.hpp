#pragma once
#include "slvr_common.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

template <class ct_params_t>
class slvr_blk_2m_common : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix; // TODO: it's now in solver_common - is it needed here?
  using real_t = typename ct_params_t::real_t;
  private:

  void zero_if_uninitialised(int e)
  {
    if (!finite(sum(this->state(e)(this->ijk)))) 
      this->state(e)(this->ijk) = 0;
  }

  protected:

  bool get_rain() { return opts.acnv; }
  void set_rain(bool val) {
    opts.acnv = val; 
    opts.RH_max = val ? 44 : 1.01;
  };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rc);
    zero_if_uninitialised(ix::rr);
    zero_if_uninitialised(ix::nc);
    zero_if_uninitialised(ix::nr);

    parent_t::hook_ante_loop(nt); // forcings after adjustments

    // recording parameters
    if(this->rank==0)
    {
      this->record_aux_const("double-moment bulk microphysics", -44);  
      this->record_aux_const("acti", opts.acti);  
      this->record_aux_const("cond", opts.cond);  
      this->record_aux_const("accr", opts.accr);  
      this->record_aux_const("acnv", opts.acnv);  
      this->record_aux_const("sedi", opts.sedi);  
      this->record_aux_const("acnv_A", opts.acnv_A);  
      this->record_aux_const("acnv_b", opts.acnv_b);  
      this->record_aux_const("acnv_c", opts.acnv_c);  
      //TODO - how to record this?
      //this->record_aux_const("blk2m_mean_rd", opts.dry_distros.mean_rd);
      //this->record_aux_const("blk2m_sdev_rd", opts.dry_distros.sdev_rd);
      //this->record_aux_const("blk2m_N_stp", opts.dry_distros.N_stp);
      //this->record_aux_const("blk2m_chem_b", opts.dry_distros.chem_b);
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
        dot_th = rhs.at(ix::th)(this->ijk),
        dot_rv = rhs.at(ix::rv)(this->ijk),
        dot_rc = rhs.at(ix::rc)(this->ijk),
        dot_rr = rhs.at(ix::rr)(this->ijk),
        dot_nc = rhs.at(ix::nc)(this->ijk),
        dot_nr = rhs.at(ix::nr)(this->ijk),
        rc     = this->state(ix::rc)(this->ijk),
        rr     = this->state(ix::rr)(this->ijk),
        nc     = this->state(ix::nc)(this->ijk),
        nr     = this->state(ix::nr)(this->ijk);

      const auto
        rhod   = (*this->mem->G)(this->ijk),
        th     = this->state(ix::th)(this->ijk),
        rv     = this->state(ix::rv)(this->ijk);

      libcloudphxx::blk_2m::rhs_cellwise<real_t>(
        opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
        rhod,     th,     rv,     rc,     nc,     rr,     nr,
        this->dt
      );
    }
    this->mem->barrier(); // TODO: if needed, move to adv+rhs
  }

  libcloudphxx::blk_2m::opts_t<real_t> opts;

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    libcloudphxx::blk_2m::opts_t<real_t> cloudph_opts;
  };

  // ctor
  slvr_blk_2m_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  {}  
};

template <class ct_params_t, class enableif = void>
class slvr_blk_2m 
{};

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?

// 2D version 
template <class ct_params_t>
class slvr_blk_2m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_2m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_2m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_2m( 
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
	dot_rr = rhs.at(parent_t::ix::rr)(i, this->j),
	dot_nr = rhs.at(parent_t::ix::nr)(i, this->j);
      const auto 
        rhod   = (*this->mem->G)(i, this->j),
	rr     = this->state(parent_t::ix::rr)(i, this->j),
	nr     = this->state(parent_t::ix::nr)(i, this->j);
      libcloudphxx::blk_2m::rhs_columnwise<real_t>(
        this->opts, dot_rr, dot_nr,
        rhod, rr, nr,
        this->dt,
        this->params.dz);
    }
  }
};

// 3D version 
template <class ct_params_t>
class slvr_blk_2m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_blk_2m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_2m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_2m( 
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
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {
        auto 
  	dot_rr = rhs.at(parent_t::ix::rr)(i, j, this->k),
  	dot_nr = rhs.at(parent_t::ix::nr)(i, j, this->k);
        const auto 
          rhod   = (*this->mem->G)(i, j, this->k),
  	rr     = this->state(parent_t::ix::rr)(i, j, this->k),
  	nr     = this->state(parent_t::ix::nr)(i, j, this->k);
        libcloudphxx::blk_2m::rhs_columnwise<real_t>(
          this->opts, dot_rr, dot_nr,
          rhod, rr, nr,
          this->dt,
          this->params.dz);
      }
    }
  }
};

#pragma once
#include "slvr_sgs.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>


template <class ct_params_t>
class slvr_blk_2m_common : public std::conditional_t<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles,
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
  using clock = typename parent_t::clock;
  private:

  // a 2D/3D array with copy of the environmental total pressure of dry air
  typename parent_t::arr_t &p_e;

  protected:

  // accumulated water falling out of domain
  real_t liquid_puddle;

  void get_puddle() override
  {
    //storing puddle
    for(int i=0; i < this->n_puddle_scalars; ++i)
    {
      this->puddle[static_cast<cmn::output_t>(i)] = (i==8 ? liquid_puddle : 0);
    }
  }

  void diag()
  {
    parent_t::diag();
    
    // TODO: recording precipitation flux
    // this->record_aux_dsc("precip_rate",precipitation_rate);    

    /*
    assert(this->rank == 0);
    parent_t::tbeg = parent_t::clock::now();

    // recording puddle
    for(int i=0; i < 10; ++i)
    {
       this->f_puddle << i << " " << (i == 8 ? this->puddle : 0) << "\n";
    }
    this->f_puddle << "\n";

    parent_t::tend = parent_t::clock::now();
    parent_t::tdiag += std::chrono::duration_cast<std::chrono::milliseconds>( parent_t::tend - parent_t::tbeg );
    */
  }

  void rc_src();
  void nc_src();
  void rr_src();
  void nr_src();

  bool get_rain() { return params.cloudph_opts.acnv; }
  void set_rain(bool val) {
    params.cloudph_opts.acnv = val ? params.flag_acnv : false;
    params.cloudph_opts.RH_max = val ? 44 : 1.01;
  };

  virtual typename parent_t::arr_t get_rc(typename parent_t::arr_t&) final
    {
  return this->state(ix::rc);
    }

  void hook_mixed_rhs_ante_loop()
  {}
  void hook_mixed_rhs_ante_step()
  {
    const auto nc     = this->state(ix::nc)(this->ijk);
    update_rhs(this->rhs, this->dt, 0);
    negcheck(nc, "nc before apply rhs ante step");
    this->apply_rhs(this->dt);
  }
  void hook_mixed_rhs_post_step()
  {
    const auto nc     = this->state(ix::nc)(this->ijk);
    update_rhs(this->rhs, this->dt, 1);
    this->apply_rhs(this->dt);
  }

  void hook_ante_loop(int nt)
  {
    params.flag_acnv = params.cloudph_opts.acnv;

    // fill with zeros
    this->state(ix::rc)(this->ijk) = 0;
    this->state(ix::rr)(this->ijk) = 0;
    this->state(ix::nc)(this->ijk) = 0;
    this->state(ix::nr)(this->ijk) = 0;

    // init the p_e array
    p_e(this->ijk).reindex(this->zero) = (*params.p_e)(this->vert_idx);

    parent_t::hook_ante_loop(nt); // forcings after adjustments

    // recording parameters
    if(this->rank==0)
    {
      this->record_aux_const("double-moment bulk microphysics", -44);
      this->record_aux_const("acti", params.cloudph_opts.acti);
      this->record_aux_const("cond", params.cloudph_opts.cond);
      this->record_aux_const("accr", params.cloudph_opts.accr);
      this->record_aux_const("acnv", params.flag_acnv);
      this->record_aux_const("sedi", params.cloudph_opts.sedi);
      this->record_aux_const("acnv_A", params.cloudph_opts.acnv_A);
      this->record_aux_const("acnv_b", params.cloudph_opts.acnv_b);
      this->record_aux_const("acnv_c", params.cloudph_opts.acnv_c);
      this->record_aux_const("user_params rc_src", params.user_params.rc_src);
      this->record_aux_const("user_params rr_src", params.user_params.rr_src);
      this->record_aux_const("user_params nc_src", params.user_params.nc_src);
      this->record_aux_const("user_params nr_src", params.user_params.nr_src);
      this->record_aux_const("rt_params rc_src", params.rc_src);
      this->record_aux_const("rt_params rr_src", params.rr_src);
      this->record_aux_const("rt_params nc_src", params.nc_src);
      this->record_aux_const("rt_params nr_src", params.nr_src);
      /* gives error
      this->record_aux_const("blk2m_mean_rd", params.cloudph_opts.dry_distros.mean_rd);
      this->record_aux_const("blk2m_sdev_rd", params.cloudph_opts.dry_distros.sdev_rd);
      this->record_aux_const("blk2m_N_stp", params.cloudph_opts.dry_distros.N_stp);
      this->record_aux_const("blk2m_chem_b", params.cloudph_opts.dry_distros.chem_b);
      */
    }
  }


  void hook_ante_step()
  {
    parent_t::hook_ante_step();

    negtozero(this->mem->advectee(ix::rv)(this->ijk), "rv after first half of rhs");
    negtozero(this->mem->advectee(ix::rc)(this->ijk), "rc after first half of rhs");
    negtozero(this->mem->advectee(ix::rr)(this->ijk), "rr after first half of rhs");
    negtozero(this->mem->advectee(ix::nc)(this->ijk), "nc after first half of rhs");
    negtozero(this->mem->advectee(ix::nr)(this->ijk), "nr after first half of rhs");
    this->mem->barrier();
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    // store rl for buoyancy
    this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk);

    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

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
        dot_rr = rhs.at(ix::rr)(this->ijk),
        dot_nc = rhs.at(ix::nc)(this->ijk),
        dot_nr = rhs.at(ix::nr)(this->ijk);
       const auto
        rc     = this->state(ix::rc)(this->ijk),
        rr     = this->state(ix::rr)(this->ijk),
        nc     = this->state(ix::nc)(this->ijk),
        nr     = this->state(ix::nr)(this->ijk),
        rhod   = (*this->mem->G)(this->ijk),
        th     = this->state(ix::th)(this->ijk),
        rv     = this->state(ix::rv)(this->ijk),
        &p_e_arg = p_e(this->ijk); //TODO: use const pressure in blk_2m
        nancheck(nc, "nc before blk_2m rhs_cellwise call");
        negtozero(this->state(ix::nc)(this->ijk), "nc before blk_2m rhs_cellwise call");
        negcheck(nc, "nc before blk_2m rhs_cellwise call");
        nancheck(nr, "nr before blk_2m rhs_cellwise call");
        negtozero(this->state(ix::nr)(this->ijk), "nr before blk_2m rhs_cellwise call");
        negcheck(nr, "nr before blk_2m rhs_cellwise call");
       libcloudphxx::blk_2m::rhs_cellwise<real_t>(
        params.cloudph_opts, dot_th,  dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
        rhod, th,     rv,     rc,     nc,     rr,     nr,
        this->dt, true, p_e_arg
      );
        nancheck(nc, "nc after blk_2m rhs_cellwise call");
        negcheck(nc, "nc after blk_2m rhs_cellwise call");
        nancheck(nc, "nr after blk_2m rhs_cellwise call");
        negcheck(nc, "nr after blk_2m rhs_cellwise call");
    }

    // forcing
    switch (at)
    {
      // for eulerian integration or used to init trapezoidal integration
      case (0):
      {
        // ---- cloud water sources ----
        rc_src();
        //this->common_water_src(ix::rc, params.rc_src);
        rhs.at(ix::rc)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rc)(this->ijk);
        
        nc_src();
        //this->common_water_src(ix::nc, params.nc_src);
        rhs.at(ix::nc)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::nc)(this->ijk);

        // ---- rain water sources ----
        rr_src();
        //this->common_water_src(ix::rr, params.rr_src);
        rhs.at(ix::rr)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rr)(this->ijk);
        
        nr_src();
        //this->common_water_src(ix::nr, params.nr_src);
        rhs.at(ix::nr)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::nr)(this->ijk);

        // when using explicit turbulence model add subgrid forces to rc and rr
        // (th and rv were already applied in slvr_sgs)
        if (ct_params_t::sgs_scheme != libmpdataxx::solvers::iles)
        {
          this->sgs_scalar_forces({ix::rc, ix::rr, ix::nc, ix::nr});
          nancheck(rhs.at(ix::rc)(this->ijk), "RHS of rc after sgs_scalar_forces");
          nancheck(rhs.at(ix::rr)(this->ijk), "RHS of rr after sgs_scalar_forces");
          nancheck(rhs.at(ix::nc)(this->ijk), "RHS of nc after sgs_scalar_forces");
          nancheck(rhs.at(ix::nr)(this->ijk), "RHS of nr after sgs_scalar_forces");
        }

        break;
      }

      case (1):
      // trapezoidal rhs^n+1
      {
        break;
      }
    }
    this->mem->barrier();
    if(this->rank == 0)
    {
      nancheck(rhs.at(ix::rc)(this->domain), "RHS of rc after rhs_update");
      nancheck(rhs.at(ix::rr)(this->domain), "RHS of rr after rhs_update");
      nancheck(rhs.at(ix::nc)(this->domain), "RHS of nc after rhs_update");
      nancheck(rhs.at(ix::nr)(this->domain), "RHS of nr after rhs_update");
      this->tend = clock::now();
      this->tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( this->tend - this->tbeg );
    }
  }

  public:

  struct rt_params_t : parent_t::rt_params_t
  {
    libcloudphxx::blk_2m::opts_t<real_t> cloudph_opts;
    bool flag_acnv; // do we want autoconversion after spinup
  };

  protected:

  // per-thread copy of params
  // TODO: but slvr_common also has a copy of it's params....
  rt_params_t params;

  public:

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 1); // p_e
  }

  // ctor
  slvr_blk_2m_common(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    params(p),
    liquid_puddle(0),
    p_e(args.mem->tmp[__FILE__][0][0])
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
  using clock = typename parent_t::clock;

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
    parent_t::update_rhs(rhs, dt, at);

    this->mem->barrier();
    if(at == 0)
    {
      if(this->rank == 0)
        this->tbeg = clock::now();

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
        this->liquid_puddle += -libcloudphxx::blk_2m::rhs_columnwise<real_t>(
          this->params.cloudph_opts, dot_rr, dot_nr, rhod, rr, nr, this->dt, this->params.dz
        );
      }

      this->mem->barrier();
      if(this->rank == 0)
      {
        nancheck(rhs.at(parent_t::ix::rr)(this->domain), "RHS of rr after rhs_update");
        nancheck(rhs.at(parent_t::ix::nr)(this->domain), "RHS of nr after rhs_update");
        this->tend = clock::now();
        this->tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( this->tend - this->tbeg );
      }
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
  using clock = typename parent_t::clock;

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
    parent_t::update_rhs(rhs, dt, at);

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
          dot_rr = rhs.at(parent_t::ix::rr)(i, j, this->k),
          dot_nr = rhs.at(parent_t::ix::nr)(i, j, this->k);
          const auto
          rhod   = (*this->mem->G)(i, j, this->k),
          rr     = this->state(parent_t::ix::rr)(i, j, this->k),
          nr     = this->state(parent_t::ix::nr)(i, j, this->k);
          this->liquid_puddle += -libcloudphxx::blk_2m::rhs_columnwise<real_t>(
            this->params.cloudph_opts, dot_rr, dot_nr, rhod, rr, nr, this->dt, this->params.dz
          );
        }

      this->mem->barrier();
      if(this->rank == 0)
      {
        nancheck(rhs.at(parent_t::ix::rr)(this->domain), "RHS of rr after rhs_update");
        nancheck(rhs.at(parent_t::ix::nr)(this->domain), "RHS of nr after rhs_update");
        this->tend = clock::now();
        this->tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( this->tend - this->tbeg );
      }
    }
  }
};

#pragma once
#include "../slvr_sgs.hpp"
#include <libcloudph++/blk_2m/options.hpp>

/**
 * @brief Common base class for 2-moment bulk microphysics solver.
 *
 * Double-moment cloud microphysics,
 * including handling of rain and cloud water mixing ratios, fluxes. Inherits from either slvr_common
 * or slvr_sgs depending on the SGS scheme.
 *
 * @tparam ct_params_t Compile-time parameter struct defining solver options.
 */
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
  
  protected:
  typename parent_t::arr_t &rr_flux;
  typename parent_t::arr_t &nr_flux;

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
    
    // recording precipitation flux
    this->record_aux_dsc("precip_rate_rr", rr_flux);    
    this->record_aux_dsc("precip_rate_nr", nr_flux);    
  }

  void rc_src(); ///< Cloud water source term
  void nc_src(); ///< Cloud droplet concentration source term
  void rr_src(); ///< Rain source term
  void nr_src(); ///< Rain drop concentration source term
  bool get_rain() { return params.cloudph_opts.acnv; }
  void set_rain(bool val) 
  {
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
    negcheck(this->mem->advectee(ix::nc)(this->ijk), "nc at start of hook_mixed_rhs_ante_step");
  //  negtozero(this->mem->advectee(ix::nc)(this->ijk), "nc at start of hook_mixed_rhs_ante_step");
    parent_t::hook_mixed_rhs_ante_step();

    nancheck(this->mem->advectee(ix::nc)(this->ijk), "nc after mixed_rhs_ante_step apply rhs");
    nancheck(this->mem->advectee(ix::nr)(this->ijk), "nr after mixed_rhs_ante_step apply rhs");
    nancheck(this->mem->advectee(ix::rc)(this->ijk), "rc after mixed_rhs_ante_step apply rhs");
    nancheck(this->mem->advectee(ix::rr)(this->ijk), "rr after mixed_rhs_ante_step apply rhs");
    negcheck2(this->mem->advectee(ix::nc)(this->ijk), this->rhs.at(ix::nc)(this->ijk), "nc after mixed_rhs_ante_step apply rhs (+ output of nc rhs)");
//    negcheck2(this->mem->advectee(ix::nr)(this->ijk), this->rhs.at(ix::nr)(this->ijk), "nr after mixed_rhs_ante_step apply rhs (+ output of nr rhs)");
    negcheck2(this->mem->advectee(ix::rc)(this->ijk), this->rhs.at(ix::rc)(this->ijk), "rc after mixed_rhs_ante_step apply rhs (+ output of rc rhs)");
    //negcheck2(this->mem->advectee(ix::rr)(this->ijk), this->rhs.at(ix::rr)(this->ijk), "rr after mixed_rhs_ante_step apply rhs (+ output of rr rhs)");
    negtozero(this->mem->advectee(ix::nr)(this->ijk), "nr after mixed_rhs_ante_step apply rhs");
    negtozero(this->mem->advectee(ix::rr)(this->ijk), "rr after mixed_rhs_ante_step apply rhs");
  }

  void hook_mixed_rhs_post_step()
  {
    negtozero(this->mem->advectee(ix::nc)(this->ijk), "nc at start of hook_mixed_rhs_ante_step");
    negtozero(this->mem->advectee(ix::nr)(this->ijk), "nr at start of hook_mixed_rhs_ante_step");
    negtozero(this->mem->advectee(ix::rr)(this->ijk), "rr at start of hook_mixed_rhs_ante_step");
    //negcheck(this->mem->advectee(ix::nc)(this->ijk), "nc at start of mixed_rhs_post_step");

    parent_t::hook_mixed_rhs_post_step();

    nancheck(this->mem->advectee(ix::nc)(this->ijk), "nc after mixed_rhs_post_step apply rhs");
    nancheck(this->mem->advectee(ix::nr)(this->ijk), "nr after mixed_rhs_post_step apply rhs");
    nancheck(this->mem->advectee(ix::rc)(this->ijk), "rc after mixed_rhs_post_step apply rhs");
    nancheck(this->mem->advectee(ix::rr)(this->ijk), "rr after mixed_rhs_post_step apply rhs");
    negcheck2(this->mem->advectee(ix::nc)(this->ijk), this->rhs.at(ix::nc)(this->ijk), "nc after mixed_rhs_post_step apply rhs (+ output of nc rhs)");
    negcheck2(this->mem->advectee(ix::nr)(this->ijk), this->rhs.at(ix::nr)(this->ijk), "nr after mixed_rhs_post_step apply rhs (+ output of nr rhs)");
    negcheck2(this->mem->advectee(ix::rc)(this->ijk), this->rhs.at(ix::rc)(this->ijk), "rc after mixed_rhs_post_step apply rhs (+ output of rc rhs)");
    negcheck2(this->mem->advectee(ix::rr)(this->ijk), this->rhs.at(ix::rr)(this->ijk), "rr after mixed_rhs_post_step apply rhs (+ output of rr rhs)");
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
      this->record_aux_const("microphysics", "double-moment bulk");
      this->record_aux_const("acti", "blk_2m", params.cloudph_opts.acti);
      this->record_aux_const("cond", "blk_2m", params.cloudph_opts.cond);
      this->record_aux_const("accr", "blk_2m", params.cloudph_opts.accr);
      this->record_aux_const("acnv", "blk_2m", params.flag_acnv);
      this->record_aux_const("sedi", "blk_2m", params.cloudph_opts.sedi);
      this->record_aux_const("acnv_A", "blk_2m", params.cloudph_opts.acnv_A);
      this->record_aux_const("acnv_b", "blk_2m", params.cloudph_opts.acnv_b);
      this->record_aux_const("acnv_c", "blk_2m", params.cloudph_opts.acnv_c);
      this->record_aux_const("rc_src", "rt_params", params.rc_src);
      this->record_aux_const("rr_src", "rt_params", params.rr_src);
      this->record_aux_const("nc_src", "rt_params", params.nc_src);
      this->record_aux_const("nr_src", "rt_params", params.nr_src);
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


  void hook_post_step()
  {
    parent_t::hook_post_step(); 
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
    parent_t::alloc_tmp_sclr(mem, __FILE__, 3); // p_e, rr_flux, nr_flux
  }

  /**
 * @brief Update the right-hand-side of the prognostic equations.
 * @param rhs RHS arrays for all scalars
 * @param dt Time step
 * @param at Current step index
 */
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  );

  /**
 * @brief Constructor.
 * @param args Constructor arguments for parent class
 * @param p Runtime parameters
 */
  slvr_blk_2m_common(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    params(p),
    liquid_puddle(0),
    p_e(args.mem->tmp[__FILE__][0][0]),
    rr_flux(args.mem->tmp[__FILE__][0][1]),
    nr_flux(args.mem->tmp[__FILE__][0][2])
  {}
};

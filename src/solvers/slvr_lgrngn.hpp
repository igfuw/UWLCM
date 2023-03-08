#pragma once
#include "slvr_sgs.hpp"
#include "../detail/outmom.hpp"
#include <libcloudph++/lgrngn/factory.hpp>
#include <libmpdata++/formulae/refined_grid.hpp>


#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
class slvr_lgrngn : public std::conditional_t<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles,
                                              slvr_common<ct_params_t>,
                                              slvr_sgs<ct_params_t>
                                             >
{
  using parent_t = std::conditional_t<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles,
                                    slvr_common<ct_params_t>,
                                    slvr_sgs<ct_params_t>
                                   >;

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  using arr_sub_t = typename parent_t::arr_sub_t;

  private:

#if defined(UWLCM_TIMING)
  setup::clock::time_point tbeg, tend;
#endif

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  // helpers for calculating RHS from condensation, probably some of the could be avoided e.g. if step_cond returnd deltas and not changed fields 
  // or if change in theta was calculated from change in rv  
  typename parent_t::arr_t &rv_pre_cond,
                           &rv_post_cond,
                           &th_pre_cond,
                           &th_post_cond,
                           &dth,
                           &drv,
                           &r_c,  // temp storate for r_c to be used in SMG, separate storage for it allows more concurrency (like r_l)
                           &tmp_ref;

  libmpdataxx::arrvec_t<typename parent_t::arr_t> courants, // courant number on the refined grid
                                                  uvw_ref,
						  negref_dbg_arrs; // arrays to be outputted when debugging negative values in refinement

  const std::vector<std::string> negref_dbg_arr_names;

//  void diag_rl()
//  {
//    // fill with rl values from superdroplets
//    if(this->rank == 0) 
//    {
//      prtcls->diag_all();
//      prtcls->diag_wet_mom(3);
//      auto r_l_indomain = this->r_l(this->domain); // rl refrences subdomain of r_l
//      r_l_indomain = typename parent_t::arr_t(prtcls->outbuf(), r_l_indomain.shape(), blitz::duplicateData); // copy in data from outbuf; total liquid third moment of wet radius per kg of dry air [m^3 / kg]
//    }
//    this->mem->barrier();
//
//    nancheck(this->r_l(this->ijk), "rl after copying from diag_wet_mom(3)");
//    this->r_l(this->ijk) *= 4./3. * 1000. * 3.14159; // get mixing ratio [kg/kg]
//    this->mem->barrier();
//
//    // average values of rl in edge cells
//    this->avg_edge_sclr(this->r_l, this->ijk); // in case of cyclic bcond, rl on edges needs to be the same
//  }

  void diag_rx(typename parent_t::arr_t &rx)
  {
    auto &rx_ref(this->tmp_ref);
    if(this->rank == 0) 
    {
      prtcls->diag_wet_mom(3);
//      auto rc = r_c(this->domain);
      rx_ref(this->domain_ref) = typename parent_t::arr_t(prtcls->outbuf(), this->shape(this->domain_ref), blitz::neverDeleteData);
 //     std::cerr << "rx_ref(this->domain_ref) after copy from libcloud: " << rx_ref(this->domain_ref) << std::endl;
 //     std::cerr << "rx_ref after copy from libcloud: " << rx_ref << std::endl;
    }
    this->mem->barrier();


    this->xchng_ref(rx_ref, this->ijk_ref);
 //   std::cerr << "rx_ref(ijk_ref) after xchng_ref: " << rx_ref(this->ijk_ref) << std::endl;
 //   std::cerr << "rx_ref after xchng_ref: " << rx_ref << std::endl;
    libmpdataxx::formulae::refined::spatial_average_ref2reg<real_t>(rx_ref, this->ijk_r2r, this->mem->n_ref/2, this->mem->distmem.grid_size_ref, true);
 //   std::cerr << "rx_ref(ijk_ref) after spatial average: " << rx_ref(this->ijk_ref) << std::endl;
    rx(this->ijk) = rx_ref(this->ijk_r2r);
   // std::cerr << "rx(ijk) after copy from refined: " << rx(this->ijk) << std::endl;

    nancheck(rx(this->ijk), "r_c after copying from diag_wet_mom(3) in diag_rc");
    rx(this->ijk) *= 4./3. * 1000. * 3.14159; // get mixing ratio [kg/kg]
    this->mem->barrier();

    this->avg_edge_sclr(rx, this->ijk); // in case of cyclic bcond, rc on edges needs to be the same
  }

  void diag_rc()
  {
    if(this->rank == 0) 
    {
      prtcls->diag_wet_rng(.5e-6, 25.e-6);
    }
    this->mem->barrier();
    diag_rx(r_c);
  }

  void diag_rl()
  {
    if(this->rank == 0) 
    {
      prtcls->diag_all();
    }
    this->mem->barrier();
    diag_rx(this->r_l);
  }

  void get_puddle() override
  {
    this->puddle = prtcls->diag_puddle();
  }

  void diag();

  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.data(), 
      arr.stride().data()
    );
  }

  std::string aux_name(
    const std::string pfx, 
    const int rng,
    const int mom
  )
  { 
    std::ostringstream tmp;
    tmp << pfx << "_rng" << std::setw(3) << std::setfill('0') << rng << "_mom" << mom;
    return tmp.str();
  }

  protected:

  bool get_rain() { return params.cloudph_opts.coal; }
  void set_rain(bool val) 
  { 
    params.cloudph_opts.coal = val ? params.flag_coal : false;
    params.cloudph_opts.RH_max = val ? 44 : 1.01; // TODO: specify it somewhere else, dup in blk_2m
  };
  
  
  virtual typename parent_t::arr_t get_rc(typename parent_t::arr_t& tmp) final
  {
    return r_c;
  }

  void hook_ante_loop(int nt);
  void hook_ante_step();
  void hook_ante_delayed_step();
  void hook_mixed_rhs_ante_step();

  void hook_mixed_rhs_ante_loop()
  {
    diag_rl(); // init r_l
    if(ct_params_t::sgs_scheme == libmpdataxx::solvers::smg)
      diag_rc(); // ditto
  } 

#if defined(STD_FUTURE_WORKS)
  #if defined(UWLCM_TIMING)
    std::future<setup::timer> ftr;
  #else
    std::future<void> ftr;
  #endif
#endif
  
  void record_all()
  {
    assert(this->rank == 0);

#if defined(UWLCM_TIMING)
        tbeg = setup::clock::now();
#endif
#if defined(STD_FUTURE_WORKS)
    if (this->timestep > 0 && params.async)
    {
      assert(ftr.valid());
#if defined(UWLCM_TIMING)
      parent_t::tasync_gpu += ftr.get();
#else
      ftr.get();
#endif
    }
#endif
#if defined(UWLCM_TIMING)
        tend = setup::clock::now();
        parent_t::tasync_wait_in_record_all += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
#endif
    parent_t::record_all();
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    libcloudphxx::lgrngn::backend_t backend = libcloudphxx::lgrngn::undefined;
    bool async = true;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    libcloudphxx::lgrngn::opts_init_t<real_t> cloudph_opts_init;
    outmom_t<real_t> out_dry, out_wet;
    bool flag_coal; // do we want coal after spinup
    real_t gccn; // multiplicity of gccn
  };

  private:

  // per-thread copy of params
  // TODO: but slvr_common also has a copy of it's params....
  rt_params_t params;

  public:

  // ctor
  slvr_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p),
    rv_pre_cond(args.mem->tmp[__FILE__][0][0]),
    rv_post_cond(args.mem->tmp[__FILE__][0][1]),
    th_pre_cond(args.mem->tmp[__FILE__][0][2]),
    th_post_cond(args.mem->tmp[__FILE__][0][3]),
    tmp_ref(args.mem->tmp[__FILE__][0][4]),
    r_c(args.mem->tmp[__FILE__][1][0]),
    dth(args.mem->tmp[__FILE__][2][0]),
    drv(args.mem->tmp[__FILE__][2][1]),
    courants(args.mem->tmp[__FILE__][3]),
    negref_dbg_arr_names({"c_j", "d_j", "f_j", "rev"})
  {
    r_c = 0.;
    uvw_ref.push_back(this->mem->never_delete(&this->mem->psi_ref.at(this->ix_r2r.at(ix::u))));
    uvw_ref.push_back(this->mem->never_delete(&this->mem->psi_ref.at(this->ix_r2r.at(ix::v))));
    uvw_ref.push_back(this->mem->never_delete(&this->mem->psi_ref.at(this->ix_r2r.at(ix::w))));

    negref_dbg_arrs.push_back(this->mem->never_delete(&this->c_j));
    negref_dbg_arrs.push_back(this->mem->never_delete(&this->d_j));
    negref_dbg_arrs.push_back(this->mem->never_delete(&this->f_j));
    negref_dbg_arrs.push_back(this->mem->never_delete(&this->mem->psi.at(ix::rv)[0]));
    
    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr_ref(mem, __FILE__, 5); // th and rv pre and post cond, tmp_ref
    parent_t::alloc_tmp_sclr(mem, __FILE__, 1);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2);
    parent_t::alloc_tmp_vctr_ref(mem, __FILE__); // courants (refined)
  }

};

#pragma once

#include "slvr_common.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
class slvr_lgrngn : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>; 

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  blitz::Array<real_t, parent_t::n_dims-1> k_i; // TODO: make it's size in x direction smaller to match thread's domain

  // global arrays, shared among threads, TODO: in fact no need to share them?
  typename parent_t::arr_t &tmp1,
                           &r_l,
                           &F,       // forcings helper
                           &alpha,   // 'explicit' rhs part - does not depend on the value at n+1
                           &beta;    // 'implicit' rhs part - coefficient of the value at n+1
  // helper methods
  void diag()
  {
    assert(this->rank == 0);

    // recording super-droplet concentration per grid cell 
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());

    // recording relative humidity
    prtcls->diag_RH();
    this->record_aux("RH", prtcls->outbuf());

    // recording precipitation rate per grid cel
    prtcls->diag_all();
    prtcls->diag_precip_rate();
    this->record_aux("precip_rate", prtcls->outbuf());

    // recording 1st mom of rw of gccns
    prtcls->diag_dry_rng(2e-6, 1);
    prtcls->diag_wet_mom(1);
    this->record_aux("gccn_rw_mom1", prtcls->outbuf());

    // recording 0th mom of rw of gccns
    prtcls->diag_dry_rng(2e-6, 1);
    prtcls->diag_wet_mom(0);
    this->record_aux("gccn_rw_mom0", prtcls->outbuf());

    // recording 1st mom of rw of non-gccns
    prtcls->diag_dry_rng(0., 2e-6);
    prtcls->diag_wet_mom(1);
    this->record_aux("non_gccn_rw_mom1", prtcls->outbuf());

    // recording 0th mom of rw of gccns
    prtcls->diag_dry_rng(0., 2e-6);
    prtcls->diag_wet_mom(0);
    this->record_aux("non_gccn_rw_mom0", prtcls->outbuf());

    // recording 1st mom of rd of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_dry_mom(1);
    this->record_aux("actrw_rd_mom1", prtcls->outbuf());

    // recording 0th mom of rd of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_dry_mom(0);
    this->record_aux("actrw_rd_mom0", prtcls->outbuf());
   
    // recording 1st mom of rd of activated drops
    prtcls->diag_RH_ge_Sc();
    prtcls->diag_dry_mom(1);
    this->record_aux("actRH_rd_mom1", prtcls->outbuf());

    // recording 0th mom of rd of activated drops
    prtcls->diag_RH_ge_Sc();
    prtcls->diag_dry_mom(0);
    this->record_aux("actRH_rd_mom0", prtcls->outbuf());
   
    // recording requested statistical moments
    {
      // dry
      int rng_num = 0;
      for (auto &rng_moms : params.out_dry)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_dry_mom(mom);
          this->record_aux(aux_name("rd", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
    {
      // wet
      int rng_num = 0;
      for (auto &rng_moms : params.out_wet)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_wet_mom(mom);
          this->record_aux(aux_name("rw", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
  } 

  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(), 
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
    params.cloudph_opts.RH_max = val ? 44 : 1.06; // 0.5% limit during spinup // TODO: specify it somewhere else, dup in blk_2m
  };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    params.flag_coal = params.cloudph_opts.coal;

    parent_t::hook_ante_loop(nt); 

    // TODO: barrier?
    if (this->rank == 0) 
    {
      assert(params.backend != -1);
      assert(params.dt != 0); 

      // async does not make sense without CUDA
      if (params.backend != libcloudphxx::lgrngn::CUDA && params.backend != libcloudphxx::lgrngn::multi_CUDA) params.async = false;

      params.cloudph_opts_init.dt = params.dt; // advection timestep = microphysics timestep

      params.cloudph_opts_init.nx = this->mem->grid_size[0].length();
      params.cloudph_opts_init.dx = this->di;
      params.cloudph_opts_init.x0 = this->di / 2;
      params.cloudph_opts_init.x1 = (params.cloudph_opts_init.nx - .5) * this->di;

      if(parent_t::n_dims == 2) // 2D
      {
        params.cloudph_opts_init.nz = this->mem->grid_size[1].length();
        params.cloudph_opts_init.dz = this->dj;
        params.cloudph_opts_init.z0 = this->dj / 2;
        params.cloudph_opts_init.z1 = (params.cloudph_opts_init.nz - .5) * this->dj;

        params.cloudph_opts_init.n_sd_max = params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * params.cloudph_opts_init.sd_conc;
        if(params.backend == libcloudphxx::lgrngn::multi_CUDA)
          params.cloudph_opts_init.n_sd_max *= 1.5; // more space for copied SDs
      }
      else // 3D
      {
        params.cloudph_opts_init.ny = this->mem->grid_size[1].length();
        params.cloudph_opts_init.dy = this->dj;
        params.cloudph_opts_init.y0 = this->dj / 2;
        params.cloudph_opts_init.y1 = (params.cloudph_opts_init.ny - .5) * this->dj;

        params.cloudph_opts_init.nz = this->mem->grid_size[2].length();
        params.cloudph_opts_init.dz = this->dk;
        params.cloudph_opts_init.z0 = this->dk / 2;
        params.cloudph_opts_init.z1 = (params.cloudph_opts_init.nz - .5) * this->dk;

        params.cloudph_opts_init.n_sd_max = params.cloudph_opts_init.nx * params.cloudph_opts_init.ny * params.cloudph_opts_init.nz * params.cloudph_opts_init.sd_conc;
        if(params.backend == libcloudphxx::lgrngn::multi_CUDA)
          params.cloudph_opts_init.n_sd_max *= 1.5; // more space for copied SDs
      }

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
        (libcloudphxx::lgrngn::backend_t)params.backend, 
        params.cloudph_opts_init
      ));

	prtcls->init(
	  make_arrinfo(this->mem->advectee(ix::th)),
	  make_arrinfo(this->mem->advectee(ix::rv)),
	  make_arrinfo(*params.rhod)
	); 

      // writing diagnostic data for the initial condition
      diag();
    }
    // TODO: barrier?
  }

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();
    // kinematic momentum flux  = -u_fric^2 * u_i / |U| * exponential decay
    for (int ji = this->j.first(); ji <= this->j.last(); ++ji)
    {
      F(this->i, ji) = - pow(setup::u_fric,2) * this->state(ix::vip_i)(this->i, 0) / sqrt(
                          pow2(this->state(ix::vip_i)(this->i, 0)))
                          * (*params.hgt_fctr_vctr)(this->i, ji);
    }

    // du/dt = sum of kinematic momentum fluxes * dt
    int nz = this->mem->grid_size[1].length(); //76
    blitz::Range notop(0, nz-2);
    this->vip_rhs[0](this->i, notop) = (F(this->i, notop) - F(this->i, notop+1)) / this->dj * this->dt;
    this->vip_rhs[0](this->i, this->j.last()) = (F(this->i, this->j.last())) / this->dj * this->dt;
    // top and bottom cells are two times lower
    this->vip_rhs[0](this->i, 0) *= 2; 
    this->vip_rhs[0](this->i, this->j.last()) *= 2; 
  }

  void buoyancy(typename parent_t::arr_t &th, typename parent_t::arr_t &rv);
  void radiation(typename parent_t::arr_t &rv);
  void rv_src();
  void th_src(typename parent_t::arr_t &rv);
  void w_src(typename parent_t::arr_t &th, typename parent_t::arr_t &rv);
  void surf_sens();
  void surf_latent();
  void subsidence(const int&);

  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  )   
  {   
    parent_t::update_rhs(rhs, dt, at);
    using ix = typename ct_params_t::ix;

    const auto &ijk = this->ijk;

    // forcing
    switch (at) 
    {   
      // for eulerian integration or used to init trapezoidal integration
      case (0): 
      {   
        // ---- water vapor sources ----
        rv_src();
        rhs.at(ix::rv)(ijk) += alpha(ijk) + beta(ijk) * this->state(ix::rv)(ijk); 
        
        // ---- potential temp sources ----
        th_src(this->state(ix::rv));
        rhs.at(ix::th)(ijk) += alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk); 

        // vertical velocity sources
        w_src(this->state(ix::th), this->state(ix::rv));
        rhs.at(ix::w)(ijk) += alpha(ijk);

        // horizontal velocity sources 
        // large-scale vertical wind
        for(auto type : std::set<int>{ix::u})
        {
          subsidence(type);
          rhs.at(type)(ijk) += F(ijk);
        }
        break;
      }   
      case (1): 
      // trapezoidal rhs^n+1
      {   
        // ---- water vapor sources ----
        rv_src();
        rhs.at(ix::rv)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::rv)(ijk)) / (1. - 0.5 * this->dt * beta(ijk)); 
        
        // ---- potential temp sources ----
        beta(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
        th_src(beta);
        rhs.at(ix::th)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk)) / (1. - 0.5 * this->dt * beta(ijk)); 

        // vertical velocity sources
        // temporarily use beta to store the th^n+1 estimate
        beta(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
        // temporarily use F to store the rv^n+1 estimate
        F(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
        w_src(beta, F);
        rhs.at(ix::w)(ijk) += alpha(ijk);

        // horizontal velocity sources 
        // large-scale vertical wind
        for(auto type : std::set<int>{ix::u})
        {
          subsidence(type);
          rhs.at(type)(ijk) += F(ijk);
        }
        break;
      }
    }  
  }


#if defined(STD_FUTURE_WORKS)
  std::future<real_t> ftr;
#endif

  void hook_ante_step()
  {
    parent_t::hook_ante_step(); // includes output
    this->mem->barrier();
    if (this->rank == 0) 
    {
      // assuring previous async step finished ...
#if defined(STD_FUTURE_WORKS)
      if (
        params.async && 
        this->timestep != 0 && // ... but not in first timestep ...
        ((this->timestep ) % this->outfreq != 0) // ... and not after diag call, note: timestep is updated after ante_step
      ) {
        assert(ftr.valid());
        this->prec_vol += ftr.get();
      } else assert(!ftr.valid()); 
#endif

      // store liquid water content to be used in update_rhs (if done in update_rhs, it fails on async runs)
      int nx = this->mem->grid_size[0].length(); //76
      int nz = this->mem->grid_size[1].length(); //76
      prtcls->diag_all();
      prtcls->diag_wet_mom(3);
      auto rl = r_l(blitz::Range(0,nx-1), blitz::Range(0,nz-1)); 
      rl = typename parent_t::arr_t(prtcls->outbuf(), blitz::shape(nx, nz), blitz::duplicateData); // copy in data from outbuf; total liquid third moment of wet radius per kg of dry air [m^3 / kg]
      rl = rl * 4./3. * 1000. * 3.14159; // get mixing ratio [kg/kg]
      // in radiation parametrization we integrate mixing ratio * this->rhod
      rl = rl * (*params.rhod);
      rl = rl * setup::heating_kappa;

      {
        using libmpdataxx::arakawa_c::h;
        // temporarily Cx & Cz are multiplied by this->rhod ...
        auto 
          Cx = this->mem->GC[0](
            this->mem->grid_size[0]^h, 
            this->mem->grid_size[1]
          ).reindex({0,0}).copy(),
          Cz = this->mem->GC[1](
            this->mem->grid_size[0], 
            this->mem->grid_size[1]^h
          ).reindex({0,0}).copy();

        // ... and now dividing them by this->rhod (TODO: z=0 is located at k=1/2)
        {
          blitz::Range all = blitz::Range::all();
          Cx(blitz::Range(1,nx), all) /= *params.rhod;
          Cz(all, blitz::Range(1,nz)) /= *params.rhod;
          Cx(0, all) /= (*params.rhod)(0, all);
          Cz(all, 0) /= (*params.rhod)(all, 0);
        }
        // running synchronous stuff
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::th)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          make_arrinfo(*params.rhod),
          make_arrinfo(Cx), // ix::u ?
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
          make_arrinfo(Cz) // ix:w ?
        );
        // artificially remove negative rv...
        // this->mem->advectee(ix::rv) = where(this->mem->advectee(ix::rv) < 0., 0., this->mem->advectee(ix::rv));
      } 
      // running asynchronous stuff
      {
        using libcloudphxx::lgrngn::particles_t;
        using libcloudphxx::lgrngn::CUDA;
        using libcloudphxx::lgrngn::multi_CUDA;
#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(!ftr.valid());
          if(params.backend == CUDA)
            ftr = std::async(
              std::launch::async, 
              &particles_t<real_t, CUDA>::step_async, 
              dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
              params.cloudph_opts
            );
          else if(params.backend == multi_CUDA)
            ftr = std::async(
              std::launch::async, 
              &particles_t<real_t, multi_CUDA>::step_async, 
              dynamic_cast<particles_t<real_t, multi_CUDA>*>(prtcls.get()),
              params.cloudph_opts
            );
          assert(ftr.valid());
        } else 
#endif
          this->prec_vol += prtcls->step_async(params.cloudph_opts);
      }
    }
  }
  // 
  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes output
    this->mem->barrier();

    if (this->rank == 0) 
    {
      // performing diagnostics
      if (this->timestep % this->outfreq == 0) 
      { 
#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(ftr.valid());
          this->prec_vol += ftr.get();
        }
#endif
        diag();
      }
    }
    this->mem->barrier();
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    int backend = -1;
    bool async = true;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    libcloudphxx::lgrngn::opts_init_t<real_t> cloudph_opts_init;
    outmom_t<real_t> out_dry, out_wet;
    bool flag_coal; // do we want coal after spinup
  };

  private:

  // per-thread copy of params
  rt_params_t params;

  public:

  // ctor
  slvr_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p),
    tmp1(args.mem->tmp[__FILE__][0][0]),
    r_l(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    F(args.mem->tmp[__FILE__][0][1])
  {
    int nx = this->mem->grid_size[0].length();
    k_i.resize(nx);
    r_l = 0.;

    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 5); // tmp1, r_l, alpha, beta, F
  }
};

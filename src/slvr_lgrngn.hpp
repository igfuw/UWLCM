#pragma once
#include <chrono>
#include "slvr_dim.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
class slvr_lgrngn : public slvr_dim<ct_params_t>
{
  using parent_t = slvr_dim<ct_params_t>; 
  using clock = std::chrono::high_resolution_clock; // TODO: option to disable timing, as it may affect performance a little?

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  // timing fields
  clock::time_point tbeg, tend, tbeg1, tend1, tbeg_loop;
  std::chrono::milliseconds tdiag, tupdate, tsync, tasync, tasync_wait, tloop, tvip_rhs; 

  // array with index of inversion
  blitz::Array<real_t, parent_t::n_dims-1> k_i;

  // global arrays, shared among threads, TODO: in fact no need to share them?
  typename parent_t::arr_t &tmp1,
                           &tmp2,
                           &r_l,
                           &F,       // forcings helper
                           &alpha,   // 'explicit' rhs part - does not depend on the value at n+1
                           &beta;    // 'implicit' rhs part - coefficient of the value at n+1
  // helper methods
  void diag()
  {
    assert(this->rank == 0);
    tbeg = clock::now();

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
    tend = clock::now();
    tdiag += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
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
    //params.cloudph_opts.adve = val;
    params.w_src = val;
    params.cloudph_opts.coal = val ? params.flag_coal : false;
    params.cloudph_opts.RH_max = val ? 44 : 1.06; // 0.5% limit during spinup // TODO: specify it somewhere else, dup in blk_2m
  };

  // deals with nitial supersaturation
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

      // temporary array of densities - prtcls cant be init'd with 1D profile
      typename parent_t::arr_t rhod(this->mem->advectee(ix::th).shape());
      rhod = (*params.rhod)(this->vert_idx);

	prtcls->init(
	  make_arrinfo(this->mem->advectee(ix::th)),
	  make_arrinfo(this->mem->advectee(ix::rv)),
	  make_arrinfo(rhod)
	); 

      // writing diagnostic data for the initial condition
      tbeg_loop = clock::now();
      diag();
    }
    // TODO: barrier?
  }

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();
/*
    this->mem->barrier();
    if(this->rank == 0)
      tbeg = clock::now();
    // kinematic momentum flux  = -u_fric^2 * u_i / |U| * exponential decay
    typename parent_t::arr_sub_t U_ground(this->shape(this->hrzntl_subdomain));
    U_ground = this->calc_U_ground();

    // loop over horizontal dimensions
    for(int it = 0; it < parent_t::n_dims-1; ++it)
    {
      F(this->ijk).reindex(this->zero) = 
        -pow(setup::u_fric,2) *  // const, cache it
        this->vip_ground[it](blitz::tensor::i, blitz::tensor::j) /              // u_i at z=0
        U_ground(blitz::tensor::i, blitz::tensor::j) *  // |U| at z=0
        (*params.hgt_fctr_vctr)(this->vert_idx);                                       // hgt_fctr

      // du/dt = sum of kinematic momentum fluxes * dt
      this->vert_grad_fwd(F, this->vip_rhs[it], params.dz);
      this->vip_rhs[it](this->ijk) *= -params.dt;
    }
    this->mem->barrier();
    if(this->rank == 0)
    {
      tend = clock::now();
      tvip_rhs += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
    }
*/
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
    this->mem->barrier();
    if(this->rank == 0)
      tbeg = clock::now();

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
        if(params.w_src)
        {
          w_src(this->state(ix::th), this->state(ix::rv));
          rhs.at(ix::w)(ijk) += alpha(ijk);
        }

        // horizontal velocity sources 
        // large-scale vertical wind
        if(params.uv_src)
        {
          for(auto type : this->hori_vel)
          {
            subsidence(type);
            rhs.at(type)(ijk) += F(ijk);
          }
        }
        break;
      }   
      case (1): 
      // trapezoidal rhs^n+1
      {   
        // ---- water vapor sources ----
        rv_src();
        rhs.at(ix::rv)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::rv)(ijk)) / (1. - 0.5 * this->dt * beta(ijk)); 
        // TODO: alpha should also take (possibly impolicit) estimate of rv^n+1 too
        //       becomes important when nudging is introduced?
                                                                                                                                            
  
        // ---- potential temp sources ----
        tmp2(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
        // todo: once rv_src beta!=0 (e.g. nudging), rv^n+1 estimate should be implicit here
        th_src(tmp2);
        rhs.at(ix::th)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk)) / (1. - 0.5 * this->dt * beta(ijk));
        // TODO: alpha should also take (possibly impolicit) estimate of th^n+1 too
        //       becomes important when nudging is introduced?

        // vertical velocity sources
        if(params.w_src)
        {
          // temporarily use beta to store the th^n+1 estimate
          beta(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
          // todo: oncethv_src beta!=0 (e.g. nudging), th^n+1 estimate should be implicit here

          // temporarily use F to store the rv^n+1 estimate
          F(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
          // todo: once rv_src beta!=0 (e.g. nudging), rv^n+1 estimate should be implicit here

          w_src(beta, F);
          rhs.at(ix::w)(ijk) += alpha(ijk);
        }

        // horizontal velocity sources 
        // large-scale vertical wind
        if(params.uv_src)
        {
          for(auto type : this->hori_vel)
          {
            subsidence(type); // TODO: in case 1 type here should be in step n+1, calc it explicitly as type + 0.5 * dt * rhs(type);
                              //       could also be calculated implicitly, but we would need implicit type^n+1 in other cells
                              //       also include absorber in type^n+1 estimate...
            rhs.at(type)(ijk) += F(ijk);
          }
        }
        break;
      }
    }  
    this->mem->barrier();
    if(this->rank == 0)
    {
      tend = clock::now();
      tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
    }
  }


#if defined(STD_FUTURE_WORKS)
  std::future<real_t> ftr;
#endif

  void hook_ante_step()
  {
//std::cout << this->timestep << std::endl;
    parent_t::hook_ante_step(); // includes output
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tbeg1 = clock::now();
      // assuring previous async step finished ...
#if defined(STD_FUTURE_WORKS)
      if (
        params.async && 
        this->timestep != 0 && // ... but not in first timestep ...
        ((this->timestep ) % this->outfreq != 0) // ... and not after diag call, note: timestep is updated after ante_step
      ) {
        assert(ftr.valid());
        tbeg = clock::now();
        this->prec_vol += ftr.get();
        tend = clock::now();
        tasync_wait += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      } else assert(!ftr.valid()); 
#endif

      // store liquid water content to be used in update_rhs (if done in update_rhs, it fails on async runs)
      prtcls->diag_all();
      prtcls->diag_wet_mom(3);
      auto rl = r_l(this->domain); // rl refrences subdomain of r_l
      rl = typename parent_t::arr_t(prtcls->outbuf(), rl.shape(), blitz::duplicateData); // copy in data from outbuf; total liquid third moment of wet radius per kg of dry air [m^3 / kg]
      rl = rl * 4./3. * 1000. * 3.14159; // get mixing ratio [kg/kg]
      {
        // temporarily Cx & Cz are multiplied by this->rhod ...
        auto 
          Cx = this->mem->GC[0](this->Cx_domain).reindex(this->zero).copy(),
          Cy = this->mem->GC[1](this->Cy_domain).reindex(this->zero).copy(),
          Cz = this->mem->GC[this->vert_dim](this->Cz_domain).reindex(this->zero).copy();

        // ... and now dividing them by this->rhod (TODO: z=0 is located at k=1/2)
        {
          Cx.reindex(this->zero) /= (*params.rhod)(this->vert_idx);
          Cy.reindex(this->zero) /= (*params.rhod)(this->vert_idx);
          Cz.reindex(this->zero) /= (*params.rhod)(this->vert_idx); // TODO: should be interpolated, since theres a shift between positions of rhod and Cz
        }

        auto Cy_arrinfo = this->n_dims == 2 ? 
          libcloudphxx::lgrngn::arrinfo_t<real_t>() : // empty for 2D run
          make_arrinfo(Cy);

        // running synchronous stuff
        tbeg = clock::now();
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::th)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
          make_arrinfo(Cx),
          Cy_arrinfo,
          make_arrinfo(Cz)
        );
if(!std::isfinite(sum(this->mem->advectee(ix::th))))
  std::cout << "nan in th: " << this->mem->advectee(ix::th);
if(!std::isfinite(sum(this->mem->advectee(ix::rv))))
  std::cout << "nan in rv: " << this->mem->advectee(ix::rv);
        tend = clock::now();
        tsync += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      } 
      // running asynchronous stuff
      {
        using libcloudphxx::lgrngn::particles_t;
        using libcloudphxx::lgrngn::CUDA;
        using libcloudphxx::lgrngn::multi_CUDA;
        tbeg = clock::now();
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
        tend = clock::now();
        tasync += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
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
      tend1 = clock::now();
      // there's no hook_post_loop, so we imitate it here to write out computation times
      if(this->timestep == params.nt-1)
      {
        tend = clock::now();
        tloop = std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg_loop );
        std::cout <<  "wall time in milliseconds: " << std::endl
          << "loop: " << tloop.count() << std::endl
          << "custom rhs update: " << tupdate.count() << " ("<< setup::real_t(tupdate.count())/tloop.count()*100 <<"%)" << std::endl
          << "custom vip_rhs: " << tvip_rhs.count() << " ("<< setup::real_t(tvip_rhs.count())/tloop.count()*100 <<"%)" << std::endl
          << "diag: " << tdiag.count() << " ("<< setup::real_t(tdiag.count())/tloop.count()*100 <<"%)" << std::endl
          << "sync: " << tsync.count() << " ("<< setup::real_t(tsync.count())/tloop.count()*100 <<"%)" << std::endl
          << "async: " << tasync.count() << " ("<< setup::real_t(tasync.count())/tloop.count()*100 <<"%)" << std::endl
          << "async_wait: " << tasync_wait.count() << " ("<< setup::real_t(tasync_wait.count())/tloop.count()*100 <<"%)" << std::endl;
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
    tmp2(args.mem->tmp[__FILE__][0][5]),
    r_l(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    F(args.mem->tmp[__FILE__][0][1])
  {
    k_i.resize(this->shape(this->hrzntl_domain));
    r_l = 0.;

    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 6); // tmp1, tmp2, r_l, alpha, beta, F
  }
};

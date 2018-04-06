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

  // helper methods
  void diag()
  {
    assert(this->rank == 0);
    parent_t::tbeg = parent_t::clock::now();

    // recording super-droplet concentration per grid cell 
    prtcls->diag_all();
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());

    // recording concentration of SDs that represent activated droplets 
    prtcls->diag_rw_ge_rc();
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc_act", prtcls->outbuf());

    // recording relative humidity
    prtcls->diag_RH();
    this->record_aux("RH", prtcls->outbuf());

    // recording precipitation rate per grid cel
    prtcls->diag_all();
    prtcls->diag_precip_rate();
    this->record_aux("precip_rate", prtcls->outbuf());

    // recording 1st mom of rw of gccns
    /*
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
    */
    // recording 0th mom of rw of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_wet_mom(0);
    this->record_aux("actrw_rw_mom0", prtcls->outbuf());

    // recording 1st mom of rw of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_wet_mom(1);
    this->record_aux("actrw_rw_mom1", prtcls->outbuf());

    // recording 2nd mom of rw of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_wet_mom(2);
    this->record_aux("actrw_rw_mom2", prtcls->outbuf());

    // recording 3rd mom of rw of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_wet_mom(3);
    this->record_aux("actrw_rw_mom3", prtcls->outbuf());

    // recording 1st mom of rd of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_dry_mom(1);
    this->record_aux("actrw_rd_mom1", prtcls->outbuf());

    // recording 0th mom of rd of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_dry_mom(0);
    this->record_aux("actrw_rd_mom0", prtcls->outbuf());
    /*
    // recording 1st mom of rd of activated drops
    prtcls->diag_RH_ge_Sc();
    prtcls->diag_dry_mom(1);
    this->record_aux("actRH_rd_mom1", prtcls->outbuf());
   
    // recording 3rd mom of rw of activated drops
    prtcls->diag_RH_ge_Sc();
    prtcls->diag_wet_mom(3);
    this->record_aux("actRH_rw_mom3", prtcls->outbuf());

    // recording 0th mom of rd of activated drops
    prtcls->diag_RH_ge_Sc();
    prtcls->diag_dry_mom(0);
    this->record_aux("actRH_rd_mom0", prtcls->outbuf());
    */
    // recording 0th wet mom of radius of rain drops (r>25um)
    prtcls->diag_wet_rng(25.e-6, 1);
    prtcls->diag_wet_mom(0);
    this->record_aux("rain_rw_mom0", prtcls->outbuf());
    
    // recording 3rd wet mom of radius of rain drops (r>25um)
    prtcls->diag_wet_rng(25.e-6, 1);
    prtcls->diag_wet_mom(3);
    this->record_aux("rain_rw_mom3", prtcls->outbuf());

    // recording 0th wet mom of radius of cloud drops (.5um< r < 25um)
    prtcls->diag_wet_rng(.5e-6, 25.e-6);
    prtcls->diag_wet_mom(0);
    this->record_aux("cloud_rw_mom0", prtcls->outbuf());

    // recording 3rd wet mom of radius of cloud drops (.5um< r < 25um)
    prtcls->diag_wet_rng(.5e-6, 25.e-6);
    prtcls->diag_wet_mom(3);
    this->record_aux("cloud_rw_mom3", prtcls->outbuf());

    // recording 3rd wet mom of radius of aerosols (r < .5um)
    prtcls->diag_wet_rng(0., .5e-6);
    prtcls->diag_wet_mom(3);
    this->record_aux("aerosol_rw_mom3", prtcls->outbuf());
   
    // recording divergence of the velocity field
    prtcls->diag_vel_div();
    this->record_aux("vel_div", prtcls->outbuf());

    // recording puddle
    auto puddle = prtcls->diag_puddle();
    for(auto elem : puddle)
    {   
       this->f_puddle << elem.first << " " << elem.second << "\n";
    }   
    this->f_puddle << "\n";
   
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
    parent_t::tend = parent_t::clock::now();
    parent_t::tdiag += std::chrono::duration_cast<std::chrono::milliseconds>( parent_t::tend - parent_t::tbeg );
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
    params.cloudph_opts.RH_max = val ? 44 : 1.01; // TODO: specify it somewhere else, dup in blk_2m
  };

  // deals with nitial supersaturation
  void hook_ante_loop(int nt)
  {
    params.flag_coal = params.cloudph_opts.coal;

    // TODO: barrier?
    this->mem->barrier();
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

        if(params.cloudph_opts_init.sd_conc)
          params.cloudph_opts_init.n_sd_max = params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * params.cloudph_opts_init.sd_conc;
        else
          params.cloudph_opts_init.n_sd_max = 1.1 * params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * 1.e8 * params.cloudph_opts_init.dx * params.cloudph_opts_init.dz / params.cloudph_opts_init.sd_const_multi; // hardcoded N_a=100/cm^3 !!
          
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

        if(params.cloudph_opts_init.sd_conc)
          params.cloudph_opts_init.n_sd_max = params.cloudph_opts_init.nx * params.cloudph_opts_init.ny *params.cloudph_opts_init.nz * params.cloudph_opts_init.sd_conc;
        else
          params.cloudph_opts_init.n_sd_max = 1.1 * params.cloudph_opts_init.nx * params.cloudph_opts_init.ny * params.cloudph_opts_init.nz * 1.e8 * params.cloudph_opts_init.dx * params.cloudph_opts_init.dy * params.cloudph_opts_init.dz / params.cloudph_opts_init.sd_const_multi; // hardcoded N_a=100/cm^3 !!

        if(params.backend == libcloudphxx::lgrngn::multi_CUDA)
          params.cloudph_opts_init.n_sd_max *= 1.5; // more space for copied SDs
      }

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
        (libcloudphxx::lgrngn::backend_t)params.backend, 
        params.cloudph_opts_init
      ));

      // temporary array of densities - prtcls cant be init'd with 1D profile
      typename parent_t::arr_t rhod(this->mem->advectee(ix::th).shape()); // TODO: replace all rhod arrays with this->mem->G
      rhod = (*params.rhod)(this->vert_idx);

	prtcls->init(
	  make_arrinfo(this->mem->advectee(ix::th)),
	  make_arrinfo(this->mem->advectee(ix::rv)),
	  make_arrinfo(rhod)
	); 

      // writing diagnostic data for the initial condition
      parent_t::tbeg_loop = parent_t::clock::now();
    }
    this->mem->barrier();
    parent_t::hook_ante_loop(nt); 
    // record microphysics options
    // TODO: divide them into groups
    // TODO: add recording of dry_distros, src_dry_distros, dry_sizes, kernel_parameters
    // TODO: record meaningful names of kernel, adve_scheme, terminal_velocity, backend
    // TODO: add git revisions of libmpdata, libcloud and uwlcm
    if (this->rank == 0) 
    {
      this->record_aux_const("super-droplet microphysics", -44);  
      this->record_aux_const("nx", params.cloudph_opts_init.nx);  
      this->record_aux_const("ny", params.cloudph_opts_init.ny);  
      this->record_aux_const("nz", params.cloudph_opts_init.nz);  
      this->record_aux_const("dx", params.cloudph_opts_init.dx);  
      this->record_aux_const("dy", params.cloudph_opts_init.dy);  
      this->record_aux_const("dz", params.cloudph_opts_init.dz);  
      this->record_aux_const("dt", params.cloudph_opts_init.dt);  
      this->record_aux_const("x0", params.cloudph_opts_init.x0);  
      this->record_aux_const("y0", params.cloudph_opts_init.y0);  
      this->record_aux_const("z0", params.cloudph_opts_init.z0);  
      this->record_aux_const("x1", params.cloudph_opts_init.x1);  
      this->record_aux_const("y1", params.cloudph_opts_init.y1);  
      this->record_aux_const("z1", params.cloudph_opts_init.z1);  
      this->record_aux_const("aerosol_independent_of_rhod", params.cloudph_opts_init.aerosol_independent_of_rhod);  
      this->record_aux_const("sd_conc", params.cloudph_opts_init.sd_conc);  
      this->record_aux_const("sd_conc_large_tail", params.cloudph_opts_init.sd_conc_large_tail);  
      this->record_aux_const("sd_const_multi", params.cloudph_opts_init.sd_const_multi);  
      this->record_aux_const("sd_const_multi_dry_sizes", params.cloudph_opts_init.sd_const_multi_dry_sizes);  
      this->record_aux_const("n_sd_max", params.cloudph_opts_init.n_sd_max);  
      this->record_aux_const("adve_scheme", params.cloudph_opts_init.adve_scheme);  
      this->record_aux_const("dev_count", params.cloudph_opts_init.dev_count);  
      this->record_aux_const("dev_id", params.cloudph_opts_init.dev_id);  
      this->record_aux_const("sstp_cond", params.cloudph_opts_init.sstp_cond);  
      this->record_aux_const("sstp_coal", params.cloudph_opts_init.sstp_coal);  
      this->record_aux_const("sstp_chem", params.cloudph_opts_init.sstp_chem);  
      this->record_aux_const("exact_sstp_cond", params.cloudph_opts_init.exact_sstp_cond);  
      this->record_aux_const("rng_seed", params.cloudph_opts_init.rng_seed);  
      this->record_aux_const("kernel", params.cloudph_opts_init.kernel);  
      this->record_aux_const("terminal_velocity", params.cloudph_opts_init.terminal_velocity);  
      this->record_aux_const("backend", params.backend);  
      this->record_aux_const("async", params.async);  
      this->record_aux_const("adve", params.cloudph_opts.adve);  
      this->record_aux_const("sedi", params.cloudph_opts.sedi);  
      this->record_aux_const("cond", params.cloudph_opts.cond);  
      this->record_aux_const("coal", params.flag_coal);  // cloudph_opts.coal could be 0 here due to spinup
      this->record_aux_const("rcyc", params.cloudph_opts.rcyc);  
      this->record_aux_const("chem_switch", params.cloudph_opts_init.chem_switch);  
      this->record_aux_const("coal_switch", params.cloudph_opts_init.coal_switch);  
      this->record_aux_const("sedi_switch", params.cloudph_opts_init.sedi_switch);  
      this->record_aux_const("src_switch", params.cloudph_opts_init.src_switch);  
      this->record_aux_const("chem_dsl", params.cloudph_opts.chem_dsl);  
      this->record_aux_const("chem_dsc", params.cloudph_opts.chem_dsc);  
      this->record_aux_const("chem_rct", params.cloudph_opts.chem_rct);  
      this->record_aux_const("chem_rho", params.cloudph_opts_init.chem_rho);  
      this->record_aux_const("opts_init RH_max", params.cloudph_opts_init.RH_max);  
      this->record_aux_const("supstp_src", params.cloudph_opts_init.supstp_src);  
      this->record_aux_const("src_sd_conc", params.cloudph_opts_init.src_sd_conc);  
      this->record_aux_const("src_z1", params.cloudph_opts_init.src_z1);  
    }

    // TODO: barrier?
  }

#if defined(STD_FUTURE_WORKS)
  std::future<void> ftr;
#endif

  void hook_ante_step()
  {
    parent_t::hook_ante_step(); // includes output
    this->mem->barrier();
    if (this->rank == 0) 
    {
      parent_t::tbeg1 = parent_t::clock::now();
      // assuring previous async step finished ...
#if defined(STD_FUTURE_WORKS)
      if (
        params.async && 
        this->timestep != 0 && // ... but not in first timestep ...
        ((this->timestep ) % this->outfreq != 0) // ... and not after diag call, note: timestep is updated after ante_step
      ) {
        assert(ftr.valid());
        parent_t::tbeg = parent_t::clock::now();
        ftr.get();
        parent_t::tend = parent_t::clock::now();
        parent_t::tasync_wait += std::chrono::duration_cast<std::chrono::milliseconds>( parent_t::tend - parent_t::tbeg );
      } else assert(!ftr.valid()); 
#endif

      // store liquid water content to be used in update_rhs (if done in update_rhs, it fails on async runs)
      prtcls->diag_all();
      prtcls->diag_wet_mom(3);
      auto rl = parent_t::r_l(this->domain); // rl refrences subdomain of r_l
      rl = typename parent_t::arr_t(prtcls->outbuf(), rl.shape(), blitz::duplicateData); // copy in data from outbuf; total liquid third moment of wet radius per kg of dry air [m^3 / kg]
      nancheck(rl, "rl after copying from diag_wet_mom(3)");
      rl = rl * 4./3. * 1000. * 3.14159; // get mixing ratio [kg/kg]

      {
        // temporarily Cx & Cz are multiplied by this->rhod ...
        auto 
          Cx = this->mem->GC[0](this->Cx_domain).reindex(this->zero).copy(),
          Cy = this->mem->GC[1](this->Cy_domain).reindex(this->zero).copy(),
          Cz = this->mem->GC[ix::w](this->Cz_domain).reindex(this->zero).copy(); 
        nancheck(Cx, "Cx after copying from mpdata");
        nancheck(Cy, "Cy after copying from mpdata");
        nancheck(Cz, "Cz after copying from mpdata");

        // ... and now dividing them by this->rhod (TODO: z=0 is located at k=1/2)
        {
          Cx.reindex(this->zero) /= (*params.rhod)(this->vert_idx);
          Cy.reindex(this->zero) /= (*params.rhod)(this->vert_idx);
          Cz.reindex(this->zero) /= (*params.rhod)(this->vert_idx); // TODO: should be interpolated, since theres a shift between positions of rhod and Cz
        }

        // running synchronous stuff
        parent_t::tbeg = parent_t::clock::now();

        // rv might be negative due to large negative RHS from SD fluctuations + large-scale subsidence?
        // turn all negative rv into rv = 0... CHEATING
        negtozero(this->mem->advectee(ix::rv), "rv before step sync");

        nancheck(this->mem->advectee(ix::th), "th before step sync");
        nancheck(this->mem->advectee(ix::rv), "rv before step sync");
        negcheck(this->mem->advectee(ix::th), "th before step sync");
        negcheck(this->mem->advectee(ix::rv), "rv before step sync");
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::th)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
          make_arrinfo(Cx),
          this->n_dims == 2 ? libcloudphxx::lgrngn::arrinfo_t<real_t>() : make_arrinfo(Cy),
          make_arrinfo(Cz)
        );
       
        // microphysics could have led to rv < 0 ?
        negtozero(this->mem->advectee(ix::rv), "rv after step sync");

        nancheck(this->mem->advectee(ix::th), "th after step sync");
        nancheck(this->mem->advectee(ix::rv), "rv after step sync");
        negcheck(this->mem->advectee(ix::th), "th after step sync");
        negcheck(this->mem->advectee(ix::rv), "rv after step sync");

        parent_t::tend = parent_t::clock::now();
        parent_t::tsync += std::chrono::duration_cast<std::chrono::milliseconds>( parent_t::tend - parent_t::tbeg );
      } 
      // running asynchronous stuff
      {
        using libcloudphxx::lgrngn::particles_t;
        using libcloudphxx::lgrngn::CUDA;
        using libcloudphxx::lgrngn::multi_CUDA;
        parent_t::tbeg = parent_t::clock::now();
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
          prtcls->step_async(params.cloudph_opts);
        parent_t::tend = parent_t::clock::now();
        parent_t::tasync += std::chrono::duration_cast<std::chrono::milliseconds>( parent_t::tend - parent_t::tbeg );
      }
    }
    this->mem->barrier();
  }
  
  void record_all()
  {
    // plain (no xdmf) hdf5 output
    parent_t::output_t::parent_t::record_all();
    // UWLCM output
#if defined(STD_FUTURE_WORKS)
    if (this->timestep > 0 && params.async)
    {
      assert(ftr.valid());
      ftr.get();
    }
#endif
    diag();
    // xmf markup
    this->write_xmfs();
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
    params(p)
  {

    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

};

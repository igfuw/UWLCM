#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_ante_loop(int nt)
{
  params.flag_coal = params.cloudph_opts.coal;

  this->mem->barrier();

  this->generate_stretching_parameters(std::random_device{}(), libmpdataxx::formulae::fractal::stretch_params::d_distro_t::LES_th_supersaturated);
  this->reconstruct_refinee(ix::th);
  this->generate_stretching_parameters(std::random_device{}(), libmpdataxx::formulae::fractal::stretch_params::d_distro_t::LES_rv_supersaturated);
  this->reconstruct_refinee(ix::rv);
//  this->mem->refinee(this->ix_r2r.at(ix::th)) = 300;
//  this->mem->refinee(this->ix_r2r.at(ix::rv)) = 1e-3;
  negtozero(this->mem->refinee(this->ix_r2r.at(ix::rv))(this->ijk_ref), "refined rv in hook_ante_loop");

  if (this->rank == 0) 
  {
    assert(params.backend != libcloudphxx::lgrngn::undefined);
    assert(params.dt != 0); 

    if(params.gccn > 0)
      params.cloudph_opts.src = true;

    params.cloudph_opts.rlx = false;

    // async does not make sense without CUDA
    if (params.backend != libcloudphxx::lgrngn::CUDA && params.backend != libcloudphxx::lgrngn::multi_CUDA) params.async = false;

    params.cloudph_opts_init.dt = params.dt; // advection timestep = microphysics timestep

    // NOTE: this assumes that refinement is done! TODO: what if there is no refinement? i.e. n_ref_iter==0?
    // TODO: DRY
    params.cloudph_opts_init.nx = this->mem->grid_size_ref[0].length();
    params.cloudph_opts_init.dx = this->di / this->mem->n_ref;

    params.cloudph_opts_init.x0 = params.cloudph_opts_init.dx / 2;
    params.cloudph_opts_init.x1 = (params.cloudph_opts_init.nx - .5) *  params.cloudph_opts_init.dx;

    int n_sd_from_dry_sizes = 0;
    for (auto const& krcm : params.cloudph_opts_init.dry_sizes)
      for (auto const& rcm : krcm.second)
        n_sd_from_dry_sizes += rcm.second.second;

    // src_dry_distros is used only in the first step to add GCCN below some level
    int n_sd_from_src_dry_distros = params.cloudph_opts_init.src_sd_conc * params.cloudph_opts_init.src_z1 / params.cloudph_opts_init.z1 + 0.5;
      
    const int n_sd_per_cell = params.cloudph_opts_init.sd_conc + n_sd_from_dry_sizes + n_sd_from_src_dry_distros;

    if(parent_t::n_dims == 2) // 2D
    {
      params.cloudph_opts_init.nz = this->mem->grid_size_ref[1].length();
      params.cloudph_opts_init.dz = this->dj / this->mem->n_ref;
      params.cloudph_opts_init.z0 = params.cloudph_opts_init.dz / 2;
      params.cloudph_opts_init.z1 = (params.cloudph_opts_init.nz - .5) * params.cloudph_opts_init.dz;

      if(params.cloudph_opts_init.sd_conc)
      {
        if(params.cloudph_opts_init.sd_conc_large_tail)
          params.cloudph_opts_init.n_sd_max = 1.2 * params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * n_sd_per_cell; /// 1.2 to make space for large tail
        else
          params.cloudph_opts_init.n_sd_max = params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * n_sd_per_cell;
      }
      else
        params.cloudph_opts_init.n_sd_max = 1.2 * params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * 1.e8 * params.cloudph_opts_init.dx * params.cloudph_opts_init.dz / params.cloudph_opts_init.sd_const_multi; // NOTE: hardcoded N_a=100/cm^3 !!
        
      if(params.backend == libcloudphxx::lgrngn::multi_CUDA || this->mem->distmem.size()>1)
        params.cloudph_opts_init.n_sd_max *= 1.4; // more space for copied SDs
    }
    else // 3D
    {
      params.cloudph_opts_init.ny = this->mem->grid_size_ref[1].length();
      params.cloudph_opts_init.dy = this->dj / this->mem->n_ref;
      params.cloudph_opts_init.y0 = params.cloudph_opts_init.dy / 2;
      params.cloudph_opts_init.y1 = (params.cloudph_opts_init.ny - .5) * params.cloudph_opts_init.dy;

      params.cloudph_opts_init.nz = this->mem->grid_size_ref[2].length();
      params.cloudph_opts_init.dz = this->dk / this->mem->n_ref;
      params.cloudph_opts_init.z0 = params.cloudph_opts_init.dz / 2;
      params.cloudph_opts_init.z1 = (params.cloudph_opts_init.nz - .5) * params.cloudph_opts_init.dz;

      if(params.cloudph_opts_init.sd_conc)
      {
        if(params.cloudph_opts_init.sd_conc_large_tail)
          params.cloudph_opts_init.n_sd_max = 1.2 * params.cloudph_opts_init.nx * params.cloudph_opts_init.ny * params.cloudph_opts_init.nz * n_sd_per_cell; /// 1.2 to make space for large tail
        else
          params.cloudph_opts_init.n_sd_max =       params.cloudph_opts_init.nx * params.cloudph_opts_init.ny * params.cloudph_opts_init.nz * n_sd_per_cell; 
      }
      else
        params.cloudph_opts_init.n_sd_max = 1.2 * params.cloudph_opts_init.nx * params.cloudph_opts_init.ny * params.cloudph_opts_init.nz * 1.e8 * params.cloudph_opts_init.dx * params.cloudph_opts_init.dy * params.cloudph_opts_init.dz / params.cloudph_opts_init.sd_const_multi; // hardcoded N_a=100/cm^3 !!

      if(params.backend == libcloudphxx::lgrngn::multi_CUDA || this->mem->distmem.size()>1)
        params.cloudph_opts_init.n_sd_max *= 1.3; // more space for copied SDs
    }

    params.cloudph_opts_init.rlx_sd_per_bin /= this->mem->distmem.size();

    // space for SD created via relaxation, impossible to know exactly how many will be added, because it depends on washout of SD...
    int n_sd_from_rlx_dry_distros = params.cloudph_opts_init.rlx_sd_per_bin * params.cloudph_opts_init.rlx_bins * params.cloudph_opts_init.nz * 100; // room for 100 rounds of full relaxation... 
    params.cloudph_opts_init.n_sd_max += n_sd_from_rlx_dry_distros;

    prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
      (libcloudphxx::lgrngn::backend_t)params.backend, 
      params.cloudph_opts_init
    ));

    // temporary array of densities - prtcls cant be init'd with 1D profile
    typename parent_t::arr_t rhod_ref(this->mem->refinee(this->ix_r2r.at(ix::th)).shape()); // TODO: replace all rhod arrays with this->mem->G
    rhod_ref = (params.profs_ref.rhod)(this->vert_idx);

    // temporary array of pressure - prtcls cant be init'd with 1D profile
    typename parent_t::arr_t p_e_ref(this->mem->refinee(this->ix_r2r.at(ix::th)).shape());  
    p_e_ref = (params.profs_ref.p_e)(this->vert_idx);

    prtcls->init(
      make_arrinfo(this->mem->refinee(this->ix_r2r.at(ix::th))),
      make_arrinfo(this->mem->refinee(this->ix_r2r.at(ix::rv))),
      make_arrinfo(rhod_ref)
      ,make_arrinfo(p_e_ref)
    ); 
  }
  this->mem->barrier();
  parent_t::hook_ante_loop(nt); 
  // record microphysics options
  // TODO: divide them into groups
  // TODO: add recording of dry_distros, src_dry_distros, dry_sizes, kernel_parameters
  if (this->rank == 0) 
  {
    this->record_aux_const("microphysics", "super-droplets");  
    this->record_aux_const("nx", "lgrngn", params.cloudph_opts_init.nx);  
    this->record_aux_const("ny", "lgrngn", params.cloudph_opts_init.ny);  
    this->record_aux_const("nz", "lgrngn", params.cloudph_opts_init.nz);  
    this->record_aux_const("dx", "lgrngn", params.cloudph_opts_init.dx);  
    this->record_aux_const("dy", "lgrngn", params.cloudph_opts_init.dy);  
    this->record_aux_const("dz", "lgrngn", params.cloudph_opts_init.dz);  
    this->record_aux_const("dt", "lgrngn", params.cloudph_opts_init.dt);  
    this->record_aux_const("x0", "lgrngn", params.cloudph_opts_init.x0);  
    this->record_aux_const("y0", "lgrngn", params.cloudph_opts_init.y0);  
    this->record_aux_const("z0", "lgrngn", params.cloudph_opts_init.z0);  
    this->record_aux_const("x1", "lgrngn", params.cloudph_opts_init.x1);  
    this->record_aux_const("y1", "lgrngn", params.cloudph_opts_init.y1);  
    this->record_aux_const("z1", "lgrngn", params.cloudph_opts_init.z1);  
    this->record_aux_const("aerosol_independent_of_rhod", "lgrngn", params.cloudph_opts_init.aerosol_independent_of_rhod);  
    this->record_aux_const("sd_conc", "lgrngn", params.cloudph_opts_init.sd_conc);  
    this->record_aux_const("sd_conc_large_tail", "lgrngn", params.cloudph_opts_init.sd_conc_large_tail);  
    this->record_aux_const("sd_const_multi", "lgrngn", params.cloudph_opts_init.sd_const_multi);  
    this->record_aux_const("n_sd_max", "lgrngn", params.cloudph_opts_init.n_sd_max);  
    this->record_aux_const("dev_count", "lgrngn", params.cloudph_opts_init.dev_count);  
    this->record_aux_const("dev_id", "lgrngn", params.cloudph_opts_init.dev_id);  
    this->record_aux_const("sstp_cond", "lgrngn", params.cloudph_opts_init.sstp_cond);  
    this->record_aux_const("sstp_coal", "lgrngn", params.cloudph_opts_init.sstp_coal);  
    this->record_aux_const("sstp_chem", "lgrngn", params.cloudph_opts_init.sstp_chem);  
    this->record_aux_const("exact_sstp_cond", "lgrngn", params.cloudph_opts_init.exact_sstp_cond);  
    this->record_aux_const("diag_incloud_time", "lgrngn", params.cloudph_opts_init.diag_incloud_time);  
    this->record_aux_const("rng_seed", "lgrngn", params.cloudph_opts_init.rng_seed);  
    this->record_aux_const("rng_seed_init", "lgrngn", params.cloudph_opts_init.rng_seed_init);  
    this->record_aux_const("rng_seed_init_switch", "lgrngn", params.cloudph_opts_init.rng_seed_init_switch);  
    this->record_aux_const("async", "lgrngn", params.async);  
    this->record_aux_const("adve", "lgrngn", params.cloudph_opts.adve);  
    this->record_aux_const("sedi", "lgrngn", params.cloudph_opts.sedi);  
    this->record_aux_const("subs", "lgrngn", params.cloudph_opts.subs);  
    this->record_aux_const("cond", "lgrngn", params.cloudph_opts.cond);  
    this->record_aux_const("coal", "lgrngn", params.flag_coal);  // cloudph_opts.coal could be 0 here due to spinup
    this->record_aux_const("rcyc", "lgrngn", params.cloudph_opts.rcyc);  
    this->record_aux_const("src", "lgrngn", params.cloudph_opts.src);  
    this->record_aux_const("rlx", "lgrngn", params.cloudph_opts.rlx);  
    this->record_aux_const("gccn", "lgrngn", params.gccn);  
    this->record_aux_const("turb_adve", "lgrngn", params.cloudph_opts.turb_adve);  
    this->record_aux_const("turb_cond", "lgrngn", params.cloudph_opts.turb_cond);  
    this->record_aux_const("turb_coal", "lgrngn", params.cloudph_opts.turb_coal);  
    this->record_aux_const("chem_switch", "lgrngn", params.cloudph_opts_init.chem_switch);  
    this->record_aux_const("coal_switch", "lgrngn", params.cloudph_opts_init.coal_switch);  
    this->record_aux_const("sedi_switch", "lgrngn", params.cloudph_opts_init.sedi_switch);  
    this->record_aux_const("subs_switch", "lgrngn", params.cloudph_opts_init.subs_switch);  
    this->record_aux_const("rlx_switch", "lgrngn", params.cloudph_opts_init.rlx_switch);  
    this->record_aux_const("turb_adve_switch", "lgrngn", params.cloudph_opts_init.turb_adve_switch);  
    this->record_aux_const("turb_cond_switch", "lgrngn", params.cloudph_opts_init.turb_cond_switch);  
    this->record_aux_const("turb_coal_switch", "lgrngn", params.cloudph_opts_init.turb_coal_switch);  
    this->record_aux_const("chem_dsl", "lgrngn", params.cloudph_opts.chem_dsl);  
    this->record_aux_const("chem_dsc", "lgrngn", params.cloudph_opts.chem_dsc);  
    this->record_aux_const("chem_rct", "lgrngn", params.cloudph_opts.chem_rct);  
    this->record_aux_const("chem_rho", "lgrngn", params.cloudph_opts_init.chem_rho);  
    this->record_aux_const("opts_init RH_max", "lgrngn", params.cloudph_opts_init.RH_max);  
    this->record_aux_const("supstp_src", "lgrngn", params.cloudph_opts_init.supstp_src);  
    this->record_aux_const("supstp_rlx", "lgrngn", params.cloudph_opts_init.supstp_rlx);  
    this->record_aux_const("src_sd_conc", "lgrngn", params.cloudph_opts_init.src_sd_conc);  
    this->record_aux_const("src_z1", "lgrngn", params.cloudph_opts_init.src_z1);
    this->record_aux_const("rlx_bins", "lgrngn", params.cloudph_opts_init.rlx_bins);  
    this->record_aux_const("rlx_sd_per_bin", "lgrngn", params.cloudph_opts_init.rlx_sd_per_bin);  
    this->record_aux_const("rlx_timescale", "lgrngn", params.cloudph_opts_init.rlx_timescale);  
    this->record_aux_const("adve_scheme", "lgrngn", libcloudphxx::lgrngn::as_name.at(params.cloudph_opts_init.adve_scheme));  
    this->record_aux_const("backend", "lgrngn", libcloudphxx::lgrngn::backend_name.at(params.backend));  
    this->record_aux_const("kernel", "lgrngn", libcloudphxx::lgrngn::kernel_name.at(params.cloudph_opts_init.kernel));  
    this->record_aux_const("src_type", "lgrngn", libcloudphxx::lgrngn::src_name.at(params.cloudph_opts_init.src_type));  
    this->record_aux_const("terminal_velocity", "lgrngn", libcloudphxx::lgrngn::vt_name.at(params.cloudph_opts_init.terminal_velocity));  
    this->record_aux_const("RH_formula", "lgrngn", libcloudphxx::lgrngn::RH_formula_name.at(params.cloudph_opts_init.RH_formula));  
    this->record_aux_const("rd_min", "lgrngn", params.cloudph_opts_init.rd_min);  
    this->record_aux_const("rd_max", "lgrngn", params.cloudph_opts_init.rd_max);  
    this->record_aux_const("relax_ccn", "user_params", params.user_params.relax_ccn);  
  }
  this->mem->barrier();
}

#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_param_t>::hook_ante_loop(int nt)
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

    if(this->mem->distmem.rank() == 0)
      params.cloudph_opts_init.x0 = this->di / 2;
    else
      params.cloudph_opts_init.x0 = 0.;

    if(this->mem->distmem.rank() == this->mem->distmem.size()-1)
      params.cloudph_opts_init.x1 = (params.cloudph_opts_init.nx - .5) * this->di;
    else
      params.cloudph_opts_init.x1 =  params.cloudph_opts_init.nx       * this->di;

    int n_sd_from_dry_sizes = 0;
    for (auto const& krcm : params.cloudph_opts_init.dry_sizes)
      for (auto const& rcm : krcm.second)
        n_sd_from_dry_sizes += rcm.second.second;
      
    const int n_sd_per_cell = params.cloudph_opts_init.sd_conc + n_sd_from_dry_sizes;

    if(parent_t::n_dims == 2) // 2D
    {
      params.cloudph_opts_init.nz = this->mem->grid_size[1].length();
      params.cloudph_opts_init.dz = this->dj;
      params.cloudph_opts_init.z0 = this->dj / 2;
      params.cloudph_opts_init.z1 = (params.cloudph_opts_init.nz - .5) * this->dj;

      if(params.cloudph_opts_init.sd_conc)
      {
        if(params.cloudph_opts_init.sd_conc_large_tail)
          params.cloudph_opts_init.n_sd_max = 1.2 * params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * n_sd_per_cell; /// 1.2 to make space for large tail
        else
          params.cloudph_opts_init.n_sd_max = params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * n_sd_per_cell;
      }
      else
        params.cloudph_opts_init.n_sd_max = 1.2 * params.cloudph_opts_init.nx * params.cloudph_opts_init.nz * 1.e8 * params.cloudph_opts_init.dx * params.cloudph_opts_init.dz / params.cloudph_opts_init.sd_const_multi; // hardcoded N_a=100/cm^3 !!
        
      if(params.backend == libcloudphxx::lgrngn::multi_CUDA || this->mem->distmem.size()>1)
        params.cloudph_opts_init.n_sd_max *= 1.4; // more space for copied SDs
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

    prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
      (libcloudphxx::lgrngn::backend_t)params.backend, 
      params.cloudph_opts_init
    ));

    // temporary array of densities - prtcls cant be init'd with 1D profile
    typename parent_t::arr_t rhod(this->mem->advectee(ix::th).shape()); // TODO: replace all rhod arrays with this->mem->G
    rhod = (*params.rhod)(this->vert_idx);

    // temporary array of pressure - prtcls cant be init'd with 1D profile
    typename parent_t::arr_t p_e(this->mem->advectee(ix::th).shape()); 
    p_e = (*params.p_e)(this->vert_idx);

    prtcls->init(
      make_arrinfo(this->mem->advectee(ix::th)),
      make_arrinfo(this->mem->advectee(ix::rv)),
      make_arrinfo(rhod)
      ,make_arrinfo(p_e)
    ); 
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
    this->record_aux_const("n_sd_max", params.cloudph_opts_init.n_sd_max);  
    this->record_aux_const("adve_scheme", params.cloudph_opts_init.adve_scheme);  
    this->record_aux_const("dev_count", params.cloudph_opts_init.dev_count);  
    this->record_aux_const("dev_id", params.cloudph_opts_init.dev_id);  
    this->record_aux_const("sstp_cond", params.cloudph_opts_init.sstp_cond);  
    this->record_aux_const("sstp_coal", params.cloudph_opts_init.sstp_coal);  
    this->record_aux_const("sstp_chem", params.cloudph_opts_init.sstp_chem);  
    this->record_aux_const("exact_sstp_cond", params.cloudph_opts_init.exact_sstp_cond);  
    this->record_aux_const("diag_incloud_time", params.cloudph_opts_init.diag_incloud_time);  
    this->record_aux_const("rng_seed", params.cloudph_opts_init.rng_seed);  
    this->record_aux_const("kernel", params.cloudph_opts_init.kernel);  
    this->record_aux_const("terminal_velocity", params.cloudph_opts_init.terminal_velocity);  
    this->record_aux_const("RH_formula", params.cloudph_opts_init.RH_formula);  
    this->record_aux_const("backend", params.backend);  
    this->record_aux_const("async", params.async);  
    this->record_aux_const("adve", params.cloudph_opts.adve);  
    this->record_aux_const("sedi", params.cloudph_opts.sedi);  
    this->record_aux_const("subs", params.cloudph_opts.subs);  
    this->record_aux_const("cond", params.cloudph_opts.cond);  
    this->record_aux_const("coal", params.flag_coal);  // cloudph_opts.coal could be 0 here due to spinup
    this->record_aux_const("rcyc", params.cloudph_opts.rcyc);  
    this->record_aux_const("out_dry_spec", params.out_dry_spec);  
    this->record_aux_const("out_wet_spec", params.out_wet_spec);  
    this->record_aux_const("gccn", params.gccn);  
    this->record_aux_const("turb_adve", params.cloudph_opts.turb_adve);  
    this->record_aux_const("turb_cond", params.cloudph_opts.turb_cond);  
    this->record_aux_const("turb_coal", params.cloudph_opts.turb_coal);  
    this->record_aux_const("chem_switch", params.cloudph_opts_init.chem_switch);  
    this->record_aux_const("coal_switch", params.cloudph_opts_init.coal_switch);  
    this->record_aux_const("sedi_switch", params.cloudph_opts_init.sedi_switch);  
    this->record_aux_const("subs_switch", params.cloudph_opts_init.subs_switch);  
    this->record_aux_const("src_switch", params.cloudph_opts_init.src_switch);  
    this->record_aux_const("turb_adve_switch", params.cloudph_opts_init.turb_adve_switch);  
    this->record_aux_const("turb_cond_switch", params.cloudph_opts_init.turb_cond_switch);  
    this->record_aux_const("turb_coal_switch", params.cloudph_opts_init.turb_coal_switch);  
    this->record_aux_const("chem_dsl", params.cloudph_opts.chem_dsl);  
    this->record_aux_const("chem_dsc", params.cloudph_opts.chem_dsc);  
    this->record_aux_const("chem_rct", params.cloudph_opts.chem_rct);  
    this->record_aux_const("chem_rho", params.cloudph_opts_init.chem_rho);  
    this->record_aux_const("opts_init RH_max", params.cloudph_opts_init.RH_max);  
    this->record_aux_const("supstp_src", params.cloudph_opts_init.supstp_src);  
    this->record_aux_const("src_sd_conc", params.cloudph_opts_init.src_sd_conc);  
    this->record_aux_const("src_z1", params.cloudph_opts_init.src_z1);
  }
  this->mem->barrier();
}

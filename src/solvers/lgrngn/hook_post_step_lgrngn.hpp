#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_post_step()
{
  parent_t::hook_post_step(); 
  
  // storej positions for precoal diag
  if (this->rank == 0)
  {
    if(params.user_params.precoal_start >= 0 && (this->timestep >= params.user_params.precoal_start) && (int(this->timestep - params.user_params.precoal_start) % int(params.user_params.precoal_outfreq) == 0))
      prtcls->store_ijk(this->timestep * params.user_params.dt);
  }
  this->mem->barrier();

  // write mean precoal distance at the end of the simulation to a file
  // TODO: won't work with MPI
  if (this->rank == 0)
  {
    if(params.user_params.precoal_start >= 0 && this->timestep == params.user_params.nt) 
    {
      std::ofstream os(this->outdir+"/precoal_mean_time.dat");
      os << "# time before collision [s] | mean distance between droplets [m]" << std::endl;
      std::vector<real_t> precoal_distance = prtcls->diag_precoal_distance();
      for(int i=0; i<precoal_distance.size(); ++i)
      {
        os << (i+0.5) * (params.cloudph_opts_init.precoal_stats_tmax / 100.) << " " << precoal_distance.at(i) << std::endl; // assuming libcloud config.precoal_stats_bins=100
      }
    //  os << std::endl;
    }
  }
  this->mem->barrier();
}

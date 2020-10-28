#pragma once
#include "../slvr_lgrngn.hpp"
#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_ante_step()
{
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
#if defined(UWLCM_TIMING)
      tbeg = parent_t::clock::now();
#endif
      ftr.get();
#if defined(UWLCM_TIMING)
      tend = parent_t::clock::now();
      parent_t::tasync_wait += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
#endif
    } else assert(!ftr.valid()); 
#endif
  }
  this->mem->barrier();
  parent_t::hook_ante_step(); // includes RHS, which in turn launches sync_in and step_cond
  negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after at the end of hook_ante_step");
}

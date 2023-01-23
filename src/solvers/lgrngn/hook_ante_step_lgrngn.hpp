#pragma once
#include "../slvr_lgrngn.hpp"
#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_ante_step()
{
  parent_t::hook_ante_step(); // includes RHS, which in turn launches sync_in and step_cond
  negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after at the end of hook_ante_step");
}

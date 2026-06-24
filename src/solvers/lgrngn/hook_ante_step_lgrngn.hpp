#pragma once
#include "../slvr_lgrngn.hpp"
#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif


/**
 * @brief Performs tasks before each simulation timestep in the Lagrangian microphysics solver.
 *
 * @details
 * This function is called at the beginning of each timestep. It executes the parent class
 * hook. It performs a sanity check to ensure
 * that the water vapor field (`rv`) has no negative values.
 */
template <class ct_params_t>
void slvr_lgrngn<ct_params_t>::hook_ante_step()
{
  parent_t::hook_ante_step(); // includes RHS, which in turn launches sync_in and step_cond
  negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after at the end of hook_ante_step");
}

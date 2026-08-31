#pragma once
#include "../slvr_lgrngn.hpp"
// #include <iostream>
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
  // if (this->rank == 0)// && this->timestep == 150)
  // {
  //   const auto T_freeze = prtcls->get_attr("T_freeze");
  //   // for (const auto value : T_freeze)
  //     // std::cout << value << std::endl;
  //   if (!T_freeze.empty())
  //     std::cout << *std::max_element(T_freeze.begin(), T_freeze.end()) << std::endl;

  //   const auto rd2_insol = prtcls->get_attr("rd2_insol");
  //   // for (const auto value : rd2_insol)
  //     // std::cout << value << std::endl;
  //   if (!rd2_insol.empty())
  //     std::cout << *std::max_element(rd2_insol.begin(), rd2_insol.end()) << std::endl;
  // }
  negcheck(this->mem->advectee(ix::rv)(this->ijk), "rv after at the end of hook_ante_step");
}

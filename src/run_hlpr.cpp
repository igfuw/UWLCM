/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "detail/setup.hpp"
#include "detail/concurr_types.hpp"
#include "detail/ct_params.hpp"

#include "cases/detail/api_test.hpp"
#include "cases/DYCOMS.hpp"
#include "cases/RICO11.hpp"
#include "cases/MoistThermalGrabowskiClark99.hpp"
#include "cases/DryThermalGMD2015.hpp"
#include "cases/LasherTrapp2001.hpp"

#include "opts/opts_lgrngn.hpp"
#include "opts/opts_blk_1m.hpp"
#include "opts/opts_blk_2m.hpp"

#include "detail/panic.hpp"

#include "solvers/slvr_lgrngn.hpp"
#include "solvers/slvr_blk_1m.hpp"
#include "solvers/slvr_blk_2m.hpp"

#include "forcings/calc_forces_blk_1m.hpp"
#include "forcings/calc_forces_blk_2m.hpp"
#include "forcings/calc_forces_common.hpp"

#include "run_hlpr.hpp"

// dimension-independent model run logic - the same for any microphysics
template <class solver_t, int n_dims>
void run(const int (&nps)[n_dims], const user_params_t &user_params)
{
  auto nz = nps[n_dims - 1];
  
  using concurr_any_t = libmpdataxx::concurr::any<
    typename solver_t::real_t, 
    n_dims
  >;

  using concurr_openmp_cyclic_t = typename concurr_openmp_cyclic<solver_t, n_dims>::type;
  using concurr_openmp_rigid_t = typename concurr_openmp_rigid<solver_t, n_dims>::type;
  using concurr_openmp_cyclic_gndsky_t = typename concurr_openmp_cyclic_gndsky<solver_t, n_dims>::type;
  using concurr_openmp_rigid_gndsky_t = typename concurr_openmp_rigid_gndsky<solver_t, n_dims>::type;
  
  using rt_params_t = typename solver_t::rt_params_t;
  using ix = typename solver_t::ix;
  
  struct case_ct_params_t
  {
    using rt_params_t = typename solver_t::rt_params_t;
    using ix = typename solver_t::ix;
    enum {enable_sgs = solver_t::ct_params_t_::sgs_scheme != libmpdataxx::solvers::iles};
  };

  using case_t = setup::CasesCommon<
    case_ct_params_t, n_dims
  >;

  using case_ptr_t = std::unique_ptr<
    case_t
  >;

  case_ptr_t case_ptr; 

  // setup choice
  if (user_params.model_case == "moist_thermal")
    case_ptr.reset(new setup::moist_thermal::MoistThermalGrabowskiClark99<case_ct_params_t, n_dims>()); 
  else if (user_params.model_case == "dry_thermal")
    case_ptr.reset(new setup::dry_thermal::DryThermal<case_ct_params_t, n_dims>()); 
  else if (user_params.model_case == "dycoms_rf01")
    case_ptr.reset(new setup::dycoms::Dycoms<case_ct_params_t, 1, n_dims>()); 
  else if (user_params.model_case == "dycoms_rf02")
    case_ptr.reset(new setup::dycoms::Dycoms<case_ct_params_t, 2, n_dims>()); 
  else if (user_params.model_case == "lasher_trapp")
    case_ptr.reset(new setup::LasherTrapp::LasherTrapp2001<case_ct_params_t, n_dims>());
  else if (user_params.model_case == "rico11")
    case_ptr.reset(new setup::rico::Rico11<case_ct_params_t, n_dims>());
  // special versions for api test - they have much lower aerosol concentrations to avoid multiplicity overflow
  else if (user_params.model_case == "moist_thermal_api_test")
    case_ptr.reset(new setup::api_test<setup::moist_thermal::MoistThermalGrabowskiClark99<case_ct_params_t, n_dims>>()); 
  else if (user_params.model_case == "dry_thermal_api_test")
    case_ptr.reset(new setup::api_test<setup::dry_thermal::DryThermal<case_ct_params_t, n_dims>>()); 
  else if (user_params.model_case == "dycoms_rf01_api_test")
    case_ptr.reset(new setup::api_test<setup::dycoms::Dycoms<case_ct_params_t, 1, n_dims>>()); 
  else if (user_params.model_case == "dycoms_rf02_api_test")
    case_ptr.reset(new setup::api_test<setup::dycoms::Dycoms<case_ct_params_t, 2, n_dims>>()); 
  else if (user_params.model_case == "lasher_trapp_api_test")
    case_ptr.reset(new setup::api_test<setup::LasherTrapp::LasherTrapp2001<case_ct_params_t, n_dims>>());
  else if (user_params.model_case == "rico11_api_test")
    case_ptr.reset(new setup::api_test<setup::rico::Rico11<case_ct_params_t, n_dims>>());
  else
    throw std::runtime_error("wrong case choice");

  // instantiation of structure containing simulation parameters
  rt_params_t p;

  // copy force constants
  p.ForceParameters = case_ptr->ForceParameters;

  // copy functions used to update surface fluxes
  p.update_surf_flux_sens = std::bind(&case_t::update_surf_flux_sens, case_ptr.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
  p.update_surf_flux_lat  = std::bind(&case_t::update_surf_flux_lat,  case_ptr.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);
  p.update_surf_flux_uv   = std::bind(&case_t::update_surf_flux_uv,   case_ptr.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7, std::placeholders::_8);

  // copy user_params for output
  p.user_params = user_params;

  // output and simulation parameters
  for (int d = 0; d < n_dims; ++d)
  {
    p.grid_size[d] = nps[d];
  }

  case_ptr->setopts(p, nps, user_params);

  // reference profiles shared among threads
  setup::profiles_t profs(nz); 
  // rhod needs to be bigger, cause it divides vertical courant number, TODO: should have a halo both up and down, not only up like now; then it should be interpolated in courant calculation

  // assign their values
  case_ptr->set_profs(profs, nz, user_params);
  // pass them to rt_params
  setup::copy_profiles(profs, p);

  // set micro-specific options, needs to be done after copy_profiles
  setopts_micro<solver_t>(p, user_params, case_ptr);

  // set outvars
  p.outvars.insert({ix::rv, {"rv", "[kg kg-1]"}});
  p.outvars.insert({ix::th, {"th", "[K]"}});
  p.outvars.insert({ix::u, {"u", "[m/s]"}});
  p.outvars.insert({ix::w, {"w", "[m/s]"}});
  if (n_dims > 2)
  {
    // WARNING: assumes certain ordering of variables to avoid tedious template programming !
    p.outvars.insert({1, {"v", "[m/s]"}});
  }

  // solver instantiation
  std::unique_ptr<concurr_any_t> concurr;

  if(user_params.model_case == "dry_thermal" || user_params.model_case == "dry_thermal_api_test")
  {
    concurr.reset(new concurr_openmp_cyclic_t(p));
  }
  else if(user_params.model_case == "lasher_trapp" || user_params.model_case == "lasher_trapp_api_test")
  {
    //concurr.reset(new concurr_openmp_rigid_gndsky_t(p));     // rigid horizontal boundaries
    concurr.reset(new concurr_openmp_cyclic_gndsky_t(p)); // cyclic horizontal boundaries, as in the ICMW2020 case
  }
  else
  {
    concurr.reset(new concurr_openmp_cyclic_gndsky_t(p));
  }
  
  case_ptr->intcond(*concurr.get(), profs.rhod, profs.th_e, profs.rv_e, profs.rl_e, profs.p_e, user_params.rng_seed);

  // setup panic pointer and the signal handler
  panic = concurr->panic_ptr();
  set_sigaction();
 
  // timestepping
  concurr->advance(user_params.nt);
}

template<template<class...> class slvr, class ct_params_dim_micro, int n_dims>
void run_hlpr(bool piggy, bool sgs, const std::string &type, const int (&nps)[n_dims], const user_params_t &user_params)
{
  if(piggy) // piggybacker
  {
#if !defined(UWLCM_DISABLE_PIGGYBACKER)
    struct ct_params_piggy : ct_params_dim_micro { enum { piggy = 1 }; };
    run<slvr<ct_params_piggy>>(nps, user_params);
#endif
  }
  else // driver
  {
#if !defined(UWLCM_DISABLE_DRIVER)
    struct ct_params_piggy : ct_params_dim_micro { enum { piggy = 0 }; };
    if (sgs)
    {
  #if !defined(UWLCM_DISABLE_SGS)
      struct ct_params_sgs : ct_params_piggy
      {
        enum { sgs_scheme = libmpdataxx::solvers::smg };
        enum { stress_diff = libmpdataxx::solvers::compact };
      };
      run<slvr<ct_params_sgs>>(nps, user_params);
  #endif
    }
    else
    {
  #if !defined(UWLCM_DISABLE_ILES)
      struct ct_params_iles : ct_params_piggy {};
      run<slvr<ct_params_iles>>(nps, user_params);
  #endif
    }
#endif
  }
}

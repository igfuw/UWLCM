/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "detail/setup.hpp"
#include "detail/concurr_types.hpp"
#include "detail/panic.hpp"

#include "cases/DYCOMS.hpp"
#include "cases/RICO11.hpp"
#include "cases/MoistThermalGrabowskiClark99.hpp"
#include "cases/DryThermalGMD2015.hpp"
#include "cases/CumulusCongestus.hpp"
#include "cases/DryPBL.hpp"

#if defined(UWLCM_TIMING)
  #include "detail/exec_timer.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_LGRNGN) || !defined(UWLCM_DISABLE_3D_LGRNGN)
  #include "opts/opts_lgrngn.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_BLK_1M) || !defined(UWLCM_DISABLE_3D_BLK_1M)
  #include "opts/opts_blk_1m.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_BLK_2M) || !defined(UWLCM_DISABLE_3D_BLK_2M)
  #include "opts/opts_blk_2m.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_NONE) || !defined(UWLCM_DISABLE_3D_NONE)
  #include "opts/opts_dry.hpp"
#endif

#include "common_headers.hpp"

// dimension-independent model run logic - the same for any microphysics
template <class solver_t, int n_dims>
void run(const int (&nps)[n_dims], const user_params_t &user_params)
{
  const int nz = nps[n_dims - 1],
            nz_ref = (nz - 1) * pow(2, user_params.n_fra_iter) + 1; // TODO: use grid_size_ref or refine_grid_size from libmpdata++
  
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

  using case_t = cases::CasesCommon<
    case_ct_params_t, n_dims
  >;

  using case_ptr_t = std::unique_ptr<
    case_t
  >;

  case_ptr_t case_ptr; 

  // setup choice
  if (user_params.model_case == "moist_thermal")
    case_ptr.reset(new cases::moist_thermal::MoistThermalGrabowskiClark99<case_ct_params_t, n_dims>(user_params.X, user_params.Y, user_params.Z)); 
  else if (user_params.model_case == "dry_thermal")
    case_ptr.reset(new cases::dry_thermal::DryThermal<case_ct_params_t, n_dims>(user_params.X, user_params.Y, user_params.Z)); 
  else if (user_params.model_case == "dycoms_rf01")
    case_ptr.reset(new cases::dycoms::Dycoms<case_ct_params_t, 1, n_dims>(user_params.X, user_params.Y, user_params.Z)); 
  else if (user_params.model_case == "dycoms_rf02")
    case_ptr.reset(new cases::dycoms::Dycoms<case_ct_params_t, 2, n_dims>(user_params.X, user_params.Y, user_params.Z)); 
  else if (user_params.model_case == "cumulus_congestus")
    case_ptr.reset(new cases::CumulusCongestus::CumulusCongestus<case_ct_params_t, n_dims>(user_params.X, user_params.Y, user_params.Z));
  else if (user_params.model_case == "rico11")
    case_ptr.reset(new cases::rico::Rico11<case_ct_params_t, n_dims>(user_params.X, user_params.Y, user_params.Z));
  else if (user_params.model_case == "dry_pbl")
    case_ptr.reset(new cases::pbl::DryPBL<case_ct_params_t, n_dims>(user_params.X, user_params.Y, user_params.Z));
  else
    throw std::runtime_error("wrong case choice");

  // instantiation of structure containing simulation parameters
  rt_params_t p;

  // copy force constants
  p.ForceParameters = case_ptr->ForceParameters;

  // copy user_params
  p.user_params = user_params;

  // some runtime parameters defined in libmpdata++ are passed via user_params
  p.outdir = user_params.outdir;
  p.outfreq = user_params.outfreq;
  p.dt = user_params.dt;
  p.n_fra_iter = user_params.n_fra_iter;

  // output and simulation parameters
  for (int d = 0; d < n_dims; ++d)
  {
    p.grid_size[d] = nps[d];
  }

  case_ptr->setopts(p, nps, user_params);

  // reference profiles (on normal and refined grids). 
  // They are shared among threads! When params are copied by each thread's solver, 
  // blitz::array copy constructor gives reference to the same data
  p.profs.init(nz);
  p.profs_ref.init(nz_ref);
  // NOTE for Anelastic: env profiles on both grids agree very well
  //                     reference profiles differ a little, e.g. for rico11 nz=51 there's up to 0.02% difference in rhod for n_fra_iter=1 (and 0.03% for n_fra_iter=2)
  //                     but thats acceptable (?) (rhod doesnt affect RH)
  case_ptr->set_profs(p.profs, nz, user_params);
  case_ptr->set_profs(p.profs_ref, nz_ref, user_params);
  // pass them to rt_params
//  p.profs     = profs;
//  p.profs_ref = profs_ref;
//  detail::copy_profiles(profs, p.profs);
//  detail::copy_profiles(profs_ref, p);

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


  // copy functions used to update surface fluxes
  // NOTE: some parameters (di, dj, dz for now) are passed by value, hence this needs to be done after their values are set
  //       and also means that these parameters cannot change during simulation
  p.update_surf_flux_sens = std::bind(&case_t::update_surf_flux_sens, case_ptr.get(), 
		                      std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 
				      std::placeholders::_4, std::placeholders::_5, p.di, p.dj, p.dz / 2);
  p.update_surf_flux_lat  = std::bind(&case_t::update_surf_flux_lat , case_ptr.get(), 
		                      std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 
				      std::placeholders::_4, std::placeholders::_5, p.di, p.dj, p.dz / 2);

  p.update_surf_flux_uv   = std::bind(&case_t::update_surf_flux_uv,   case_ptr.get(), 
		                      std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, 
				      std::placeholders::_4, std::placeholders::_5, p.di, p.dj, p.dz / 2,
				      std::placeholders::_6);

  // copy functions used to update large-scale forcings
  p.update_rv_LS = std::bind(&case_t::update_rv_LS, case_ptr.get(), std::ref(p.profs.rv_LS), std::placeholders::_1, std::placeholders::_2, p.dz);
  p.update_th_LS = std::bind(&case_t::update_th_LS, case_ptr.get(), std::ref(p.profs.th_LS), std::placeholders::_1, std::placeholders::_2, p.dz);

  // solver instantiation
  std::unique_ptr<concurr_any_t> concurr;

  if(user_params.model_case == "dry_thermal")
  {
    concurr.reset(new concurr_openmp_cyclic_t(p));
  }
  else if(user_params.model_case == "cumulus_congestus")
  {
    //concurr.reset(new concurr_openmp_rigid_gndsky_t(p));     // rigid horizontal boundaries
    concurr.reset(new concurr_openmp_cyclic_gndsky_t(p)); // cyclic horizontal boundaries, as in the ICMW2020 case
  }
  else
  {
    concurr.reset(new concurr_openmp_cyclic_gndsky_t(p));
  }
  
  case_ptr->intcond(*concurr.get(), p.profs.rhod, p.profs.th_e, p.profs.rv_e, p.profs.rl_e, p.profs.p_e, user_params.rng_seed_init);

  // setup panic pointer and the signal handler
  panic = concurr->panic_ptr();
  set_sigaction();
 
  // timestepping
  concurr->advance(user_params.nt);
}

#if defined(UWLCM_TIMING)
  template<class slvr>
  using timer = exec_timer<slvr>;
#else
  template<class slvr>
  using timer = slvr;
#endif


template<template<class...> class slvr, class ct_params_dim_micro, int n_dims>
void run_hlpr(bool piggy, bool sgs, const std::string &type, const int (&nps)[n_dims], const user_params_t &user_params)
{
  if(piggy) // piggybacker
  {
#if !defined(UWLCM_DISABLE_PIGGYBACKER)
    struct ct_params_piggy : ct_params_dim_micro { enum { piggy = 1 }; };
    run<timer<slvr<ct_params_piggy>>>(nps, user_params);
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
      run<timer<slvr<ct_params_sgs>>>(nps, user_params);
  #endif
    }
    else
    {
  #if !defined(UWLCM_DISABLE_ILES)
      struct ct_params_iles : ct_params_piggy {};
      run<timer<slvr<ct_params_iles>>>(nps, user_params);
  #endif
    }
#endif
  }
}

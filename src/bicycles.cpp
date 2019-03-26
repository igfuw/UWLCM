/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "setup.hpp"
#include "concurr_types.hpp"

#include "cases/DYCOMS.hpp"
#include "cases/MoistThermalGrabowskiClark99.hpp"
#include "cases/DryThermalGMD2015.hpp"
#include "cases/LasherTrapp2001.hpp"

#include "opts_lgrngn.hpp"
#include "opts_blk_1m.hpp"

#include "panic.hpp"
#include <map>

// dimension-independent model run logic - the same for any microphysics
template <class solver_t, int n_dims>
void run(const int (&nps)[n_dims], const user_params_t &user_params)
{
  auto nz = nps[n_dims - 1];
  
  using concurr_any_t = concurr::any<
    typename solver_t::real_t, 
    n_dims
  >;

  using concurr_openmp_cyclic_t = typename concurr_openmp_cyclic<solver_t, n_dims>::type;
  using concurr_openmp_rigid_t = typename concurr_openmp_rigid<solver_t, n_dims>::type;
  using concurr_openmp_cyclic_rigid_t = typename concurr_openmp_cyclic_rigid<solver_t, n_dims>::type;
  
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
  else
    throw std::runtime_error("wrong setup choice");

  // instantiation of structure containing simulation parameters
  rt_params_t p;

  // copy force constants
  p.ForceParameters = case_ptr->ForceParameters;

  // copy functions used to update surface fluxes
  p.update_surf_flux_sens = std::bind(&case_t::update_surf_flux_sens, case_ptr.get(), std::placeholders::_1,std::placeholders:: _2, std::placeholders::_3);
  p.update_surf_flux_lat = std::bind(&case_t::update_surf_flux_lat, case_ptr.get(), std::placeholders::_1,std::placeholders:: _2, std::placeholders::_3);

  // copy user_params for output
  p.user_params = user_params;

  // output and simulation parameters
  for (int d = 0; d < n_dims; ++d)
  {
    p.grid_size[d] = nps[d];
  }

  case_ptr->setopts(p, nps, user_params);
  setopts_micro<solver_t>(p, user_params, case_ptr);

  // reference profiles shared among threads
  setup::profiles_t profs(nz); 
  // rhod needs to be bigger, cause it divides vertical courant number, TODO: should have a halo both up and down, not only up like now; then it should be interpolated in courant calculation

  // assign their values
  case_ptr->env_prof(profs, nz, user_params);
  // pass them to rt_params
  setup::copy_profiles(profs, p);

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

  if(user_params.model_case == "dry_thermal")
  {
    concurr.reset(new concurr_openmp_cyclic_t(p));
  }
  else if(user_params.model_case == "lasher_trapp")
  {
    concurr.reset(new concurr_openmp_rigid_t(p));
  }
  else
  {
    concurr.reset(new concurr_openmp_cyclic_rigid_t(p));
  }
  
  case_ptr->intcond(*concurr.get(), profs.rhod, profs.th_e, profs.rv_e, profs.rl_e, profs.p_e, user_params.rng_seed);

  // setup panic pointer and the signal handler
  panic = concurr->panic_ptr();
  set_sigaction();
 
  // timestepping
  concurr->advance(user_params.nt);
}

// libmpdata++'s compile-time parameters
struct ct_params_common : ct_params_default_t
{
  using real_t = setup::real_t;
  enum { rhs_scheme = solvers::mixed }; 
  enum { prs_scheme = solvers::cr };
  enum { vip_vab = solvers::expl };
};

struct ct_params_2D_sd : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 4 };
  struct ix { enum {
    u, w, th, rv, 
    vip_i=u, vip_j=w, vip_den=-1
  }; };
  enum { delayed_step = opts::bit(ix::th) | opts::bit(ix::rv) };
};

struct ct_params_3D_sd : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 5 };
  struct ix { enum {
    u, v, w, th, rv, 
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
  enum { delayed_step = opts::bit(ix::th) | opts::bit(ix::rv) };
};

struct ct_params_2D_blk_1m : ct_params_common
{
  enum { n_dims = 2 };
  enum { n_eqns = 6 };
  struct ix { enum {
    u, w, th, rv, rc, rr,
    vip_i=u, vip_j=w, vip_den=-1
  }; };
};

struct ct_params_3D_blk_1m : ct_params_common
{
  enum { n_dims = 3 };
  enum { n_eqns = 7 };
  struct ix { enum {
    u, v, w, th, rv, rc, rr, 
    vip_i=u, vip_j=v, vip_k=w, vip_den=-1
  }; };
};

// function used to modify ct_params before running
template<template<class...> class slvr, class ct_params_dim_micro, int n_dims>
void run_hlpr(bool piggy, bool sgs, const std::string &type, const int (&nps)[n_dims], const user_params_t &user_params)
{
  if(!piggy) // no piggybacking
  {
    struct ct_params_piggy : ct_params_dim_micro { enum { piggy = 0 }; };
    if(type == "moist_thermal") // use abs option in moist_thermal
    {
      struct ct_params_final : ct_params_piggy { enum { opts = opts::nug | opts::abs }; };
      run<slvr<ct_params_final>>(nps, user_params);
    }
    else // default is the iga | fct option
    {
      struct ct_params_opts : ct_params_piggy { enum { opts = opts::nug | opts::iga | opts::fct }; };
      if (sgs)
      {
        struct ct_params_final : ct_params_opts
        {
          enum { sgs_scheme = solvers::smg };
          enum { stress_diff = solvers::compact };
        };
        run<slvr<ct_params_final>>(nps, user_params);
      }
      else
      {
        struct ct_params_final : ct_params_opts {};
        run<slvr<ct_params_final>>(nps, user_params);
      }
    }
  }
  else // piggybacking
  {
    struct ct_params_piggy : ct_params_dim_micro { enum { piggy = 1 }; };
    if(type == "moist_thermal") // use abs option in moist_thermal
    {
      struct ct_params_final : ct_params_piggy { enum { opts = opts::nug | opts::abs }; };
      run<slvr<ct_params_final>>(nps, user_params);
    }
    else // default is the iga | fct option
    {
      struct ct_params_final : ct_params_piggy { enum { opts = opts::nug | opts::iga | opts::fct }; };
      run<slvr<ct_params_final>>(nps, user_params);
    }
  }
}


// all starts here with handling general options 
int main(int argc, char** argv)
{
  omp_set_nested(1); // to allow openmp calls from libcloudphxx multi_CUDA backend
  // making argc and argv global
  ac = argc;
  av = argv;

  {
    // note: all options should have default values here to make "--micro=? --help" work
    opts_main.add_options()
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("case", po::value<std::string>()->required(), "one of: dry_thermal, moist_thermal, dycoms, lasher_trapp")
      ("nx", po::value<int>()->default_value(76) , "grid cell count in horizontal")
      ("ny", po::value<int>()->default_value(0) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(76) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(3600) , "timestep count")
      ("rng_seed", po::value<int>()->default_value(-1) , "rng seed, negative for random")
      ("dt", po::value<setup::real_t>()->required() , "timestep length")
      ("z_rlx_sclr", po::value<setup::real_t>()->default_value(25) , "scalars surface flux charasteristic heihjt")
      ("outdir", po::value<std::string>(), "output file name (netCDF-compatible HDF5)")
      ("outfreq", po::value<int>(), "output rate (timestep interval)")
      ("spinup", po::value<int>()->default_value(2400) , "number of initial timesteps during which rain formation is to be turned off")
      ("serial", po::value<bool>()->default_value(false), "force advection and microphysics to be computed on single thread")
      ("th_src", po::value<bool>()->default_value(true) , "temp src")
      ("rv_src", po::value<bool>()->default_value(true) , "water vap source")
      ("rc_src", po::value<bool>()->default_value(true) , "cloud water source (in blk_1m)")
      ("rr_src", po::value<bool>()->default_value(true) , "rain water source (in blk_1m)")
      ("uv_src", po::value<bool>()->default_value(true) , "horizontal vel src")
      ("w_src", po::value<bool>()->default_value(true) , "vertical vel src")
      ("piggy", po::value<bool>()->default_value(false) , "is it a piggybacking run")
      ("sgs", po::value<bool>()->default_value(false) , "is subgrid-scale turbulence model on")
      ("sgs_delta", po::value<setup::real_t>()->default_value(-1) , "subgrid-scale turbulence model delta")
      ("help", "produce a help message (see also --micro X --help)")
    ;
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); // ignores unknown

    // hendling the "help" option
    if (ac == 1 || (vm.count("help") && !vm.count("micro"))) 
    {
      std::cout << opts_main;
      exit(EXIT_SUCCESS);
    }

    // checking if all required options present
    po::notify(vm); 

    // instantiating user params container
    user_params_t user_params;
    
    if (!vm.count("help"))
    {
      if (!vm.count("outdir")) throw po::required_option("outdir");
      if (!vm.count("outfreq")) throw po::required_option("outfreq");
      user_params.outdir = vm["outdir"].as<std::string>();
      user_params.outfreq = vm["outfreq"].as<int>();
    }

    int 
      nx = vm["nx"].as<int>(),
      ny = vm["ny"].as<int>(),
      nz = vm["nz"].as<int>();

    user_params.nt = vm["nt"].as<int>(),
    user_params.spinup = vm["spinup"].as<int>();
 
    // handling rng_seed
    user_params.rng_seed = vm["rng_seed"].as<int>();
    if(user_params.rng_seed < 0) //if negative, get random seed
    {
      std::random_device rd; 
      user_params.rng_seed = rd();
    }
    std::cout << "rng seed: " << user_params.rng_seed << std::endl;
   
    //handling timestep length
    user_params.dt = vm["dt"].as<setup::real_t>();

    //handling z_rlx_sclr
    user_params.z_rlx_sclr = vm["z_rlx_sclr"].as<setup::real_t>();

    // handling serial-advection-forcing flag
    if(vm["serial"].as<bool>()) setenv("OMP_NUM_THREADS", "1", 1);

    // handling sources flags
    user_params.th_src = vm["th_src"].as<bool>();
    user_params.rv_src = vm["rv_src"].as<bool>();
    user_params.rc_src = vm["rc_src"].as<bool>();
    user_params.rr_src = vm["rr_src"].as<bool>();
    user_params.uv_src = vm["uv_src"].as<bool>();
    user_params.w_src = vm["w_src"].as<bool>();

    bool piggy = vm["piggy"].as<bool>();
    bool sgs = vm["sgs"].as<bool>();
    user_params.sgs_delta = vm["sgs_delta"].as<setup::real_t>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();

    // handling the "case" option
    user_params.model_case = vm["case"].as<std::string>();

    // run the simulation
    if (micro == "lgrngn" && ny == 0) // 2D super-droplet
      run_hlpr<slvr_lgrngn, ct_params_2D_sd>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);

    else if (micro == "lgrngn" && ny > 0) // 3D super-droplet
      run_hlpr<slvr_lgrngn, ct_params_3D_sd>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);

    else if (micro == "blk_1m" && ny == 0) // 2D one-moment
      run_hlpr<slvr_blk_1m, ct_params_2D_blk_1m>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);

    else if (micro == "blk_1m" && ny > 0) // 3D one-moment
      run_hlpr<slvr_blk_1m, ct_params_3D_blk_1m>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);

  // TODO: not only micro can be wrong
    else throw 
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro" 
      );
  }
}

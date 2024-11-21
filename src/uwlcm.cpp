/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <iostream>
#include <random>

#ifdef _OPENMP
# include <omp.h>
#endif

#include "opts/opts_common.hpp"

#include "run_hlpr.hpp"

#include "detail/ct_params.hpp"

#if !defined(UWLCM_DISABLE_2D_LGRNGN) || !defined(UWLCM_DISABLE_3D_LGRNGN)
//  #include "opts/opts_lgrngn.hpp"
  #include "solvers/slvr_lgrngn.hpp"
  #include "solvers/lgrngn/diag_lgrngn.hpp" 
  #include "solvers/lgrngn/hook_ante_delayed_step_lgrngn.hpp"
  #include "solvers/lgrngn/hook_ante_loop_lgrngn.hpp"
  #include "solvers/lgrngn/hook_ante_step_lgrngn.hpp" 
  #include "solvers/lgrngn/hook_mixed_rhs_ante_step_lgrngn.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_BLK_1M) || !defined(UWLCM_DISABLE_3D_BLK_1M)
//  #include "opts/opts_blk_1m.hpp"
  #include "solvers/slvr_blk_1m.hpp"
  #include "solvers/blk_1m/calc_forces_blk_1m_common.hpp"
  #include "solvers/blk_1m/update_rhs_blk_1m_common.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_BLK_2M) || !defined(UWLCM_DISABLE_3D_BLK_2M)
//  #include "opts/opts_blk_2m.hpp"
  #include "solvers/slvr_blk_2m.hpp"
  #include "solvers/blk_2m/calc_forces_blk_2m_common.hpp"
  #include "solvers/blk_2m/update_rhs_blk_2m_common.hpp"
#endif

#if !defined(UWLCM_DISABLE_2D_NONE) || !defined(UWLCM_DISABLE_3D_NONE)
  #include "solvers/slvr_dry.hpp"
#endif

#include "solvers/common/calc_forces_common.hpp"

#include <map>

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
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn, none")
      ("case", po::value<std::string>()->required(), "one of: dry_thermal, moist_thermal, dycoms_rf01, dycoms_rf02, cumulus_congestus_icmw20, cumulus_congestus_icmw24, rico11, dry_pbl")
      ("X", po::value<setup::real_t>()->default_value(-1) , "domain size in X [m] (set negative for case default)")
      ("Y", po::value<setup::real_t>()->default_value(-1) , "domain size in Y [m] (set negative for case default)")
      ("Z", po::value<setup::real_t>()->default_value(-1) , "domain size in Z [m] (set negative for case default)")
      ("nx", po::value<int>()->required() , "grid cell count in horizontal")
      ("ny", po::value<int>()->default_value(0) , "grid cell count in horizontal, 0 for 2D simulation")
      ("nz", po::value<int>()->required() , "grid cell count in vertical")
      ("nt", po::value<int>()->required() , "timestep count")
      ("rng_seed", po::value<int>()->default_value(0) , "rng seed for randomness post initialization (currently only in Lagrangian microphysics), 0 for random")
      ("rng_seed_init", po::value<int>()->default_value(0) , "rng seed for initial conditions (perturbations of th and rv and initialization of Lagrangian microphysics), 0 for rng_seed_init=rng_seed")
      ("dt", po::value<setup::real_t>()->required() , "timestep length [s]")
      ("outdir", po::value<std::string>()->required(), "output directory name (netCDF-compatible HDF5)")
      ("outfreq", po::value<int>()->required(), "output rate (timestep interval)")
      ("outstart", po::value<int>()->default_value(0), "output starts after this many timesteps")
      ("outwindow", po::value<int>()->default_value(1), "number of consecutive timesteps output is done, starts at outfreq (doesnt affect output of droplet spectra from lagrangian microphysics)")
      ("spinup", po::value<int>()->default_value(0) , "number of initial timesteps during which rain formation is to be turned off")
      ("serial", po::value<bool>()->default_value(false), "force CPU component of the model (dynamics and bulk microphysics) to be computed on single thread")
      ("window", po::value<bool>()->default_value(false), "moving-window simulation, i.e. mean horizontal velocity substracted from advectors")
//      ("th_src", po::value<bool>()->default_value(true) , "temp src")
//      ("rv_src", po::value<bool>()->default_value(true) , "water vap source")
//      ("rc_src", po::value<bool>()->default_value(true) , "cloud water source (in blk_1/2m)")
//      ("rr_src", po::value<bool>()->default_value(true) , "rain water source (in blk_1/2m)")
//      ("nc_src", po::value<bool>()->default_value(true) , "cloud water concentration source (in blk_2m)")
//      ("nr_src", po::value<bool>()->default_value(true) , "rain water concentration source (in blk_2m)")
//      ("uv_src", po::value<bool>()->default_value(true) , "horizontal vel src")
//      ("w_src", po::value<bool>()->default_value(true) , "vertical vel src")
      ("piggy", po::value<bool>()->default_value(false) , "do piggybacking from a velocity field stored on a disk")
      ("sgs", po::value<bool>()->default_value(false) , "turn Eulerian SGS model on/off")
      ("sgs_delta", po::value<setup::real_t>()->default_value(-1) , "subgrid-scale turbulence model length scale [m]. If negative, sgs_delta = dz")
      ("help", "produce a help message (see also --micro X --help)")
      ("relax_th_rv", po::value<bool>()->default_value(false) , "relax per-level mean theta and rv to a desired (case-specific) profile")

      // aerosol distribution params
      // default values are realistic params, except n1_stp=n2_stp=-1
      // if n1_stp<0 and n2_stp<0, the case-default aerosols distribution is used,
      // concentration of this case-default distribution can be changed by setting the case_n_stp_multiplier
      ("mean_rd1", po::value<setup::real_t>()->default_value(1.0e-6) , "aerosol distirbution lognormal mode 1: mean radius [m] (lgrngn and blk_2m microphysics only)")
      ("sdev_rd1", po::value<setup::real_t>()->default_value(1.2) , "aerosol distirbution lognormal mode 1: geometric standard deviation (lgrngn and blk_2m microphysics only)")
      ("n1_stp", po::value<setup::real_t>()->default_value(-1.0) , "aerosol distirbution lognormal mode 1: concentration at STP [1/m^3] (lgrngn and blk_2m microphysics only). If n1_stp<0 and n2_stp<0, case-specific aerosol distribution is used.")
      ("kappa1", po::value<setup::real_t>()->default_value(0.61) , "aerosol distirbution lognormal mode 1: hygroscopicity parameter (lgrngn and blk_2m microphysics only)")
      ("mean_rd2", po::value<setup::real_t>()->default_value(1.0e-6) , "aerosol distirbution lognormal mode 2: mean radius [m] (lgrngn and blk_2m microphysics only)")
      ("sdev_rd2", po::value<setup::real_t>()->default_value(1.2) , "aerosol distirbution lognormal mode 2: geometric standard deviation (lgrngn and blk_2m microphysics only)")
      ("n2_stp", po::value<setup::real_t>()->default_value(-1.0) , "aerosol distirbution lognormal mode 2: concentration at STP [1/m^3] (lgrngn and blk_2m microphysics only). If n1_stp<0 and n2_stp<0, case-specific aerosol distribution is used.")
      ("kappa2", po::value<setup::real_t>()->default_value(0.61) , "aerosol distirbution lognormal mode 2: hygroscopicity parameter (lgrngn and blk_2m microphysics only)")
      ("case_n_stp_multiplier", po::value<setup::real_t>()->default_value(1.0) , "if case-specific aerosol distribution is used, multiply the case-default aerosols concentration by this value.")
//      ("relax_ccn", po::value<bool>()->default_value(false) , "add CCN if per-level mean of CCN concentration is lower than (case-specific) desired concentration")
    ;
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(opts_main).allow_unregistered().run(), vm); // ignores unknown

    // handling the "help" option
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
      user_params.outstart = vm["outstart"].as<int>();
      user_params.outwindow = vm["outwindow"].as<int>();
    }

    int
      nx = vm["nx"].as<int>(),
      ny = vm["ny"].as<int>(),
      nz = vm["nz"].as<int>();

    user_params.nt = vm["nt"].as<int>(),
    user_params.spinup = vm["spinup"].as<int>();

    // handling rng_seed
    user_params.rng_seed = vm["rng_seed"].as<int>();
    while(user_params.rng_seed == 0) //if = 0, get random seed
    {
      std::random_device rd;
      user_params.rng_seed = rd();
    }
    user_params.rng_seed_init = vm["rng_seed_init"].as<int>();
    if(user_params.rng_seed_init == 0) 
      user_params.rng_seed_init = user_params.rng_seed;

    user_params.dt = vm["dt"].as<setup::real_t>();
    user_params.X = vm["X"].as<setup::real_t>();
    user_params.Y = vm["Y"].as<setup::real_t>();
    user_params.Z = vm["Z"].as<setup::real_t>();

    user_params.window = vm["window"].as<bool>();

    // handling serial-advection-forcing flag
    if(vm["serial"].as<bool>()) setenv("OMP_NUM_THREADS", "1", 1);

    // handling sources flags
//    user_params.th_src = vm["th_src"].as<bool>();
//    user_params.rv_src = vm["rv_src"].as<bool>();
//    user_params.rc_src = vm["rc_src"].as<bool>();
//    user_params.rr_src = vm["rr_src"].as<bool>();
//    user_params.nc_src = vm["nc_src"].as<bool>();
//    user_params.nr_src = vm["nr_src"].as<bool>();
//    user_params.uv_src = vm["uv_src"].as<bool>();
//    user_params.w_src = vm["w_src"].as<bool>();

//    user_params.relax_ccn = vm["relax_ccn"].as<bool>();
    user_params.relax_th_rv = vm["relax_th_rv"].as<bool>();

    bool piggy = vm["piggy"].as<bool>();
    bool sgs = vm["sgs"].as<bool>();
    user_params.sgs_delta = vm["sgs_delta"].as<setup::real_t>();
    
// sanity check if desired options were compiled
#if defined(UWLCM_DISABLE_PIGGYBACKER)
    if(piggy)  throw std::runtime_error("UWLCM: Piggybacker option was disabled at compile time");
#endif
#if defined(UWLCM_DISABLE_DRIVER)
    if(!piggy)  throw std::runtime_error("UWLCM: Driver option was disabled at compile time");
#endif
#if defined(UWLCM_DISABLE_SGS)
    if(sgs)  throw std::runtime_error("UWLCM: SGS option was disabled at compile time");
#endif
#if defined(UWLCM_DISABLE_ILES)
    if(!sgs)  throw std::runtime_error("UWLCM: ILES option was disabled at compile time");
#endif

    if(piggy && sgs) throw std::runtime_error("UWLCM: SGS does not work in a piggybacker run");

    // set aerosol params to user_params data structure
    user_params.mean_rd1 = vm["mean_rd1"].as<setup::real_t>() * si::metres;
    user_params.sdev_rd1 = vm["sdev_rd1"].as<setup::real_t>();
    user_params.n1_stp = vm["n1_stp"].as<setup::real_t>() / si::cubic_metres;
    user_params.kappa1 = vm["kappa1"].as<setup::real_t>();
    user_params.mean_rd2 = vm["mean_rd2"].as<setup::real_t>() * si::metres;
    user_params.sdev_rd2 = vm["sdev_rd2"].as<setup::real_t>();
    user_params.n2_stp = vm["n2_stp"].as<setup::real_t>() / si::cubic_metres;
    user_params.kappa2 = vm["kappa2"].as<setup::real_t>();

    user_params.case_n_stp_multiplier = vm["case_n_stp_multiplier"].as<setup::real_t>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();

    // handling the "case" option
    user_params.model_case = vm["case"].as<std::string>();

    if(micro != "none" && user_params.model_case == "dry_thermal")
      throw std::runtime_error("UWLCM: The dry_thermal case needs micro set to 'none'");

    if(micro != "none" && user_params.model_case == "dry_pbl")
      throw std::runtime_error("UWLCM: The dry_pbl case needs micro set to 'none'");

    // run the simulation
    if (micro == "lgrngn" && ny == 0) // 2D super-droplet
#if !defined(UWLCM_DISABLE_2D_LGRNGN)
      run_hlpr<slvr_lgrngn, ct_params_2D_lgrngn>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 2D Lagrangian option was disabled at compile time");
#endif

    else if (micro == "lgrngn" && ny > 0) // 3D super-droplet
#if !defined(UWLCM_DISABLE_3D_LGRNGN)
      run_hlpr<slvr_lgrngn, ct_params_3D_lgrngn>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 3D Lagrangian option was disabled at compile time");
#endif

    else if (micro == "blk_1m" && ny == 0) // 2D one-moment
#if !defined(UWLCM_DISABLE_2D_BLK_1M)
      run_hlpr<slvr_blk_1m, ct_params_2D_blk_1m>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 2D Bulk 1-moment option was disabled at compile time");
#endif

    else if (micro == "blk_1m" && ny > 0) // 3D one-moment
#if !defined(UWLCM_DISABLE_3D_BLK_1M)
      run_hlpr<slvr_blk_1m, ct_params_3D_blk_1m>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 3D Bulk 1-moment option was disabled at compile time");
#endif

    else if (micro == "blk_2m" && ny == 0) // 2D two-moment
#if !defined(UWLCM_DISABLE_2D_BLK_2M)
      run_hlpr<slvr_blk_2m, ct_params_2D_blk_2m>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 2D Bulk 2-moment option was disabled at compile time");
#endif

    else if (micro == "blk_2m" && ny > 0) // 3D two-moment
#if !defined(UWLCM_DISABLE_3D_BLK_2M)      
      run_hlpr<slvr_blk_2m, ct_params_3D_blk_2m>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 3D Bulk 2-moment option was disabled at compile time");
#endif

    else if (micro == "none" && ny == 0) // 2D two-moment
#if !defined(UWLCM_DISABLE_2D_NONE)
      run_hlpr<slvr_dry, ct_params_2D_dry>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 2D none micro option was disabled at compile time");
#endif

    else if (micro == "none" && ny > 0) // 3D two-moment
#if !defined(UWLCM_DISABLE_3D_NONE)      
      run_hlpr<slvr_dry, ct_params_3D_dry>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);
#else
      throw std::runtime_error("UWLCM: 3D none micro option was disabled at compile time");
#endif

    // TODO: not only micro can be wrong
    else throw
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro"
      );
  }
}

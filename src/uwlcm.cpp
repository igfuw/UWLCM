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

#include "run_hlpr.hpp"

#include "detail/ct_params.hpp"

#include "opts/opts_common.hpp"

#include "solvers/slvr_lgrngn.hpp"
#include "solvers/slvr_blk_1m.hpp"
#include "forcings/calc_forces_blk_1m.hpp"
#include "forcings/calc_forces_common.hpp"

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
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("case", po::value<std::string>()->required(), "one of: dry_thermal, moist_thermal, dycoms, lasher_trapp")
      ("nx", po::value<int>()->default_value(76) , "grid cell count in horizontal")
      ("ny", po::value<int>()->default_value(0) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(76) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(3600) , "timestep count")
      ("rng_seed", po::value<int>()->default_value(0) , "rng seed, 0 for random")
      ("dt", po::value<setup::real_t>()->required() , "timestep length")
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
    while(user_params.rng_seed == 0) //if = 0, get random seed
    {
      std::random_device rd;
      user_params.rng_seed = rd();
    }
    std::cout << "rng seed: " << user_params.rng_seed << std::endl;

    //handling timestep length
    user_params.dt = vm["dt"].as<setup::real_t>();

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

    // sanity check if desired options were compiled
#if defined(UWLCM_DISABLE_PIGGYBACKER)
    if(piggy)  throw std::runtime_error("Piggybacker option was disabled at compile time");
#endif
#if defined(UWLCM_DISABLE_DRIVER)
    if(!piggy)  throw std::runtime_error("Driver option was disabled at compile time");
#endif
#if defined(UWLCM_DISABLE_SGS)
    if(sgs)  throw std::runtime_error("SGS option was disabled at compile time");
#endif
#if defined(UWLCM_DISABLE_ILES)
    if(!sgs)  throw std::runtime_error("ILES option was disabled at compile time");
#endif

    user_params.sgs_delta = vm["sgs_delta"].as<setup::real_t>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();

    // handling the "case" option
    user_params.model_case = vm["case"].as<std::string>();

    if(micro != "blk_1m" && user_params.model_case == "dry_thermal")
      throw std::runtime_error("The DryThermal case needs micro set to blk_1m (although there's no microphysics done at all");

    // run the simulation
    if (micro == "lgrngn" && ny == 0) // 2D super-droplet
#if !defined(UWLCM_DISABLE_2D_LGRNGN)
      run_hlpr<slvr_lgrngn, ct_params_2D_lgrngn>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);
#else
      throw std::runtime_error("2D Lagrangian option was disabled at compile time");
#endif

    else if (micro == "lgrngn" && ny > 0) // 3D super-droplet
#if !defined(UWLCM_DISABLE_3D_LGRNGN)
      run_hlpr<slvr_lgrngn, ct_params_3D_lgrngn>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);
#else
      throw std::runtime_error("3D Lagrangian option was disabled at compile time");
#endif

    else if (micro == "blk_1m" && ny == 0) // 2D one-moment
#if !defined(UWLCM_DISABLE_2D_BLK_1M)
      run_hlpr<slvr_blk_1m, ct_params_2D_blk_1m>(piggy, sgs, user_params.model_case, {nx, nz}, user_params);
#else
      throw std::runtime_error("2D Bulk 1-moment option was disabled at compile time");
#endif

    else if (micro == "blk_1m" && ny > 0) // 3D one-moment
#if !defined(UWLCM_DISABLE_3D_BLK_1M)
      run_hlpr<slvr_blk_1m, ct_params_3D_blk_1m>(piggy, sgs, user_params.model_case, {nx, ny, nz}, user_params);
#else
      throw std::runtime_error("3D Bulk 1-moment option was disabled at compile time");
#endif

  // TODO: not only micro can be wrong
    else throw
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro"
      );
  }
}

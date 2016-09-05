/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/bcond/cyclic_3d.hpp>
#include <libmpdata++/bcond/open_3d.hpp>
#include <libmpdata++/concurr/boost_thread.hpp> // not to conflict with OpenMP used via Thrust in libcloudph++
#include <libmpdata++/concurr/serial.hpp> // not to conflict with OpenMP used via Thrust in libcloudph++
#include "setup.hpp"
#include "opts_lgrngn.hpp"
#include "panic.hpp"
#include <map>

// simulation parameters container
// TODO: write them to rt_params directly in main()
struct user_params_t
{
  int nt, outfreq, spinup, rng_seed;
  setup::real_t dt, z_rlx_sclr;
  std::string outdir;
  bool serial, relax_th_rv;
};

// copy external profiles into rt_parameters
// TODO: more elegant way
template<class params_t>
void copy_profiles(setup::arr_1D_t &th_e, setup::arr_1D_t &rv_e, setup::arr_1D_t &th_ref, setup::arr_1D_t &rhod, setup::arr_1D_t &w_LS, setup::arr_1D_t &hgt_v, setup::arr_1D_t &hgt_s, params_t &p)
{
  p.hgt_fctr_sclr = new setup::arr_1D_t(hgt_s.dataFirst(), hgt_s.shape(), blitz::neverDeleteData);
  p.hgt_fctr_vctr = new setup::arr_1D_t(hgt_v.dataFirst(), hgt_v.shape(), blitz::neverDeleteData);
  p.th_e = new setup::arr_1D_t(th_e.dataFirst(), th_e.shape(), blitz::neverDeleteData);
  p.rv_e = new setup::arr_1D_t(rv_e.dataFirst(), rv_e.shape(), blitz::neverDeleteData);
  p.th_ref = new setup::arr_1D_t(th_ref.dataFirst(), th_ref.shape(), blitz::neverDeleteData);
  p.rhod = new setup::arr_1D_t(rhod.dataFirst(), rhod.shape(), blitz::neverDeleteData);
  p.w_LS = new setup::arr_1D_t(w_LS.dataFirst(), w_LS.shape(), blitz::neverDeleteData);
}

// 2D model run logic - the same for any microphysics
template <class solver_t>
void run(int nx, int nz, const user_params_t &user_params)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::rt_params_t p;

  // output and simulation parameters
  p.grid_size = {nx, nz};

  setup::setopts(p, nx, nz, user_params);
  setopts_micro<solver_t>(p, user_params);

  // reference profiles shared among threads
  setup::arr_1D_t th_e(nz), rv_e(nz), th_ref(nz), rhod(nz+1), w_LS(nz), hgt_fctr_vctr(nz), hgt_fctr_sclr(nz); 
  // rhod needs to be bigger, cause it divides vertical courant number, TODO: should have a halo both up and down, not only up like now; then it should be interpolated in courant calculation

  // assign their values
  setup::env_prof(th_e, rv_e, th_ref, rhod, w_LS, hgt_fctr_vctr, hgt_fctr_sclr, nz, user_params);
  // pass them to rt_params
  copy_profiles(th_e, rv_e, th_ref, rhod, w_LS, hgt_fctr_vctr, hgt_fctr_sclr, p);

  // solver instantiation
  std::unique_ptr<
    concurr::any<
      typename solver_t::real_t, 
      solver_t::n_dims
    >
  > slv;
  if (user_params.serial)
  {
    using concurr_t = concurr::serial<
      solver_t, 
      bcond::cyclic, bcond::cyclic,
      bcond::rigid,  bcond::rigid 
    >;
    slv.reset(new concurr_t(p));

    // initial condition
    setup::intcond<2>(*static_cast<concurr_t*>(slv.get()), rhod, user_params.rng_seed);
  }
  else
  {
    using concurr_t = concurr::boost_thread<
      solver_t, 
      bcond::cyclic, bcond::cyclic,
      bcond::rigid,  bcond::rigid 
    >;
    slv.reset(new concurr_t(p));

    // initial condition
    setup::intcond<2>(*static_cast<concurr_t*>(slv.get()), rhod, user_params.rng_seed);
  }

  // setup panic pointer and the signal handler
  panic = slv->panic_ptr();
  set_sigaction();
 
  // timestepping
  slv->advance(user_params.nt);
}

// 3D model run logic - the same for any microphysics; TODO: still a lot of common code with 2D run
template <class solver_t>
void run(int nx, int ny, int nz, const user_params_t &user_params)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::rt_params_t p;

  // output and simulation parameters
  p.grid_size = {nx, ny, nz};

  setup::setopts(p, nx, ny, nz, user_params);
  setopts_micro<solver_t>(p, user_params);

  // reference profiles shared among threads
  setup::arr_1D_t th_e(nz), rv_e(nz), th_ref(nz), rhod(nz+1), w_LS(nz), hgt_fctr_vctr(nz), hgt_fctr_sclr(nz);
  // assign their values
  setup::env_prof(th_e, rv_e, th_ref, rhod, w_LS, hgt_fctr_vctr, hgt_fctr_sclr, nz, user_params);
  // pass them to rt_params
  copy_profiles(th_e, rv_e, th_ref, rhod, w_LS, hgt_fctr_vctr, hgt_fctr_sclr, p);

  // solver instantiation
  std::unique_ptr<
    concurr::any<
      typename solver_t::real_t, 
      solver_t::n_dims
    >
  > slv;
  if (user_params.serial)
  {
    using concurr_t = concurr::serial<
      solver_t, 
      bcond::cyclic, bcond::cyclic,
      bcond::cyclic, bcond::cyclic,
      bcond::rigid,  bcond::rigid 
    >;
    slv.reset(new concurr_t(p));

    // initial condition
    setup::intcond<3>(*static_cast<concurr_t*>(slv.get()), rhod, user_params.rng_seed);
  }
  else
  {
    using concurr_t = concurr::boost_thread<
      solver_t, 
      bcond::cyclic, bcond::cyclic,
      bcond::cyclic, bcond::cyclic,
      bcond::rigid,  bcond::rigid 
    >;
    slv.reset(new concurr_t(p));

    // initial condition
    setup::intcond<3>(*static_cast<concurr_t*>(slv.get()), rhod, user_params.rng_seed);
  }

  // setup panic pointer and the signal handler
  panic = slv->panic_ptr();
  set_sigaction();
 
  // timestepping
  slv->advance(user_params.nt);
}

// libmpdata++'s compile-time parameters
struct ct_params_common : ct_params_default_t
{
  using real_t = setup::real_t;
  enum { opts = opts::nug | opts::iga | opts::fct };  // TODO: reenable nug once it works in 3D
  enum { rhs_scheme = solvers::trapez }; 
  enum { prs_scheme = solvers::cr };
  enum { vip_vab = solvers::expl };
};


// all starts here with handling general options 
int main(int argc, char** argv)
{
  // making argc and argv global
  ac = argc;
  av = argv;

  {
    // note: all options should have default values here to make "--micro=? --help" work
    opts_main.add_options()
      ("micro", po::value<std::string>()->required(), "one of: blk_1m, blk_2m, lgrngn")
      ("nx", po::value<int>()->default_value(76) , "grid cell count in horizontal")
      ("ny", po::value<int>()->default_value(0) , "grid cell count in horizontal")
      ("nz", po::value<int>()->default_value(76) , "grid cell count in vertical")
      ("nt", po::value<int>()->default_value(3600) , "timestep count")
      ("rng_seed", po::value<int>()->default_value(-1) , "rng seed, negative for random")
      ("dt", po::value<setup::real_t>()->required() , "timestep length")
      ("z_rlx_sclr", po::value<setup::real_t>()->default_value(10) , "scalars surface flux charasteristic heihjt")
      ("outdir", po::value<std::string>(), "output file name (netCDF-compatible HDF5)")
      ("outfreq", po::value<int>(), "output rate (timestep interval)")
      ("spinup", po::value<int>()->default_value(2400) , "number of initial timesteps during which rain formation is to be turned off")
      ("adv_serial", po::value<bool>()->default_value(false), "force advection to be computed on single thread")
      ("relax_th_rv", po::value<bool>()->default_value(false) , "relaxation of th and rv")
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
    user_params.serial = vm["adv_serial"].as<bool>();

    // handling relaxation flag
    user_params.relax_th_rv = vm["relax_th_rv"].as<bool>();

    // handling the "micro" option
    std::string micro = vm["micro"].as<std::string>();

    if (micro == "lgrngn" && ny == 0) // 2D super-droplet
    {
      struct ct_params_t : ct_params_common
      {
        enum { n_dims = 2 };
    	enum { n_eqns = 4 };
        struct ix { enum {
          u, w, th, rv, 
          vip_i=u, vip_j=w, vip_den=-1
        }; };
      };
      run<slvr_lgrngn<ct_params_t>>(nx, nz, user_params);
    }
    else if (micro == "lgrngn" && ny > 0) // 3D super-droplet
    {
      struct ct_params_t : ct_params_common
      {
        enum { n_dims = 3 };
    	enum { n_eqns = 5 };
        struct ix { enum {
          u, v, w, th, rv, 
          vip_i=u, vip_j=v, vip_k=w, vip_den=-1
        }; };
      };
      run<slvr_lgrngn<ct_params_t>>(nx, ny, nz, user_params);
    }
    else throw
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro" 
      );
  }
}

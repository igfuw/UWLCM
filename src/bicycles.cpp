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

// model run logic - the same for any microphysics
template <class solver_t>
void run(int nx, int nz, int nt, setup::real_t dt, const std::string &outdir, const int &outfreq, int spinup, bool serial, bool relax_th_rv, bool gccn, bool onishi, bool pristine, setup::real_t eps, setup::real_t ReL, setup::real_t z_rlx_sclr, int rng_seed)
{
  // instantiation of structure containing simulation parameters
  typename solver_t::rt_params_t p;

  // output and simulation parameters
  p.grid_size = {nx, nz};
  p.outdir = outdir;
  p.outfreq = outfreq;
  p.spinup = spinup;
  p.relax_th_rv = relax_th_rv;
  p.prs_tol=1e-6;
  p.dt = dt;
  setopts_micro<solver_t>(p, nx, nz, nt, gccn, onishi, pristine, eps, ReL);
  //std::cout << "params.rhod po setopts micro " << p.rhod << " "  << *p.rhod << std::endl;
  setup::setopts(p, nx, nz);

  // --------
  // reference profiles init, they will be passed to solvers through rt_params_t
  blitz::secondIndex k;
  // env profiles of th and rv (for buoyancy)
  setup::arr_1D_t th_e(nx, nz), rv_e(nx, nz), th_ref(nx, nz), rhod(nx, nz);
  setup::env_prof(th_e, rv_e, th_ref, rhod, nz);
  p.th_e = new typename setup::arr_1D_t(th_e.dataFirst(), th_e.shape(), blitz::neverDeleteData);
  p.rv_e = new typename setup::arr_1D_t(rv_e.dataFirst(), rv_e.shape(), blitz::neverDeleteData);
  p.th_ref = new typename setup::arr_1D_t(th_ref.dataFirst(), th_ref.shape(), blitz::neverDeleteData);
  p.rhod = new typename setup::arr_1D_t(rhod.dataFirst(), rhod.shape(), blitz::neverDeleteData);
  // subsidence rate
  typename setup::arr_1D_t w_LS(nx, nz);
  w_LS = setup::w_LS_fctr()(k * p.dz);
  p.w_LS = new typename setup::arr_1D_t(w_LS.dataFirst(), w_LS.shape(), blitz::neverDeleteData);
  // surface sources relaxation factors
  // for vectors
  typename setup::arr_1D_t hgt_fctr_vctr(nx, nz);
  setup::real_t z_0 = setup::z_rlx_vctr / si::metres;
  hgt_fctr_vctr = exp(- (k-0.5) * p.dz / z_0); // z=0 at k=1/2
  hgt_fctr_vctr(blitz::Range::all(),0) = 1;
  p.hgt_fctr_vctr = new typename setup::arr_1D_t(hgt_fctr_vctr.dataFirst(), hgt_fctr_vctr.shape(), blitz::neverDeleteData);
  // for scalars
  typename setup::arr_1D_t hgt_fctr_sclr(nx, nz);
  z_0 = z_rlx_sclr;
  hgt_fctr_sclr = exp(- (k-0.5) * p.dz / z_0);
  hgt_fctr_sclr(blitz::Range::all(),0) = 1;
  p.hgt_fctr_sclr = new typename setup::arr_1D_t(hgt_fctr_sclr.dataFirst(), hgt_fctr_sclr.shape(), blitz::neverDeleteData);

  // --------
  // solver instantiation
  std::unique_ptr<
    concurr::any<
      typename solver_t::real_t, 
      solver_t::n_dims
    >
  > slv;
  if (serial)
  {
    using concurr_t = concurr::serial<
      solver_t, 
      bcond::cyclic, bcond::cyclic,
      bcond::rigid,  bcond::rigid 
    >;
    slv.reset(new concurr_t(p));

    // initial condition
    setup::intcond(*static_cast<concurr_t*>(slv.get()), rhod, rng_seed);
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
    setup::intcond(*static_cast<concurr_t*>(slv.get()), rhod, rng_seed);
  }


  // setup panic pointer and the signal handler
  panic = slv->panic_ptr();
  set_sigaction();
 
  // timestepping
  slv->advance(nt);
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
      ("relax_th_rv", po::value<bool>()->default_value(true) , "relaxation of th and rv")
      ("gccn", po::value<bool>()->default_value(false) , "add GCCNs")
      ("onishi", po::value<bool>()->default_value(false) , "use the turbulent onishi kernel")
      ("pristine", po::value<bool>()->default_value(false) , "pristine conditions")
      ("eps", po::value<setup::real_t>()->default_value(0.01) , "turb dissip rate (for onishi kernel) [m^2/s^3]")
      ("ReL", po::value<setup::real_t>()->default_value(5000) , "taylor-microscale reynolds number (onishi kernel)")
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
    
    // handling outdir && outfreq
    std::string outdir; 
    int outfreq;
    if (!vm.count("help"))
    {
      if (!vm.count("outdir")) throw po::required_option("outdir");
      if (!vm.count("outfreq")) throw po::required_option("outfreq");
      outdir = vm["outdir"].as<std::string>();
      outfreq = vm["outfreq"].as<int>();
    }

    // handling nx, nz, nt options
    int 
      nx = vm["nx"].as<int>(),
      ny = vm["ny"].as<int>(),
      nz = vm["nz"].as<int>(),
      nt = vm["nt"].as<int>(),
      spinup = vm["spinup"].as<int>();
 
    // handling rng_seed
    int rng_seed = vm["rng_seed"].as<int>();
   
    //handling timestep length
    setup::real_t dt = vm["dt"].as<setup::real_t>();

    //handling z_rlx_sclr
    setup::real_t z_rlx_sclr = vm["z_rlx_sclr"].as<setup::real_t>();

    // handling serial-advection-forcing flag
    bool adv_serial = vm["adv_serial"].as<bool>();

    // handling relaxation flag
    bool relax_th_rv = vm["relax_th_rv"].as<bool>();

    // handling GCCN flag
    bool gccn = vm["gccn"].as<bool>();

    // handling pristine flag
    bool pristine = vm["pristine"].as<bool>();

    // handling onishi flag
    bool onishi = vm["onishi"].as<bool>();
   
    //handling epsilon
    setup::real_t eps = vm["eps"].as<setup::real_t>();
   
    //handling Re lambda
    setup::real_t ReL = vm["ReL"].as<setup::real_t>();

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
      run<slvr_lgrngn<ct_params_t>>(nx, nz, nt, dt, outdir, outfreq, spinup, adv_serial, relax_th_rv, gccn, onishi, pristine, eps, ReL, z_rlx_sclr, rng_seed);
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
      run<slvr_lgrngn<ct_params_t>>(nx, nz, nt, dt, outdir, outfreq, spinup, adv_serial, relax_th_rv, gccn, onishi, pristine, eps, ReL, z_rlx_sclr, rng_seed);
    }
    else throw
      po::validation_error(
        po::validation_error::invalid_option_value, micro, "micro" 
      );
  }
}

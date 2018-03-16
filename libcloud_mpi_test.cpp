#include <libcloudph++/lgrngn/factory.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#include <stdio.h>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>
#include <iostream>
#include "mpi.h"


using namespace std;
using namespace libcloudphxx::lgrngn;
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  //aerosol bimodal lognormal dist. 
  const quantity<si::length, double>
    mean_rd1 = double(30e-6) * si::metres,
    mean_rd2 = double(40e-6) * si::metres;
  const quantity<si::dimensionless, double>
    sdev_rd1 = double(1.4),
    sdev_rd2 = double(1.6);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, double>
    n1_stp = double(60e6) / si::cubic_metres,
    n2_stp = double(40e6) / si::cubic_metres;


/*
  // lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii : public libcloudphxx::common::unary_function<T>
  {
    const quantity<si::length, real_t> mean_rd1, mean_rd2;
    const quantity<si::dimensionless, real_t> sdev_rd1, sdev_rd2;
    const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_stp, n2_stp;

    log_dry_radii(
      quantity<si::length, real_t> mean_rd1,
      quantity<si::length, real_t> mean_rd2,
      quantity<si::dimensionless, real_t> sdev_rd1,
      quantity<si::dimensionless, real_t> sdev_rd2,
      quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_stp,
      quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n2_stp
    ):
    mean_rd1(mean_rd1),
    mean_rd2(mean_rd2),
    sdev_rd1(sdev_rd1),
    sdev_rd2(sdev_rd2),
    n1_stp(n1_stp),
    n2_stp(n2_stp) {}


    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }
  };
*/

// lognormal aerosol distribution
template <typename T>
struct log_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    return T(( 
        lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, double>(lnrd)) +
        lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, double>(lnrd)) 
      ) * si::cubic_metres
    );  
  }   
};  


void two_step(particles_proto_t<double> *prtcls, 
             arrinfo_t<double> th,
             arrinfo_t<double> rhod,
             arrinfo_t<double> rv,
             opts_t<double> opts)
{
    prtcls->step_sync(opts,th,rv,rhod);
    prtcls->step_async(opts);
}

void test(backend_t backend, int ndims)
{
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  opts_init_t<double> opts_init;
  opts_init.dt=3.;
  opts_init.sstp_coal = 1; 
  opts_init.kernel = kernel_t::geometric;
  opts_init.terminal_velocity = vt_t::beard76;
  opts_init.dx = 1;
  opts_init.nx = 2*(rank+1); 
  opts_init.x1 = 2*(rank+1);
  opts_init.sd_conc = 64;
  opts_init.n_sd_max = 10*opts_init.sd_conc;
  opts_init.rng_seed = 4444;// + rank;
  if(ndims>1)
  {
    opts_init.dz = 1; 
    opts_init.nz = 1; 
    opts_init.z1 = 1; 
  }
  opts_init.dev_id = rank; 
  std::cout << opts_init.dev_id << std::endl;
//  opts_init.sd_const_multi = 1;

/*
  boost::assign::ptr_map_insert<
    log_dry_radii<double> // value type
  >(  
    opts_init.dry_distros // map
  )(  
    0.001 // key
  ); 
*/

  opts_init.dry_distros.emplace( 
    0.001, // kappa
    std::make_shared<log_dry_radii<double>>() // distribution
  );

  particles_proto_t<double> *prtcls;

  printf("nx = %d\n", opts_init.nx);
  printf("ny = %d\n", opts_init.ny);
  printf("nz = %d\n", opts_init.nz);
  prtcls = factory<double>(
    backend,
  //  (backend_t)multi_CUDA, 
  //  (backend_t)CUDA, 
  //  (backend_t)serial, 
    opts_init
  );
  printf("po factory\n");
  double pth[] = {300., 300., 300., 300.};
  double prhod[] = {1.225, 1.225, 1.225, 1.225};
  double prv[] = {.01, 0.01, 0.01, 0.01};
  double pCx[] = {-1, -1, -1, -1, -1, -1, -1};
  double pCz[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  //long int strides[] = {sizeof(double)};
  long int strides[] = {1, 1};
  long int xstrides[] = {1, 1};
  long int ystrides[] = {1, 1};

  arrinfo_t<double> th(pth, strides);
  arrinfo_t<double> rhod(prhod, strides);
  arrinfo_t<double> rv(prv, strides);
  arrinfo_t<double> Cx(pCx, xstrides);
  arrinfo_t<double> Cz(pCz, ystrides);

  if(ndims==1)
    prtcls->init(th,rv,rhod, Cx);
  else if(ndims==2)
    prtcls->init(th,rv, rhod, Cx, arrinfo_t<double>(), Cz);

  opts_t<double> opts;
  opts.adve = 0;
  opts.sedi = 0;
  opts.cond = 0;
  opts.coal = 1;
//  opts.chem = 0;


  prtcls->diag_all();
  prtcls->diag_sd_conc();
  double *out = prtcls->outbuf();
  printf("---sd_conc init---\n");
  printf("%d: %lf %lf %lf %lf\n",rank, out[0], out[1], out[2], out[3]);
  MPI_Barrier(MPI_COMM_WORLD);
  

  for(int i=0;i<70;++i)
  {
//    if(rank==0)
      two_step(prtcls,th,rhod,rv,opts);
//   //   MPI_Barrier(MPI_COMM_WORLD);
  //  if(rank==1)
   //   two_step(prtcls,th,rhod,rv,opts);
//   // MPI_Barrier(MPI_COMM_WORLD);
  }
  prtcls->diag_all();
  prtcls->diag_sd_conc();
  out = prtcls->outbuf();
  printf("---sd_conc po coal---\n");
  printf("%d: %lf %lf %lf %lf\n",rank, out[0], out[1], out[2], out[3]);
  opts.coal = 0;
  opts.adve = 1;
  two_step(prtcls,th,rhod,rv,opts);
  prtcls->diag_all();
  prtcls->diag_sd_conc();
  out = prtcls->outbuf();

  MPI_Barrier(MPI_COMM_WORLD);
  printf("---sd_conc po adve---\n");
  printf("%d: %lf %lf %lf %lf\n",rank, out[0], out[1], out[2], out[3]);
}

int main(int argc, char *argv[]){
  int ndims;
  sscanf(argv[1], "%d", &ndims);
  printf("ndims %d\n", ndims);

  int provided_thread_lvl;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided_thread_lvl);
  printf("provided thread lvl: %d\n", provided_thread_lvl);
 
  test(backend_t(multi_CUDA), ndims);
  MPI_Finalize();
}

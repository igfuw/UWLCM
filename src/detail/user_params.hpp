#pragma once

#include <string>
#include "setup.hpp"

// simulation parameters container
// TODO: write them to rt_params directly in main()
struct user_params_t
{
  int nt, outfreq, spinup, rng_seed, rng_seed_init;
  setup::real_t dt;
  std::string outdir, model_case;
  bool th_src, rv_src, rc_src, rr_src, nc_src, nr_src, uv_src, w_src, rng_seed_init_switch;
  setup::real_t sgs_delta;
  quantity<si::length, setup::real_t> mean_rd1, mean_rd2;		
  quantity<si::dimensionless, setup::real_t> sdev_rd1, sdev_rd2;		
  quantity<power_typeof_helper<si::length, static_rational<-3>>::type, setup::real_t> n1_stp, n2_stp;		
  quantity<si::dimensionless, setup::real_t> kappa1, kappa2;
};

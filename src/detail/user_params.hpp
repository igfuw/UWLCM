#pragma once

#include <string>
#include "setup.hpp"

// user-defined simulation parameters container
// note: description and default values are in uwlcm.cpp, all parameters have to be handled there
struct user_params_t
{
  int nt, outfreq, outstart, outwindow, spinup, rng_seed, rng_seed_init;
  setup::real_t X, Y, Z, dt;
  std::string outdir, model_case;
  setup::real_t sgs_delta;
  quantity<si::length, setup::real_t> mean_rd1, mean_rd2;		
  quantity<si::dimensionless, setup::real_t> sdev_rd1, sdev_rd2;		
  quantity<power_typeof_helper<si::length, static_rational<-3>>::type, setup::real_t> n1_stp, n2_stp;		
  quantity<si::dimensionless, setup::real_t> kappa1, kappa2;
  quantity<si::dimensionless, setup::real_t> case_n_stp_multiplier;

  bool relax_th_rv,
       window,
       relax_ccn = false; // relevant only for lgrngn micro, hence needs a default value as otherwise it might be undefined in blk_1m/blk_2m
};

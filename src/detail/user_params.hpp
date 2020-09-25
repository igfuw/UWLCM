#pragma once

#include <string>
#include "setup.hpp"

// simulation parameters container
// TODO: write them to rt_params directly in main()
struct user_params_t
{
  int nt, outfreq, spinup, rng_seed;
  setup::real_t dt;
  std::string outdir, model_case;
  bool th_src, rv_src, rc_src, rr_src, nc_src, nr_src, uv_src, w_src;
  setup::real_t sgs_delta;
};

/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "opts_common.hpp"
#include "slvr_blk_1m.hpp"

// simulation and output parameters for micro=blk_1m
template <class solver_t, class user_params_t, class case_ptr_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params, 
  const user_params_t &user_params,
  const case_ptr_t &case_ptr,
  typename std::enable_if<std::is_same<
    decltype(solver_t::rt_params_t::cloudph_opts),
    libcloudphxx::blk_1m::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  po::options_description opts("Single-moment bulk microphysics options"); 
  opts.add_options()
    ("cond", po::value<bool>()->default_value(rt_params.cloudph_opts.cond) , "cloud water condensation (1=on, 0=off)")
    ("cevp", po::value<bool>()->default_value(rt_params.cloudph_opts.cevp) , "cloud water evaporation (1=on, 0=off)")
    ("revp", po::value<bool>()->default_value(rt_params.cloudph_opts.revp) , "rain water evaporation (1=on, 0=off)")
    ("conv", po::value<bool>()->default_value(rt_params.cloudph_opts.conv) , "autoconversion of cloud water into rain (1=on, 0=off)")
    ("accr", po::value<bool>()->default_value(rt_params.cloudph_opts.accr) , "cloud water collection by rain (1=on, 0=off)")
    ("sedi", po::value<bool>()->default_value(rt_params.cloudph_opts.sedi) , "rain water sedimentation (1=on, 0=off)")
    ("r_c0",   po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.r_c0)   , "blk 1m autoconversion treshold")
    ("k_acnv", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.k_acnv) , "blk 1m autoconversion rate")
    ("r_eps",  po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.r_eps)  , "blk 1m absolute tolerance")
//TODO: autoconv_threshold, epsilon
  ;
  po::variables_map vm;
  handle_opts(opts, vm);

  // Kessler scheme options
  rt_params.cloudph_opts.cond = vm["cond"].as<bool>();
  rt_params.cloudph_opts.cevp = vm["cevp"].as<bool>();
  rt_params.cloudph_opts.revp = vm["revp"].as<bool>();
  rt_params.cloudph_opts.conv = vm["conv"].as<bool>();
  rt_params.cloudph_opts.accr = vm["accr"].as<bool>();
  rt_params.cloudph_opts.sedi = vm["sedi"].as<bool>();
  rt_params.cloudph_opts.r_eps = 1e-6;
  rt_params.cloudph_opts.r_c0   = vm["r_c0"].as<typename solver_t::real_t>();
  rt_params.cloudph_opts.k_acnv = vm["k_acnv"].as<typename solver_t::real_t>();
  rt_params.cloudph_opts.r_eps  = vm["r_eps"].as<typename solver_t::real_t>();

  // output variables
  rt_params.outvars.insert({solver_t::ix::rc, {"rc", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::rr, {"rr", "[kg kg-1]"}});
}

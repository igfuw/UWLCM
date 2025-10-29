/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libcloudph++/blk_1m/options.hpp>
#include "opts_common.hpp"

// simulation and output parameters for micro=blk_1m
template <class solver_t, class user_params_t, class case_ptr_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params,
  const user_params_t &user_params,
  const case_ptr_t &case_ptr,
  uwlcm_blk_1m_ice_family_tag
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
    ("homA1", po::value<bool>()->default_value(rt_params.cloudph_opts.homA1) , "homogeneous nucleation of ice A from water vapor (1=on, 0=off)")
    ("homA2", po::value<bool>()->default_value(rt_params.cloudph_opts.homA2) , "homogeneous nucleation of ice A from cloud droplets (1=on, 0=off)")
    ("hetA", po::value<bool>()->default_value(rt_params.cloudph_opts.hetA) , "heterogeneous nucleation of ice A (1=on, 0=off)")
    ("hetB", po::value<bool>()->default_value(rt_params.cloudph_opts.hetB) , "heterogeneous nucleation of ice B (1=on, 0=off)")
    ("depA", po::value<bool>()->default_value(rt_params.cloudph_opts.depA) , "depositional growth of ice A (1=on, 0=off)")
    ("depB", po::value<bool>()->default_value(rt_params.cloudph_opts.depB) , "depositional growth of ice B (1=on, 0=off)")
    ("rimA", po::value<bool>()->default_value(rt_params.cloudph_opts.rimA) , "growth of ice A by riming (1=on, 0=off)")
    ("rimB", po::value<bool>()->default_value(rt_params.cloudph_opts.rimB) , "growth of ice B by riming (1=on, 0=off)")
    ("melA", po::value<bool>()->default_value(rt_params.cloudph_opts.melA) , "melting of ice A (1=on, 0=off)")
    ("melB", po::value<bool>()->default_value(rt_params.cloudph_opts.melB) , "melting of ice B (1=on, 0=off)")
    ("r_c0", po::value<setup::real_t>()->default_value(rt_params.cloudph_opts.r_c0) , "rain autoconversion threshold [kg/kg]")
    ("k_acnv", po::value<setup::real_t>()->default_value(rt_params.cloudph_opts.k_acnv) , "rain autoconversion parameter")
    ("r_eps", po::value<setup::real_t>()->default_value(1e-6) , "absolute tolerance")
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
  rt_params.cloudph_opts.homA1 = vm["homA1"].as<bool>();
  rt_params.cloudph_opts.homA2 = vm["homA2"].as<bool>();
  rt_params.cloudph_opts.hetA = vm["hetA"].as<bool>();
  rt_params.cloudph_opts.hetB = vm["hetB"].as<bool>();
  rt_params.cloudph_opts.depA = vm["depA"].as<bool>();
  rt_params.cloudph_opts.depB = vm["depB"].as<bool>();
  rt_params.cloudph_opts.rimA = vm["rimA"].as<bool>();
  rt_params.cloudph_opts.rimB = vm["rimB"].as<bool>();
  rt_params.cloudph_opts.melA = vm["melA"].as<bool>();
  rt_params.cloudph_opts.melB = vm["melB"].as<bool>();
  rt_params.cloudph_opts.r_c0 = vm["r_c0"].as<setup::real_t>();
  rt_params.cloudph_opts.k_acnv = vm["k_acnv"].as<setup::real_t>();
  rt_params.cloudph_opts.r_eps = vm["r_eps"].as<setup::real_t>();

  rt_params.cloudph_opts.th_dry = false;
  rt_params.cloudph_opts.const_p = true;
  rt_params.cloudph_opts.adj_nwtrph = true;

  // output variables
  rt_params.outvars.insert({solver_t::ix::rc, {"rc", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::rr, {"rr", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::ria, {"ria", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::rib, {"rib", "[kg kg-1]"}});

  if(rt_params.aerosol_independent_of_rhod)
    std::cerr << "UWLCM warning: aerosol_independent_of_rhod has no effect on blk_1m microphysics" << std::endl;
  if(rt_params.aerosol_conc_factor.size()!=0)
    std::cerr << "UWLCM warning: aerosol_conc_factor has no effect on blk_1m microphysics" << std::endl;
}

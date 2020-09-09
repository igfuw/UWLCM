/**
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

 #pragma once

 #include <libcloudph++/blk_2m/options.hpp>
 #include "opts_common.hpp"
 
// simulation and output parameters for micro=blk_2m
template <class solver_t, class user_params_t, class case_ptr_t>
void setopts_micro(
  typename solver_t::rt_params_t &rt_params,
  const user_params_t &user_params,
  const case_ptr_t &case_ptr,
  typename std::enable_if<std::is_same<
    decltype(solver_t::rt_params_t::cloudph_opts),
    libcloudphxx::blk_2m::opts_t<typename solver_t::real_t>
  >::value>::type* = 0
)
{
  po::options_description opts("Double-moment bulk microphysics options");
  opts.add_options()
    ("acti", po::value<bool>()->default_value(rt_params.cloudph_opts.acti) , "TODO (on/off)")
    ("cond", po::value<bool>()->default_value(rt_params.cloudph_opts.cond) , "TODO (on/off)")
    ("accr", po::value<bool>()->default_value(rt_params.cloudph_opts.accr) , "TODO (on/off)")
    ("acnv", po::value<bool>()->default_value(rt_params.cloudph_opts.acnv) , "TODO (on/off)")
    ("sedi", po::value<bool>()->default_value(rt_params.cloudph_opts.sedi) , "TODO (on/off)")
    ("acnv_A", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_A), "parameter in autoconversion rate formulae")
    ("acnv_b", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_b), "parameter in autoconversion rate formulae")
    ("acnv_c", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_c), "parameter in autoconversion rate formulae")
    ("blk2m_mean_rd", po::value<typename solver_t::real_t>()->default_value(0.02e-6), "mean aerosol dry radius [m]")
    ("blk2m_sdev_rd", po::value<typename solver_t::real_t>()->default_value(1.4),     "aerosol standard deviation")
    ("blk2m_N_stp",   po::value<typename solver_t::real_t>()->default_value(60e6),    "aerosol concentration [1/m3]")
    ("blk2m_chem_b",  po::value<typename solver_t::real_t>()->default_value(.55),     "kappa - chemical composition parameter") 
 ;
  po::variables_map vm;
  handle_opts(opts, vm);

  // Morrison and Grabowski 2007 scheme options
  rt_params.cloudph_opts.acti = vm["acti"].as<bool>();
  rt_params.cloudph_opts.cond = vm["cond"].as<bool>();
  rt_params.cloudph_opts.accr = vm["accr"].as<bool>();
  rt_params.cloudph_opts.acnv = vm["acnv"].as<bool>();
  rt_params.cloudph_opts.sedi = vm["sedi"].as<bool>();
  rt_params.cloudph_opts.acnv_A = vm["acnv_A"].as<typename solver_t::real_t>();
  rt_params.cloudph_opts.acnv_b = vm["acnv_b"].as<typename solver_t::real_t>();
  rt_params.cloudph_opts.acnv_c = vm["acnv_c"].as<typename solver_t::real_t>();

  rt_params.cloudph_opts.dry_distros.push_back({
    .mean_rd = vm["blk2m_mean_rd"].as<typename solver_t::real_t>(),
    .sdev_rd = vm["blk2m_sdev_rd"].as<typename solver_t::real_t>(),
    .N_stp   = vm["blk2m_N_stp"].as<typename solver_t::real_t>(),
    .chem_b  = vm["blk2m_chem_b"].as<typename solver_t::real_t>()
  });

  // output variables
  rt_params.outvars.insert({solver_t::ix::rc, {"rc", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::rr, {"rr", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::nc, {"nc", "[kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::nr, {"nr", "[kg-1]"}});
}

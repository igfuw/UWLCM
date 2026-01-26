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
  using thrust_real_t = setup::real_t;

  po::options_description opts("Double-moment bulk microphysics options");
  opts.add_options()
    ("acti", po::value<bool>()->default_value(rt_params.cloudph_opts.acti) , "activation (on/off)")
    ("cond", po::value<bool>()->default_value(rt_params.cloudph_opts.cond) , "condensation (on/off)")
    ("accr", po::value<bool>()->default_value(rt_params.cloudph_opts.accr) , "accretion (on/off)")
    ("acnv", po::value<bool>()->default_value(rt_params.cloudph_opts.acnv) , "autoconversion (on/off)")
    ("sedi", po::value<bool>()->default_value(rt_params.cloudph_opts.sedi) , "sedimentation (on/off)")
    ("acnv_A", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_A), "parameter in autoconversion rate formulae")
    ("acnv_b", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_b), "parameter in autoconversion rate formulae")
    ("acnv_c", po::value<typename solver_t::real_t>()->default_value(rt_params.cloudph_opts.acnv_c), "parameter in autoconversion rate formulae")
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

  if (user_params.n1_stp * si::cubic_metres >= 0) {
    rt_params.cloudph_opts.dry_distros.push_back({
      .mean_rd = user_params.mean_rd1 / si::metres,
      .sdev_rd = user_params.sdev_rd1,
      .N_stp   = user_params.n1_stp * si::cubic_metres,
      .chem_b  = user_params.kappa1
    });
  } 
  if (user_params.n2_stp * si::cubic_metres >= 0) {
    rt_params.cloudph_opts.dry_distros.push_back({
      .mean_rd = user_params.mean_rd2 / si::metres,
      .sdev_rd = user_params.sdev_rd2,
      .N_stp   = user_params.n2_stp * si::cubic_metres,
      .chem_b  = user_params.kappa2
    });
  }
  if (user_params.n1_stp * si::cubic_metres < 0 && user_params.n2_stp * si::cubic_metres < 0) {
    rt_params.cloudph_opts.dry_distros.push_back({
      .mean_rd = case_ptr->mean_rd1 / si::metres,
      .sdev_rd = case_ptr->sdev_rd1,
      .N_stp   = user_params.case_n_stp_multiplier * case_ptr->n1_stp * si::cubic_metres,
      .chem_b  = case_ptr->kappa
    });
    rt_params.cloudph_opts.dry_distros.push_back({
      .mean_rd = case_ptr->mean_rd2 / si::metres,
      .sdev_rd = case_ptr->sdev_rd2,
      .N_stp   = user_params.case_n_stp_multiplier * case_ptr->n2_stp * si::cubic_metres,
      .chem_b  = case_ptr->kappa
    });
  }

  rt_params.cloudph_opts.th_dry = false;
  rt_params.cloudph_opts.const_p = true;

  // output variables
  rt_params.outvars.insert({solver_t::ix::rc, {"rc", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::rr, {"rr", "[kg kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::nc, {"nc", "[kg-1]"}});
  rt_params.outvars.insert({solver_t::ix::nr, {"nr", "[kg-1]"}});

  if(rt_params.aerosol_independent_of_rhod)
    std::cerr << "UWLCM warning: aerosol_independent_of_rhod has no effect on blk_2m microphysics" << std::endl;
  if(rt_params.aerosol_conc_factor.size()!=0)
    std::cerr << "UWLCM warning: aerosol_conc_factor has no effect on blk_2m microphysics" << std::endl;
}

/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "detail/user_params.hpp"

template<template<class...> class slvr, class ct_params_dim_micro, int n_dims>
void run_hlpr(bool piggy, bool sgs, const std::string &type, const int (&nps)[n_dims], const user_params_t &user_params);

/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "run_hlpr.cpp"

// explicit instantiation
template
void run_hlpr<slvr_dry, ct_params_2D_dry, 2>(bool piggy, bool sgs, const std::string &type, const int (&nps)[2], const user_params_t &user_params);

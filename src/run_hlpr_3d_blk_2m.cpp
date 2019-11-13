/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "run_hlpr.cpp"

// explicit instantiation
template
void run_hlpr<slvr_blk_2m, ct_params_3D_blk_2m, 3>(bool piggy, bool sgs, const std::string &type, const int (&nps)[3], const user_params_t &user_params);

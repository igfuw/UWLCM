/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "run_hlpr.cpp"

// explicit instantiation
template
void run_hlpr<slvr_blk_1m_ice, ct_params_2D_blk_1m_ice, 2>(bool piggy, bool sgs, const std::string &type, const int (&nps)[2], const user_params_t &user_params);

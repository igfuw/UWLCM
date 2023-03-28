#pragma once
#include "../slvr_lgrngn.hpp"

template <class ct_params_t>
virtual void slvr_lgrngn<ct_params_t>::init_prtcls()
{
  assert(this->rank == 0);

  // temporary array of densities - prtcls cant be init'd with 1D profile
  typename parent_t::arr_t rhod(this->mem->advectee(ix::th).shape()); // TODO: replace all rhod arrays with this->mem->G
  rhod = (*params.rhod)(this->vert_idx);

  // temporary array of pressure - prtcls cant be init'd with 1D profile
  typename parent_t::arr_t p_e(this->mem->advectee(ix::th).shape()); 
  p_e = (*params.p_e)(this->vert_idx);

  prtcls->init(
    make_arrinfo(this->mem->advectee(ix::th)),
    make_arrinfo(this->mem->advectee(ix::rv)),
    make_arrinfo(rhod),
    make_arrinfo(p_e)
  ); 
}

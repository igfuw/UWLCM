#pragma once
#include <libcloudph++/blk_1m/rhs_cellwise.hpp>
//#include "../slvr_blk_1m.hpp"

template <class ct_params_t>
void slvr_blk_1m_ice_common<ct_params_t>::update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) 
{
  // store rl for buoyancy
  this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk) + this->state(ix::ria)(this->ijk) + this->state(ix::rib)(this->ijk);
  this->mem->barrier();

  parent_t::parent_t::update_rhs(rhs, dt, at);

  // zero-out precipitation rate, will be set in columnwise
  if(at ==0)
  {
    this->precipitation_rate(this->ijk) = 0;
    iceA_precipitation_rate(this->ijk) = 0;
    iceB_precipitation_rate(this->ijk) = 0;
  }

  this->mem->barrier();

  // cell-wise
  // TODO: rozne cell-wise na n i n+1 ?
  if(at == 0)
  {
    auto
      dot_th = rhs.at(ix::th)(this->ijk),
      dot_rv = rhs.at(ix::rv)(this->ijk),
      dot_rc = rhs.at(ix::rc)(this->ijk),
      dot_rr = rhs.at(ix::rr)(this->ijk),
      dot_ria = rhs.at(ix::ria)(this->ijk),
      dot_rib = rhs.at(ix::rib)(this->ijk);
    const auto
      th   = this->state(ix::th)(this->ijk),
      rv   = this->state(ix::rv)(this->ijk),
      rc   = this->state(ix::rc)(this->ijk),
      rr   = this->state(ix::rr)(this->ijk),
      ria   = this->state(ix::ria)(this->ijk),
      rib   = this->state(ix::rib)(this->ijk),
      rhod = (*this->mem->G)(this->ijk),
      &p_e_arg = this->p_e(this->ijk);
    libcloudphxx::blk_1m::rhs_cellwise_nwtrph_ice<real_t>(
        params.cloudph_opts,
        dot_th, dot_rv, dot_rc, dot_rr, dot_ria, dot_rib,
        rhod, p_e_arg, th, rv, rc, rr, ria, rib,
        dt
    );

    this->mem->barrier();

    nancheck(rhs.at(ix::th)(this->ijk), "RHS of th after rhs_cellwise");
    nancheck(rhs.at(ix::rv)(this->ijk), "RHS of rv after rhs_cellwise");
    nancheck2(rhs.at(ix::rc)(this->ijk), this->state(ix::rc)(this->ijk), "RHS of rc after rhs_cellwise (+ output of rc)");
    nancheck2(rhs.at(ix::rr)(this->ijk), this->state(ix::rr)(this->ijk), "RHS of rr after rhs_cellwise (+ output of rr)");
    nancheck2(rhs.at(ix::ria)(this->ijk), this->state(ix::ria)(this->ijk), "RHS of ria after rhs_cellwise (+ output of ria)");
    nancheck2(rhs.at(ix::rib)(this->ijk), this->state(ix::rib)(this->ijk), "RHS of rib after rhs_cellwise (+ output of rib)");
    //nancheck(rhs.at(ix::rr)(this->ijk), "RHS of rr after rhs_cellwise");
  }

  // forcing
  switch (at)
  {
    // for eulerian integration or used to init trapezoidal integration
    case (0):
    {
      // ---- cloud water sources ----
      parent_t::rc_src();
      rhs.at(ix::rc)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rc)(this->ijk);
      nancheck(rhs.at(ix::rc)(this->ijk), "RHS of rc after rc_src");

      // ---- rain water sources ----
      parent_t::rr_src();
      rhs.at(ix::rr)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rr)(this->ijk);
      nancheck(rhs.at(ix::rr)(this->ijk), "RHS of rr after rr_src");

      // ---- iceA sources ----
      ria_src();
      rhs.at(ix::ria)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rc)(this->ijk);
      nancheck(rhs.at(ix::ria)(this->ijk), "RHS of ria after ria_src");

      // ---- iceB sources ----
      rib_src();
      rhs.at(ix::rib)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rc)(this->ijk);
      nancheck(rhs.at(ix::rib)(this->ijk), "RHS of rib after rib_src");

  
      // when using explicit turbulence model add subgrid forces to rc and rr
      // (th and rv were already applied in slvr_sgs)
      if (ct_params_t::sgs_scheme != libmpdataxx::solvers::iles)
      {
        this->sgs_scalar_forces({ix::rc, ix::rr, ix::ria, ix::rib});
        nancheck(rhs.at(ix::rc)(this->ijk), "RHS of rc after sgs_scalar_forces");
        nancheck(rhs.at(ix::rr)(this->ijk), "RHS of rr after sgs_scalar_forces");
        nancheck(rhs.at(ix::ria)(this->ijk), "RHS of ria after sgs_scalar_forces");
        nancheck(rhs.at(ix::rib)(this->ijk), "RHS of rib after sgs_scalar_forces");
      }
      
      break;
    }

    case (1):
    // trapezoidal rhs^n+1
    {
    /*
      // ---- cloud water sources ----
      rc_src();
      rhs.at(ix::rc)(this->ijk) += this->alpha(this->ijk) + this->beta(this->ijk) * this->state(ix::rc)(this->ijk) / (1. - 0.5 * this->dt * this->beta(this->ijk));

      // ---- rain water sources ----
      rr_src();
      rhs.at(ix::rr)(this->ijk) += this->alpha(this->ijk) + this->beta(this->ijk) * this->state(ix::rr)(this->ijk) / (1. - 0.5 * this->dt * this->beta(this->ijk));

*/
      break;
    }
  }
  this->mem->barrier();
}

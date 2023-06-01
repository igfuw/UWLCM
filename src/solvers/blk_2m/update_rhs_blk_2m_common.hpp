#pragma once
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include "../slvr_blk_2m.hpp"

template <class ct_params_t>
void slvr_blk_2m_common<ct_params_t>::update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) 
{
  // store rl for buoyancy
  this->r_l(this->ijk) = this->state(ix::rc)(this->ijk) + this->state(ix::rr)(this->ijk);

  parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

  // zero-out precipitation rate, will be set in columnwise
  if(at == 0) 
  {
    rr_flux(this->ijk) = 0;
    nr_flux(this->ijk) = 0;
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
      dot_nc = rhs.at(ix::nc)(this->ijk),
      dot_rr = rhs.at(ix::rr)(this->ijk),
      dot_nr = rhs.at(ix::nr)(this->ijk);
     const auto
      th     = this->state(ix::th)(this->ijk),
      rv     = this->state(ix::rv)(this->ijk),
      rc     = this->state(ix::rc)(this->ijk),
      nc     = this->state(ix::nc)(this->ijk),
      rr     = this->state(ix::rr)(this->ijk),
      nr     = this->state(ix::nr)(this->ijk),
      rhod   = (*this->mem->G)(this->ijk),
      &p_e_arg = p_e(this->ijk); //TODO: use const pressure in blk_2m
      nancheck(nc, "nc before blk_2m rhs_cellwise call");
      negtozero(this->state(ix::nc)(this->ijk), "nc before blk_2m rhs_cellwise call");
      negcheck(nc, "nc before blk_2m rhs_cellwise call");
      nancheck(nr, "nr before blk_2m rhs_cellwise call");
      negtozero(this->state(ix::nr)(this->ijk), "nr before blk_2m rhs_cellwise call");
      negcheck(nr, "nr before blk_2m rhs_cellwise call");
     libcloudphxx::blk_2m::rhs_cellwise<real_t>(
      params.cloudph_opts, dot_th,  dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
      rhod, th,     rv,     rc,     nc,     rr,     nr,
      this->dt, true, p_e_arg
    );
      nancheck(nc, "nc after blk_2m rhs_cellwise call");
      negcheck(nc, "nc after blk_2m rhs_cellwise call");
      nancheck(nc, "nr after blk_2m rhs_cellwise call");
      negcheck(nc, "nr after blk_2m rhs_cellwise call");

//      negcheck(dot_nc, "dot_nc after blk_2m rhs_cellwise call");
  }

  // forcing
  switch (at)
  {
    // for eulerian integration or used to init trapezoidal integration
    case (0):
    {
      // ---- cloud water sources ----
      rc_src();
      //this->common_water_src(ix::rc, params.rc_src);
      rhs.at(ix::rc)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rc)(this->ijk);
      
      nc_src();
      //this->common_water_src(ix::nc, params.nc_src);
      rhs.at(ix::nc)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::nc)(this->ijk);

      // ---- rain water sources ----
      rr_src();
      //this->common_water_src(ix::rr, params.rr_src);
      rhs.at(ix::rr)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::rr)(this->ijk);
      
      nr_src();
      //this->common_water_src(ix::nr, params.nr_src);
      rhs.at(ix::nr)(this->ijk) += this->alpha(this->ijk);// + this->beta(this->ijk) * this->state(ix::nr)(this->ijk);

      // when using explicit turbulence model add subgrid forces to rc and rr
      // (th and rv were already applied in slvr_sgs)
      if (ct_params_t::sgs_scheme != libmpdataxx::solvers::iles)
      {
        this->sgs_scalar_forces({ix::rc, ix::rr, ix::nc, ix::nr});
        nancheck(rhs.at(ix::rc)(this->ijk), "RHS of rc after sgs_scalar_forces");
        nancheck(rhs.at(ix::rr)(this->ijk), "RHS of rr after sgs_scalar_forces");
        nancheck(rhs.at(ix::nc)(this->ijk), "RHS of nc after sgs_scalar_forces");
        nancheck(rhs.at(ix::nr)(this->ijk), "RHS of nr after sgs_scalar_forces");
      }

      break;
    }

    case (1):
    // trapezoidal rhs^n+1
    {
      break;
    }
  }

  // assure that we do not remove more cloud/rain water than there is!
  // this could happen due to blk_2m microphyscs combined with other forcings (e.g. subsidence)
  auto
    dot_rc = rhs.at(ix::rc)(this->ijk),
    dot_rr = rhs.at(ix::rr)(this->ijk),
    dot_nc = rhs.at(ix::nc)(this->ijk),
    dot_nr = rhs.at(ix::nr)(this->ijk);

  const auto
    rc     = this->state(ix::rc)(this->ijk),
    rr     = this->state(ix::rr)(this->ijk),
    nc     = this->state(ix::nc)(this->ijk),
    nr     = this->state(ix::nr)(this->ijk);

  dot_rc = where(dot_rc * dt <= -rc, -rc / dt, dot_rc);
  dot_rr = where(dot_rr * dt <= -rr, -rr / dt, dot_rr);
  dot_nc = where(dot_nc * dt <= -nc, -nc / dt, dot_nc);
  dot_nr = where(dot_nr * dt <= -nr, -nr / dt, dot_nr);


  this->mem->barrier();
  if(this->rank == 0)
  {
    nancheck(rhs.at(ix::rc)(this->domain), "RHS of rc after rhs_update");
    nancheck(rhs.at(ix::rr)(this->domain), "RHS of rr after rhs_update");
    nancheck(rhs.at(ix::nc)(this->domain), "RHS of nc after rhs_update");
    nancheck(rhs.at(ix::nr)(this->domain), "RHS of nr after rhs_update");
  }
}

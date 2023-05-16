#pragma once
#include "blk_2m/slvr_blk_2m_common.hpp"
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

template <class ct_params_t, class enableif = void>
class slvr_blk_2m
{};

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?

// 2D version
template <class ct_params_t>
class slvr_blk_2m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_2m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_2m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_2m(
    typename parent_t::ctor_args_t args,
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p)
  {}

  protected:
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at);

    this->mem->barrier();
    if(at == 0)
    {
      // column-wise
      for (int i = this->i.first(); i <= this->i.last(); ++i)
      {
        auto
          dot_rr = this->rr_flux(i, this->j),
          dot_nr = this->nr_flux(i, this->j);
        const auto
          rhod   = (*this->mem->G)(i, this->j),
          rr     = this->state(parent_t::ix::rr)(i, this->j),
          nr     = this->state(parent_t::ix::nr)(i, this->j);
        this->liquid_puddle += -libcloudphxx::blk_2m::rhs_columnwise<real_t>(
          this->params.cloudph_opts, dot_rr, dot_nr, rhod, rr, nr, this->dt, this->params.dz
        );
      }
      rhs.at(parent_t::ix::rr)(this->ijk) += this->rr_flux(this->ijk);
      rhs.at(parent_t::ix::nr)(this->ijk) += this->nr_flux(this->ijk);

      this->mem->barrier();
      if(this->rank == 0)
      {
        nancheck(rhs.at(parent_t::ix::rr)(this->domain), "RHS of rr after rhs_update");
        nancheck(rhs.at(parent_t::ix::nr)(this->domain), "RHS of nr after rhs_update");
      }
    }
  }
};

// 3D version
template <class ct_params_t>
class slvr_blk_2m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_blk_2m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_2m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  // ctor
  slvr_blk_2m(
    typename parent_t::ctor_args_t args,
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p)
  {}

  protected:
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at);

    this->mem->barrier();
    if(at == 0)
    {
      // column-wise
      for (int i = this->i.first(); i <= this->i.last(); ++i)
        for (int j = this->j.first(); j <= this->j.last(); ++j)
        {
          auto
          dot_rr = this->rr_flux(i, j, this->k),
          dot_nr = this->nr_flux(i, j, this->k);    
          const auto
          rhod   = (*this->mem->G)(i, j, this->k),
          rr     = this->state(parent_t::ix::rr)(i, j, this->k),
          nr     = this->state(parent_t::ix::nr)(i, j, this->k);
          this->liquid_puddle += -libcloudphxx::blk_2m::rhs_columnwise<real_t>(
            this->params.cloudph_opts, dot_rr, dot_nr, rhod, rr, nr, this->dt, this->params.dz
          );
        }
        rhs.at(parent_t::ix::rr)(this->ijk) += this->rr_flux(this->ijk);
        rhs.at(parent_t::ix::nr)(this->ijk) += this->nr_flux(this->ijk);

      this->mem->barrier();
      if(this->rank == 0)
      {
        nancheck(rhs.at(parent_t::ix::rr)(this->domain), "RHS of rr after rhs_update");
        nancheck(rhs.at(parent_t::ix::nr)(this->domain), "RHS of nr after rhs_update");
      }
    }
  }
};

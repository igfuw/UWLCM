#pragma once
#include "blk_1m/slvr_blk_1m_common.hpp"
#include <libcloudph++/blk_1m/rhs_columnwise.hpp>

template <class ct_params_t, class enableif = void>
class slvr_blk_1m
{};

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?

/**
 * \class slvr_blk_1m_2D
 * @brief Solver for 2D simulations using the bulk 1-moment scheme.
 *
 * @tparam ct_params_t Compile-time parameters (must define n_dims == 2).
 */
template <class ct_params_t>
class slvr_blk_1m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_blk_1m_common<ct_params_t>
{
  public:
  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  /**
 * @brief Constructor.
 *
 * @param args Arguments forwarded to the parent solver.
 * @param p Runtime parameters.
 */
  slvr_blk_1m(
    typename parent_t::ctor_args_t args,
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p)
  {}

  protected:
  /**
   * @brief Updates the right-hand side (RHS) tendencies for 2D case.
   * Includes condensation and precipitation source terms.
   *
   * @param rhs Reference to RHS array vector to be updated.
   * @param dt Time step size.
   * @param at Substep index (applied only for at == 0).
   */
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    this->mem->barrier();
    if(at == 0)
    {
      // column-wise
      for (int i = this->i.first(); i <= this->i.last(); ++i)
      {
        auto 
          precipitation_rate_arg = this->precipitation_rate(i, this->j);
        const auto 
          rhod   = (*this->mem->G)(i, this->j),
          rr     = this->state(parent_t::ix::rr)(i, this->j);
        this->liquid_puddle += - libcloudphxx::blk_1m::rhs_columnwise<real_t>(
          this->params.cloudph_opts, precipitation_rate_arg, rhod, rr, this->params.dz
        );
      }
      rhs.at(parent_t::ix::rr)(this->ijk) += this->precipitation_rate(this->ijk);

      nancheck(rhs.at(parent_t::ix::rr)(this->ijk), "RHS of rr after rhs_update");
      this->mem->barrier();
    }
  }
};

/**
 *\class slvr_blk_1m_3D
 * @brief Solver for 3D simulations using the bulk 1-moment scheme.
 *
 * @tparam ct_params_t Compile-time parameters (must define n_dims == 3).
 */
template <class ct_params_t>
class slvr_blk_1m<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_blk_1m_common<ct_params_t>
{

  public:
  using parent_t = slvr_blk_1m_common<ct_params_t>;
  using real_t = typename ct_params_t::real_t;

  /**
 * @brief Constructor.
 *
 * @param args Arguments forwarded to the parent solver.
 * @param p Runtime parameters.
 */
  slvr_blk_1m(
    typename parent_t::ctor_args_t args,
    const typename parent_t::rt_params_t &p
  ) :
    parent_t(args, p)
  {}

  protected:
  /**
 * @brief Updates the right-hand side (RHS) tendencies for 3D case.
 *
 * Includes condensation and precipitation source terms.
 *
 * @param rhs Reference to RHS array vector to be updated.
 * @param dt Time step size.
 * @param at Substep index (applied only for at == 0).
 */
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at); // shouldnt forcings be after condensation to be consistent with lgrngn solver?

    this->mem->barrier();
    if(at == 0)
    {
      // column-wise
      for (int i = this->i.first(); i <= this->i.last(); ++i)
        for (int j = this->j.first(); j <= this->j.last(); ++j)
        {
          auto 
          precipitation_rate_arg = this->precipitation_rate(i, j, this->k);
          const auto 
          rhod   = (*this->mem->G)(i, j, this->k),
          rr     = this->state(parent_t::ix::rr)(i, j, this->k);
          this->liquid_puddle += - libcloudphxx::blk_1m::rhs_columnwise<real_t>(
            this->params.cloudph_opts, precipitation_rate_arg, rhod, rr, this->params.dz
          );
        }
      rhs.at(parent_t::ix::rr)(this->ijk) += this->precipitation_rate(this->ijk);

      nancheck(rhs.at(parent_t::ix::rr)(this->ijk), "RHS of rr after rhs_update");
      this->mem->barrier();
    }
  }
};

#pragma once
#include "slvr_common.hpp"

template <class ct_params_t, class enableif = void>
class slvr_dim
{};

using libmpdataxx::arakawa_c::h;

// 2D version 
template <class ct_params_t>
class slvr_dim<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  protected:
  // inject dimension-independent ranges
  idx_t<2> domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]});
  rng_t horizontal_domain = this->mem->grid_size[0];
  idx_t<2> Cx_domain = idx_t<2>({this->mem->grid_size[0]^h, this->mem->grid_size[1]});
  idx_t<2> Cy_domain;
  idx_t<2> Cz_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]^h});

  void vert_grad(typename parent_t::arr_t &in, typename parent_t::arr_t &out, setup::real_t dz)
  {
    for (auto &bc : this->bcs[1]) bc->fill_halos_sclr(in, this->i, false);
    out(this->i, this->j) = ( in(this->i, this->j) - in(this->i, this->j+1)) / dz;
    // top and bottom cells are two times lower
    out(this->i, 0) *= 2; 
    out(this->i, this->j.last()) *= 2; 
  }

  enum {vert_dim = 1};

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p)
  {}
};

// 3D version
template <class ct_params_t>
class slvr_dim<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  protected:
  // inject dimension-independent ranges
  idx_t<3> domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<2> horizontal_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]});
  idx_t<3> Cx_domain = idx_t<3>({this->mem->grid_size[0]^h, this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<3> Cy_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1]^h, this->mem->grid_size[2]});
  idx_t<3> Cz_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]^h});

  enum {vert_dim = 2};

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args, 
    typename parent_t::rt_params_t const &p
  ) : 
    parent_t(args, p)
  {}
};


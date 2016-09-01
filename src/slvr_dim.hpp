#pragma once
#include "slvr_common.hpp"

template <class ct_params_t, class enableif = void>
class slvr_dim
{};

using libmpdataxx::arakawa_c::h;

// 2D version - inject dimension-independent ranges
template <class ct_params_t>
class slvr_dim<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  protected:
  idx_t<2> domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]});
  idx_t<2> Cx_domain = idx_t<2>({this->mem->grid_size[0]^h, this->mem->grid_size[1]});
  idx_t<2> Cy_domain;
  idx_t<2> Cz_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]^h});

  blitz::TinyVector<float, 2> domain_size = {
    this->mem->grid_size[0].length(),
    this->mem->grid_size[1].length()
  }; // TODO: could be replaced by mpdata's rt_params::grid_size ?

  enum {vert_dim = 1};

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p)
  {}
};

// 3D version - inject dimension-independent ranges
template <class ct_params_t>
class slvr_dim<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  protected:
  idx_t<3> domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<3> Cx_domain = idx_t<3>({this->mem->grid_size[0]^h, this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<3> Cy_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1]^h, this->mem->grid_size[2]});
  idx_t<3> Cz_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]^h});

  blitz::TinyVector<float, 3> domain_size = {
    this->mem->grid_size[0].length(),
    this->mem->grid_size[1].length(),
    this->mem->grid_size[2].length()
  }; // TODO: could be replaced by mpdata's rt_params::grid_size ?

  enum {vert_dim = 2};

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args, 
    typename parent_t::rt_params_t const &p
  ) : 
    parent_t(args, p)
  {}
};


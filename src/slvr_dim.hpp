#pragma once
#include "slvr_common.hpp"

template <class ct_params_t, class enableif = void>
class slvr_dim
{};

// 2D version - inject dimension-independent ranges
template <class ct_params_t>
class slvr_dim<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  private:
  int nx = this->mem->grid_size[0].length();
  int nz = this->mem->grid_size[1].length();

  protected:
  blitz::RectDomain<2> domain = blitz::RectDomain<2>({blitz::Range(0,nx-1), blitz::Range(0,nz-1)});
  blitz::RectDomain<2> Cx_domain = 
    blitz::RectDomain<2>({this->mem->grid_size[0]^libmpdataxx::arakawa_c::h, this->mem->grid_size[1]});
  blitz::RectDomain<2> Cz_domain = 
    blitz::RectDomain<2>({this->mem->grid_size[0], this->mem->grid_size[1]^libmpdataxx::arakawa_c::h});
  blitz::TinyVector<float, 2> domain_size = {nx, nz};

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

  private: 
  int nx = this->mem->grid_size[0].length();
  int ny = this->mem->grid_size[1].length();
  int nz = this->mem->grid_size[2].length();

  protected:
  blitz::RectDomain<3> domain = blitz::RectDomain<3>({blitz::Range(0,nx-1), blitz::Range(0,ny-1), blitz::Range(0,nz-1)});
  blitz::TinyVector<float, 3> domain_size = {nx, ny, nz};

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args, 
    typename parent_t::rt_params_t const &p
  ) : 
    parent_t(args, p)
  {}
};


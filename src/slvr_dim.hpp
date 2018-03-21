#pragma once

#include "slvr_piggy.hpp"
#include <libmpdata++/formulae/arakawa_c.hpp>

// custom 3D idxperm that accepts idx_t; todo: make it part of libmpdata?
namespace libmpdataxx
{
  namespace idxperm
  {
    template<int d>
    inline idx_t<3> pi(const rng_t &rng, const idx_t<2> &idx) { return pi<d>(rng, idx[0], idx[1]); }

    template<int d>
    inline idx_t<3> pi(const int &i, const idx_t<2> &idx) { return pi<d>(rng_t(i,i), idx); }
  };
};

template <class ct_params_t, class enableif = void>
class slvr_dim
{};

using libmpdataxx::arakawa_c::h;
using namespace libmpdataxx; // TODO: get rid of it?
using namespace arakawa_c;

// 2D version 
template <class ct_params_t>
class slvr_dim<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 2 >::type
> : public slvr_piggy<ct_params_t> 
{
  using parent_t = slvr_piggy<ct_params_t>;
  using ix = typename ct_params_t::ix;

  protected:
  using arr_sub_t = blitz::Array<setup::real_t, 1>;
  // inject dimension-independent ranges
  idx_t<2> domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]});
  rng_t hrzntl_domain = this->mem->grid_size[0];
  rng_t hrzntl_subdomain = this->i;
  idx_t<2> Cx_domain = idx_t<2>({this->mem->grid_size[0]^h, this->mem->grid_size[1]});
  idx_t<2> Cy_domain;
  idx_t<2> Cz_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]^h});

  blitz::TinyVector<int, 2> zero = blitz::TinyVector<int, 2>({0,0});
  blitz::secondIndex vert_idx;
  libmpdataxx::arrvec_t<arr_sub_t> vip_ground;
  std::set<int> hori_vel = std::set<int>{ix::u};

  void vert_grad_fwd(typename parent_t::arr_t &in, typename parent_t::arr_t &out, setup::real_t dz)
  {
    in(this->i, this->j.last() + 1) = in(this->i, this->j.last()); 
    out(this->i, this->j) = ( in(this->i, this->j+1) - in(this->i, this->j)) / dz;
    // top and bottom cells are two times lower
    out(this->i, 0) *= 2; 
    out(this->i, this->j.last()) *= 2; 
  }

  void vert_grad_cnt(typename parent_t::arr_t &in, typename parent_t::arr_t &out, setup::real_t dz)
  {
    in(this->i, this->j.last() + 1) = in(this->i, this->j.last()); 
    in(this->i, this->j.first() - 1) = in(this->i, this->j.first()); 
    out(this->i, this->j) = ( in(this->i, this->j+1) - in(this->i, this->j-1)) / 2./ dz;
    // top and bottom cells are two times lower
    out(this->i, 0) *= 2; 
    out(this->i, this->j.last()) *= 2; 
  }

  void smooth(typename parent_t::arr_t &in, typename parent_t::arr_t &out)
  {
    this->xchng_sclr(in, this->ijk);

    out(this->i, this->j) = (4 * in(this->i, this->j) + 
                                 in(this->i + 1, this->j) + in(this->i - 1, this->j) +
                                 in(this->i, this->j + 1) + in(this->i, this->j - 1) 
                            ) / 8.;
    this->mem->barrier();
  }

  auto calc_U_ground() 
    return_macro(,
    abs(this->state(ix::vip_i)(this->i, 0).reindex({0}))
  )

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p)
  {
    vip_ground.push_back(new arr_sub_t(this->state(ix::vip_i)(this->i, 0).reindex({0})));
  }
};

// 3D version
template <class ct_params_t>
class slvr_dim<
  ct_params_t,
  typename std::enable_if<ct_params_t::n_dims == 3 >::type
> : public slvr_piggy<ct_params_t> 
{
  using parent_t = slvr_piggy<ct_params_t>;
  using ix = typename ct_params_t::ix;

  protected:
  using arr_sub_t = blitz::Array<setup::real_t, 2>;
  // inject dimension-independent ranges
  idx_t<3> domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<2> hrzntl_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]});
  idx_t<2> hrzntl_subdomain = idx_t<2>({this->i, this->j});
  idx_t<3> Cx_domain = idx_t<3>({this->mem->grid_size[0]^h, this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<3> Cy_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1]^h, this->mem->grid_size[2]});
  idx_t<3> Cz_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]^h});

  blitz::TinyVector<int, 3> zero = blitz::TinyVector<int, 3>({0,0,0});
  blitz::thirdIndex vert_idx;
  libmpdataxx::arrvec_t<arr_sub_t> vip_ground;
  std::set<int> hori_vel = std::set<int>{ix::u, ix::v};

  void vert_grad_fwd(typename parent_t::arr_t &in, typename parent_t::arr_t &out, setup::real_t dz)
  {
    in(this->i, this->j, this->k.last() + 1) = in(this->i, this->j, this->k.last()); 
    out(this->i, this->j, this->k) = ( in(this->i, this->j, this->k+1) - in(this->i, this->j, this->k)) / dz;
    // top and bottom cells are two times lower
    out(this->i, this->j, 0) *= 2; 
    out(this->i, this->j, this->k.last()) *= 2; 
  }

  void vert_grad_cnt(typename parent_t::arr_t &in, typename parent_t::arr_t &out, setup::real_t dz)
  {
    in(this->i, this->j, this->k.last() + 1) = in(this->i, this->j, this->k.last()); 
    in(this->i, this->j, this->k.first() - 1) = in(this->i, this->j, this->k.first()); 
    out(this->i, this->j, this->k) = ( in(this->i, this->j, this->k+1) - in(this->i, this->j, this->k-1)) / 2./ dz;
    // top and bottom cells are two times lower
    out(this->i, this->j, 0) *= 2; 
    out(this->i, this->j, this->k.last()) *= 2; 
  }

  void smooth(typename parent_t::arr_t &in, typename parent_t::arr_t &out)
  {
    this->xchng_sclr(in, this->ijk); 
    out(this->i, this->j, this->k) = (6 * in(this->i, this->j, this->k) + 
                                      in(this->i + 1, this->j, this->k) + in(this->i - 1, this->j, this->k) +
                                      in(this->i, this->j + 1, this->k) + in(this->i, this->j - 1, this->k) +
                                      in(this->i, this->j, this->k + 1) + in(this->i, this->j, this->k - 1)
                                     ) / 12.;
    this->mem->barrier();
  }

  auto calc_U_ground() 
    return_macro(,
    sqrt(pow2(this->state(ix::vip_i)(this->i, this->j, 0).reindex({0,0})) + pow2(this->state(ix::vip_j)(this->i, this->j, 0).reindex({0,0})))
  )

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args, 
    typename parent_t::rt_params_t const &p
  ) : 
    parent_t(args, p)
  {
    vip_ground.push_back(new arr_sub_t(this->state(ix::vip_i)(this->i, this->j, 0).reindex({0,0})));
    vip_ground.push_back(new arr_sub_t(this->state(ix::vip_j)(this->i, this->j, 0).reindex({0,0})));
  }
};


#pragma once

#include "slvr_piggy.hpp"
#include <libmpdata++/formulae/arakawa_c.hpp>
#include <libmpdata++/formulae/nabla_formulae.hpp>

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
  public:
  using parent_t = slvr_piggy<ct_params_t>;
  using ix = typename ct_params_t::ix;
  using arr_sub_t = blitz::Array<setup::real_t, 1>;

  protected:
  // inject dimension-independent ranges
  idx_t<2> domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]});
  rng_t hrzntl_domain = this->mem->grid_size[0];
  rng_t hrzntl_subdomain = this->i;
  idx_t<2> Cx_domain = idx_t<2>({this->mem->grid_size[0]^h, this->mem->grid_size[1]}); // libcloudphxx requires courants with a halo of 2 in the x direction
  idx_t<2> Cy_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]^h}); // just fill in with Cz_domain to avoid some asserts
  idx_t<2> Cz_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]^h});


  blitz::TinyVector<int, 2> zero = blitz::TinyVector<int, 2>({0,0});
  blitz::TinyVector<int, 1> zero_plane = blitz::TinyVector<int, 1>({0});
  blitz::TinyVector<int, 2> origin = blitz::TinyVector<int, 2>({this->i.first(), this->j.first()});

  blitz::secondIndex vert_idx;
  const rng_t &vert_rng = this->j;
  std::set<int> hori_vel = std::set<int>{ix::u};

  idx_t<2> hrzntl_slice(int k)
  {
      return idx_t<2>({this->i, rng_t(k, k)});
  }
  
  auto hrzntl_slice(const typename parent_t::arr_t &a, int k)
  {
      return blitz::safeToReturn(a(idx_t<2>({this->i, rng_t(k, k)})) + 0);
  }

  void vert_grad_fwd(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t dz)
  {
    // extrapolate upward, top cell is two times lower
    in(this->i, this->j.last() + 1) = 1.5*in(this->i, this->j.last()) - .5 * in(this->i, this->j.last()-1); 
    out(this->i, this->j) = ( in(this->i, this->j+1) - in(this->i, this->j)) / dz;
    // top nad bottom cells are two times lower
    out(this->i, this->j.last()) *= 2; 
    out(this->i, 0) *= 2; 
  }

  void vert_grad_cnt(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t dz)
  {
    in(this->i, this->j.last() + 1) = in(this->i, this->j.last()); 
    in(this->i, this->j.first() - 1) = in(this->i, this->j.first()); 
    out(this->i, this->j) = ( in(this->i, this->j+1) - in(this->i, this->j-1)) / 2./ dz;
    // top and bottom cells are two times lower
    out(this->i, 0) *= 2; 
    // set to 0 at top level to have no subsidence there - TODO: it messes with other possible uses of this function 
    out(this->i, this->j.last()) = 0; 
  }
  
  void vert_grad_cmpct(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t dz)
  {
    in(this->i, this->j.last() + 1) = in(this->i, this->j.last());
    in(this->i, this->j.first() - 1) = in(this->i, this->j.first());
    out(this->i, this->j + h) = libmpdataxx::formulae::nabla::grad_cmpct<1>(in, this->j, this->i, dz);
  }
  
  void vert_aver_cmpct(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t coeff = 1)
  {
    // assumes filled halos
    using libmpdataxx::arakawa_c::h;
    out(this->i, this->j) = coeff * (in(this->i, this->j - h) + in(this->i, this->j + h)) / 2;
  }

  void smooth(typename parent_t::arr_t in, typename parent_t::arr_t out)
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
    abs(this->state(ix::vip_i)(this->hrzntl_slice(0))) // at 1st level, because 0-th level has no clear interpretation? 0-th is ground level, but with horizontal winds
  )

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
> : public slvr_piggy<ct_params_t> 
{
  public:
  using parent_t = slvr_piggy<ct_params_t>;
  using ix = typename ct_params_t::ix;
  using arr_sub_t = blitz::Array<setup::real_t, 2>;

  protected:
  // inject dimension-independent ranges
  idx_t<3> domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<2> hrzntl_domain = idx_t<2>({this->mem->grid_size[0], this->mem->grid_size[1]});
  idx_t<2> hrzntl_subdomain = idx_t<2>({this->i, this->j});
  idx_t<3> Cx_domain = idx_t<3>({this->mem->grid_size[0]^h, this->mem->grid_size[1], this->mem->grid_size[2]});
  idx_t<3> Cy_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1]^h, this->mem->grid_size[2]});
  idx_t<3> Cz_domain = idx_t<3>({this->mem->grid_size[0], this->mem->grid_size[1], this->mem->grid_size[2]^h});

  blitz::TinyVector<int, 3> zero = blitz::TinyVector<int, 3>({0,0,0});
  blitz::TinyVector<int, 2> zero_plane = blitz::TinyVector<int, 2>({0,0});
  blitz::TinyVector<int, 3> origin = blitz::TinyVector<int, 3>({this->i.first(), this->j.first(), this->k.first()});

  blitz::thirdIndex vert_idx;
  const rng_t &vert_rng = this->k;
  std::set<int> hori_vel = std::set<int>{ix::u, ix::v};

  idx_t<3> hrzntl_slice(int k)
  {
      return idx_t<3>({this->i, this->j, rng_t(k, k)});
  }

  auto hrzntl_slice(const typename parent_t::arr_t &a, int k)
  {
      return blitz::safeToReturn(a(idx_t<3>({this->i, this->j, rng_t(k, k)})) + 0);
  }

  void vert_grad_fwd(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t dz)
  {
    // extrapolate upward
    in(this->i, this->j, this->k.last() + 1) = 1.5*in(this->i, this->j, this->k.last()) - 0.5*in(this->i, this->j, this->k.last()-1); 
    out(this->i, this->j, this->k) = ( in(this->i, this->j, this->k+1) - in(this->i, this->j, this->k)) / dz;
    // top and bottom cells are two times lower
    out(this->i, this->j, this->k.last()) *= 2; 
    out(this->i, this->j, 0) *= 2; 
  }

  void vert_grad_cnt(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t dz)
  {
    in(this->i, this->j, this->k.last() + 1) = in(this->i, this->j, this->k.last()); 
    in(this->i, this->j, this->k.first() - 1) = in(this->i, this->j, this->k.first()); 
    out(this->i, this->j, this->k) = ( in(this->i, this->j, this->k+1) - in(this->i, this->j, this->k-1)) / 2./ dz;
    // top and bottom cells are two times lower
    out(this->i, this->j, 0) *= 2; 
    //out(this->i, this->j, this->k.last()) *= 2; 
    out(this->i, this->j, this->k.last()) = 0; 
  }
  
  void vert_grad_cmpct(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t dz)
  {
    in(this->i, this->j, this->k.last() + 1) = in(this->i, this->j, this->k.last());
    in(this->i, this->j, this->k.first() - 1) = in(this->i, this->j, this->k.first());
    out(this->i, this->j, this->k + h) = libmpdataxx::formulae::nabla::grad_cmpct<2>(in, this->k, this->i, this->j, dz);
  }

  void vert_aver_cmpct(typename parent_t::arr_t in, typename parent_t::arr_t out, setup::real_t coeff = 1)
  {
    // assumes filled halos
    using libmpdataxx::arakawa_c::h;
    out(this->i, this->j, this->k) = coeff * (in(this->i, this->j, this->k - h) + in(this->i, this->j, this->k + h)) / 2;
  }

  void smooth(typename parent_t::arr_t in, typename parent_t::arr_t out)
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
    sqrt(pow2(this->state(ix::vip_i)(this->hrzntl_slice(0))) + pow2(this->state(ix::vip_j)(this->hrzntl_slice(0))))
  )

  // ctor
  slvr_dim(
    typename parent_t::ctor_args_t args, 
    typename parent_t::rt_params_t const &p
  ) : 
    parent_t(args, p)
    {}
};


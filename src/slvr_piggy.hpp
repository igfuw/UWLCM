#pragma once
#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

template <class ct_params_t, class enableif = void>
class slvr_piggy
{};

using namespace libmpdataxx; // TODO: get rid of it?

// driver
template <class ct_params_t>
class slvr_piggy<
  ct_params_t,
  typename std::enable_if<ct_params_t::piggy == 0 >::type
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs<ct_params_t>
  >
{
  protected:
  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs<ct_params_t>
  >;  

  std::ofstream f_cour_out; // output courant number file

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 
    // open file for out courants
    f_cour_out.open(this->outdir+"/courants_out.dat"); 
  }

  void hook_ante_step()
  {
    // save courant numbers
    if(this->rank==0)
    {
      for (int d = 0; d < parent_t::n_dims; ++d)
      {
        f_cour_out << this->vips()[d];
      }
    }
    parent_t::hook_ante_step();
  }

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p) {}
};


#pragma once
#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs_fra.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include "../detail/checknan.cpp"

template <class ct_params_t, class enableif = void>
class slvr_piggy
{};

using namespace libmpdataxx; // TODO: get rid of it?

constexpr int minhalo = 2; 

// driver
template <class ct_params_t>
class slvr_piggy<
  ct_params_t,
  typename std::enable_if<ct_params_t::piggy == 0 >::type
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>
  >
{
  private:
  bool save_vel; // should velocity field be stored for piggybacking
  setup::real_t prs_tol; // store a copy for output purposes

  protected:
  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs_sgs_fra<ct_params_t, minhalo>
  >;  

  std::ofstream f_vel_out; // file for velocity field

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 
    if(this->rank==0)
    {
      // open file for out vel
      if(save_vel)
      {
        try{
          f_vel_out.open(this->outdir+"/velocity_out.dat"); 
        }
        catch(...)
        {
          throw std::runtime_error("error opening velocity output file '{outdir}/velocity_out.dat'");
        }
      }
      this->record_aux_const("save_vel", "piggy", save_vel);  
      this->record_aux_const("rt_params prs_tol", "piggy", prs_tol);  
    }
    this->mem->barrier();
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes changes of velocity field due to vip_rhs_impl_fnlz()
    this->mem->barrier();
    // save velocity field
    if(this->rank==0 && save_vel)
    {
      for (int d = 0; d < parent_t::n_dims; ++d)
      {
        f_vel_out << this->state(this->vip_ixs[d]);
      }
    }
  }

  struct rt_params_t : parent_t::rt_params_t 
  {
    bool save_vel;

    // ctor
    rt_params_t()
    {
      po::options_description opts("Driver options"); 
      opts.add_options()
        ("save_vel", po::value<bool>()->default_value(false), "should velocity field be stored (for future piggybacking)")
      ;
      opts.add_options()
        ("prs_tol", po::value<setup::real_t>()->default_value(1e-6) , "pressure solver tolerance"); // not really related to piggybacking, but convenient to put here as it is the first solver to inherit from libmpdata++
      po::variables_map vm;
      handle_opts(opts, vm);
          
      save_vel = vm["save_vel"].as<bool>();
      this->prs_tol = vm["prs_tol"].as<setup::real_t>();
    }
  };

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    rt_params_t p
  ) :
    parent_t(args, p),
    save_vel(p.save_vel),
    prs_tol(p.prs_tol)
    {}
};


// piggybacker
template <class ct_params_t>
class slvr_piggy<
  ct_params_t,
  typename std::enable_if<ct_params_t::piggy == 1 >::type
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip<ct_params_t, minhalo>
  >
{

  protected:
  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs_vip<ct_params_t, minhalo>
  >;  

  private:
  typename parent_t::arr_t in_bfr; // input buffer for velocity
  std::string vel_in;
  
  protected:

  std::ifstream f_vel_in; // input velocity file

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    if(this->rank==0)
    {
      po::options_description opts("Piggybacker options"); 
      opts.add_options()
        ("vel_in", po::value<std::string>()->required(), "file with input velocities (for piggybacking)")
      ;
      po::variables_map vm;
      handle_opts(opts, vm);
          
      vel_in = vm["vel_in"].as<std::string>();
      std::cout << "piggybacking from: " << vel_in << std::endl;

      in_bfr.resize(this->state(this->vip_ixs[0]).shape());
      // open file for in vel
      // TODO: somehow check dimensionality of the input arrays
      try{
        f_vel_in.open(vel_in); 
      }
      catch(...)
      {
        throw std::runtime_error("error opening velocities input file defined by --vel_in");
      }
      this->record_aux_const("piggybacking", "piggy", "true");
      this->record_aux_const("vel_in", "piggy", vel_in); 
    }
    this->mem->barrier();
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); // do whatever
    this->mem->barrier(); //necessary?
    // read velo, overwrite any vel rhs
    if(this->rank==0)
    {
      using ix = typename ct_params_t::ix;

      for (int d = 0; d < parent_t::n_dims; ++d)
      {
        // read in through buffer, if done directly caused data races
        f_vel_in >> in_bfr;
        this->state(this->vip_ixs[d]) = in_bfr;
//std::cout << this->state(this->vip_ixs[d]);
      }
    }
    this->mem->barrier();
  }

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p) {}
};


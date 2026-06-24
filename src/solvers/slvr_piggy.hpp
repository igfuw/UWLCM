#pragma once
#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
#include "../detail/checknan.cpp"
#include <H5Cpp.h>
#include <libmpdata++/output/hdf5.hpp>

template <class ct_params_t, class enableif = void>
class slvr_piggy
{};

using namespace libmpdataxx; // TODO: get rid of it?

constexpr int minhalo = 1; 

// driver
/**
 * \class slvr_piggy_driver
 * @brief Solver class in driver mode (piggy == 0).
 *
 * This specialization of slvr_piggy handles running the simulation and storing velocity fields
 * for piggybacking. It inherits from @ref output::hdf5_xdmf with
 * @ref solvers::mpdata_rhs_vip_prs_sgs.
 *
 * @tparam ct_params_t Compile-time parameters.
 */
template <class ct_params_t>
class slvr_piggy<
  ct_params_t,
  typename std::enable_if<ct_params_t::piggy == 0 >::type
> : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
  >
{
  private:
  bool save_vel_flag; // should velocity field be stored for piggybacking
//  setup::real_t prs_tol; // store a copy for output purposes

  protected:
  using parent_t = output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs_sgs<ct_params_t, minhalo>
  >;  

  std::unique_ptr<H5::H5File> hdfpu_vel;

  void save_vel()
  {
    if(this->rank==0 && save_vel_flag)
    {
      hdfpu_vel.reset(new H5::H5File(this->outdir + "/velocities/" + this->hdf_name(this->base_name("velocity")), H5F_ACC_TRUNC
#if defined(USE_MPI)
          , H5P_DEFAULT, this->fapl_id
#endif
        ));

        for (int d = 0; d < parent_t::n_dims; ++d)
          this->record_aux_halo_hlpr(
            this->outvars[this->vip_ixs[d]].name,
            this->state(this->vip_ixs[d]),
            *hdfpu_vel
          );
    }
  }

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 
    if(this->rank==0)
    {
      this->record_aux_const("save_vel", "piggy", save_vel_flag);  
      this->record_aux_const("rt_params prs_tol", "piggy", this->prs_tol);  

      if (this->mem->distmem.rank() == 0 && save_vel_flag)
      {
        // creating the directory for velocity output
        boost::filesystem::create_directory(this->outdir + "/velocities/");
      }
    }
    this->mem->barrier();
    save_vel();
    this->mem->barrier();
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes changes of velocity field due to vip_rhs_impl_fnlz()
    this->mem->barrier();
    save_vel();
    this->mem->barrier();
  }

  struct rt_params_t : parent_t::rt_params_t 
  {
    bool save_vel_flag;

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
      handle_opts(opts, vm, false);
      save_vel_flag = vm["save_vel"].as<bool>();
      this->prs_tol = vm["prs_tol"].as<setup::real_t>();
    }
  };

  rt_params_t params;

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    save_vel_flag(p.save_vel_flag),
    params(p)
    //prs_tol(p.prs_tol)
    {}
};


// piggybacker
/**
 * \class slvr_piggy_piggybacker
 * @brief Solver class in piggybacker mode (piggy == 1).
 *
 * This specialization of slvr_piggy handles reading precomputed velocity
 * fields from a driver run and using them in a piggyback simulation.
 * It inherits from @ref output::hdf5_xdmf with
 * @ref solvers::mpdata_rhs_vip.
 *
 * @tparam ct_params_t Compile-time parameters.
 */
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
  std::string vel_in;
  blitz::TinyVector<hsize_t, parent_t::n_dims> read_shape_h, read_offst_h;

  
  protected:

  void read_vel()
  {
    if(this->rank==0)
    {
      using ix = typename ct_params_t::ix;
      H5::H5File h5f(vel_in+ "/velocities/" + this->hdf_name(this->base_name("velocity")), H5F_ACC_RDONLY
#if defined(USE_MPI)
        // set collective reading of velocity file, needs to be used together with dxpl_id. 
        // Significantly slows down for few MPI tasks, starts to be beneficial for many tasks
        // , hence disabled for now
//        , H5P_DEFAULT, this->fapl_id 
#endif
      );
      
      for (int d = 0; d < parent_t::n_dims; ++d)
      {
        H5::DataSet dataset = h5f.openDataSet(this->outvars[this->vip_ixs[d]].name);
        H5::DataSpace dataspace = dataset.getSpace();
        dataspace.selectHyperslab(H5S_SELECT_SET, read_shape_h.data(), read_offst_h.data());
        // in 3D, velocities are transposed to kji order before saving in record_aux_halo
        // hence we read in kji and copy back to kij used in 3D UWLCM
        // TODO: all this kij - kji copying could probably be avoided by using a smart
        //       hyperslab when reading (especially important with MPI)
        if(parent_t::n_dims == 3)
        {
          typename parent_t::arr_t kji_arr(read_shape_h);
          dataset.read(kji_arr.data(), this->flttype_solver, H5::DataSpace(parent_t::n_dims, read_shape_h.data()) , dataspace
#if defined(USE_MPI)
            // see above comments to h5f
  //          , this->dxpl_id
#endif
          );
          this->state(this->vip_ixs[d]) = kji_arr;
        }
        else
        {
          dataset.read(this->state(this->vip_ixs[d]).data(), this->flttype_solver, H5::DataSpace(parent_t::n_dims, read_shape_h.data()) , dataspace
#if defined(USE_MPI)
            // see above comments to h5f
  //          , this->dxpl_id
#endif
          );
        }
      }
    }
  }

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    for (int d = 0; d < parent_t::n_dims; ++d)
      read_shape_h[d] = this->mem->grid_size[d].length() + 2 * this->halo;

    read_offst_h = 0;
    read_offst_h[0] = this->mem->grid_size[0].first();

    if(this->rank==0)
    {
      po::options_description opts("Piggybacker options"); 
      opts.add_options()
        ("vel_in", po::value<std::string>()->required(), "directory with the 'velocities' direcotry (for piggybacking)")
      ;
      po::variables_map vm;
      handle_opts(opts, vm);
          
      vel_in = vm["vel_in"].as<std::string>();
      std::cout << "piggybacking from: " << vel_in << std::endl;

      this->record_aux_const("piggybacking", "piggy", "true");
      this->record_aux_const("vel_in", "piggy", vel_in); 
    }

    read_vel();
    this->mem->barrier();
  }

  void hook_post_step() 
  {
    parent_t::hook_post_step(); // do whatever
    this->mem->barrier(); //necessary?
    // read velo, overwrite any vel rhs
    read_vel();
    this->mem->barrier();
  }

  // ctor
  slvr_piggy(
    typename parent_t::ctor_args_t args,
    typename parent_t::rt_params_t const &p
  ) :
    parent_t(args, p) {}
};


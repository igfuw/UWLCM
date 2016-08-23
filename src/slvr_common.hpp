#pragma once

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>

using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t>
class slvr_common : public 
  output::hdf5_xdmf<
    solvers::mpdata_rhs_vip_prs<ct_params_t>
  >
{
  using parent_t = output::hdf5_xdmf<solvers::mpdata_rhs_vip_prs<ct_params_t>>;

  public:
  using real_t = typename ct_params_t::real_t;

  protected:

  real_t dx, dz; // 0->dx, 1->dy ! TODO
  int spinup; // number of timesteps

  // relaxation stuff
  bool relax_th_rv;
  blitz::Array<real_t, 2> &th_eq, &rv_eq;
  blitz::Array<real_t, 2> &th_ref, &rhod;

  // surface precip stuff
  real_t prec_vol;
  std::ofstream f_prec;
  
  // spinup stuff
  virtual bool get_rain() = 0;
  virtual void set_rain(bool) = 0;

  void hook_ante_loop(int nt) 
  {
    if (get_rain() == false) 
    {
      ;
      // spinup and relaxation do not make sense without autoconversion  (TODO: issue a warning?)
      // spinup = relax_th_rv = 0;      
    }
    if (spinup > 0)
    {
      blitz::secondIndex k;
      // initially the enviromental and reference profiles are the initial profiles
     /*
      rhod = setup::rhod_fctr()(k * dz);
      th_eq = setup::th_dry_fctr()(k * dz);
      th_ref = th_eq;
      rv_eq = setup::r_t()(k * dz);
*/
/*
      std::cout << "pre-spinup profiles:" << std::endl;
      std::cout << "rhod: " << rhod;
      std::cout << "th_ref: " << th_ref;
      std::cout << "th_eq: " << th_eq;
      std::cout << "rv_eq: " << rv_eq;
*/
      set_rain(false);
    }

    parent_t::hook_ante_loop(nt); 

    // open file for output of precitpitation volume
    f_prec.open(this->outdir+"/prec_vol.dat");
    prec_vol = 0.;
  }

  void hook_ante_step()
  {
    if (spinup != 0 && spinup == this->timestep)
    {
      // turn autoconversion on only after spinup (if spinup was specified)
      set_rain(true);
    }

    using ix = typename ct_params_t::ix;
    if(spinup == this->timestep)
    {
      // save horizontal means of th and rv after spinup
      // they will be the relaxation goals
      // and also the enviormental profiles for buoyancy
      // TODO: when calculating mean, do not include first or last point (which is the same in cyclic boundaries);
      //       right now it is accounted for twice, but the concurrency-aware sum cannot exclude single point
     /*
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {  
        th_eq(this->i, j) = this->mem->sum(this->state(ix::th), this->i, rng_t(j, j), false)  /  (this->mem->grid_size[0].length());
        rv_eq(this->i, j) = this->mem->sum(this->state(ix::rv), this->i, rng_t(j, j), false)  /  (this->mem->grid_size[0].length());
      }
*/
    // calculate reference theta and rhod
    // like in Wojtek's code
    {
/*
      blitz::secondIndex k;
      int nz = this->mem->grid_size[1].length();
      // calculate average stability
      blitz::Range notopbot(1, nz-2);
      blitz::Array<setup::real_t, 1> st(nz);
      st=0;
      st(notopbot) = (th_eq(0, notopbot+1) - th_eq(0, notopbot-1)) / th_eq(0, notopbot);
      setup::real_t st_avg = blitz::sum(st) / (nz-2) / (2.*dz);
      // reference theta
      th_ref = th_eq(0,0) * exp(st_avg * k * dz);
      // virtual temp at surface
      using libcloudphxx::common::moist_air::R_d_over_c_pd;
      using libcloudphxx::common::moist_air::c_pd;
      using libcloudphxx::common::moist_air::R_d;
      using libcloudphxx::common::theta_std::p_1000;

      setup::real_t T_surf = th_eq(0, 0) *  pow(setup::p_0 / p_1000<setup::real_t>(),  R_d_over_c_pd<setup::real_t>());
      setup::real_t T_virt_surf = T_surf * (1. + 0.608 * rv_eq(0, 0));
      setup::real_t rho_surf = (setup::p_0 / si::pascals) / T_virt_surf / 287. ; // TODO: R_d instead of 287
      setup::real_t cs = 9.81 / (c_pd<setup::real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / T_surf;
      // rhod profile
      rhod = rho_surf * exp(- st_avg * k * dz) * pow(
               1. - cs * (1 - exp(- st_avg * k * dz)), (1. / R_d_over_c_pd<setup::real_t>()) - 1);
*/
//      g_factor() = rhod;

/*
      std::cout << "post-spinup profiles:" << std::endl;
      std::cout << "rhod: " << rhod;
      std::cout << "th_ref: " << th_ref;
      std::cout << "th_eq: " << th_eq;
      std::cout << "rv_eq: " << rv_eq;
*/
    }

    }

    parent_t::hook_ante_step(); 
  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); 

    // recording total precipitation volume through the lower boundary
    if(this->rank==0)
    {
      f_prec << this->timestep << " "  << prec_vol << "\n";
      prec_vol = 0.;
    }
  }


  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  )   
  {   
    parent_t::update_rhs(rhs, dt, at);
    using ix = typename ct_params_t::ix;

    // relaxation terms; added only after spinup, when get_rain returns true
//    if(relax_th_rv && get_rain())
//    {
//      // computed level-wise
//      for (int j = this->j.first(); j <= this->j.last(); ++j)
//      {  
//        const auto tau = setup::tau_rlx / si::seconds * exp(j * dz / setup::z_rlx * si::metres);
//
//        for(auto a: std::list<int>({ix::th, ix::rv}))
//        {
//          const auto &psi = this->state(a);
//          // relax horizontal mean
//          /*
//          const auto psi_mean = this->mem->sum(psi, this->i, rng_t(j, j), false)  /  (this->mem->grid_size[0].length());
//          if(a == ix::th)
//            rhs.at(a)(this->i, j) =  (th_eq(j) - psi_mean) / tau;
//          else
//            rhs.at(a)(this->i, j) =  (rv_eq(j) - psi_mean) / tau;
//          */
//          // relax each cell 
//          if(a == ix::th)
//            rhs.at(a)(this->i, j) +=  (th_eq(j) - psi(this->i, j)) / tau;
//          else
//            rhs.at(a)(this->i, j) +=  (rv_eq(j) - psi(this->i, j)) / tau;
//        }
//      }
//    }
  }

  public:


  struct rt_params_t : parent_t::rt_params_t 
  { 
    typename ct_params_t::real_t dx = 0, dz = 0;
    int spinup = 0; // number of timesteps during which autoconversion is to be turned off
    bool relax_th_rv = true;
    typename parent_t::arr_t *th_e, *rv_e, *th_ref, *rhod;
    real_t z_rlx_sclr;
  };

  // ctor
  slvr_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    dx(p.dx),
    dz(p.dz),
    spinup(p.spinup),
    relax_th_rv(p.relax_th_rv),
    th_eq(*p.th_e),
    rv_eq(*p.rv_e),
    th_ref(*p.th_ref),
    rhod(*p.rhod)
  {
    assert(dx != 0);
    assert(dz != 0);
  }  
};

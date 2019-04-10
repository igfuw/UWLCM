#pragma once
#include <random>
#include <fstream>
#include "CasesCommon.hpp"
#include "LasherTrapp2001_sounding/x7221545.adjdec2.hpp"

namespace setup 
{
  namespace LasherTrapp
  {
    namespace moist_air = libcloudphxx::common::moist_air;
    namespace const_cp = libcloudphxx::common::const_cp;
    namespace theta_std = libcloudphxx::common::theta_std;

    using libcloudphxx::common::theta_std::p_1000;
    using libcloudphxx::common::moist_air::R_d_over_c_pd;
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::moist_air::R_d;
    using libcloudphxx::common::moist_air::R_v;
    using libcloudphxx::common::const_cp::l_tri;
    using libcloudphxx::common::const_cp::p_vs;
    using libcloudphxx::common::theta_std::p_1000;

    const quantity<si::pressure, real_t> 
      p_0 = 101800 * si::pascals;
    const quantity<si::length, real_t> 
      z_0  = 0    * si::metres,
      Z    = 8000 * si::metres, // DYCOMS: 1500
      X    = 10000 * si::metres, // DYCOMS: 6400
      Y    = 10000 * si::metres; // DYCOMS: 6400
    const real_t z_abs = 7000;
    const quantity<si::length, real_t> z_rlx_vctr = 25 * si::metres;

    // env profiles of th and rv from the sounding
    arr_1D_t th_dry_env;
    arr_1D_t th_std_env;
    arr_1D_t p_env;
    arr_1D_t rv_env;


    // RH T and p to rv
    quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
    }

    template<class case_ct_params_t, int n_dims>
    class LasherTrapp2001Common : public CasesCommon<case_ct_params_t, n_dims>
    {

      protected:
      using parent_t = CasesCommon<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
        params.outdir = user_params.outdir;
        params.outfreq = user_params.outfreq;
        params.spinup = user_params.spinup;
        params.w_src = user_params.w_src;
        params.uv_src = false; // ?
        params.th_src = user_params.th_src;
        params.rv_src = user_params.rv_src;
        params.rc_src = user_params.rc_src;
        params.rr_src = user_params.rr_src;
        params.dt = user_params.dt;
        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = false;
        params.friction = true;
        params.coriolis = false;
        params.radiation = false;

        this->setopts_sgs(params);
      }
  
  
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &solver, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        // we assume here that env_prof was called already, so that *_env profiles are initialized
        int nz = solver.advectee().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (Z / si::metres) / (nz-1); 
        // copy the env profiles into 2D/3D arrays
        solver.advectee(ix::rv) = rv_env(index); 
        solver.advectee(ix::th) = th_dry_env(index); 
  
        solver.advectee(ix::u) = 0;
        solver.advectee(ix::w) = 0;  
       
        // absorbers
        solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);
        solver.vab_relaxed_state(0) = solver.advectee(ix::u);
        solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
  
        // density profile
        solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
  
        // randomly prtrb tht and rv in the lowest 1km
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-0.1, 0.1);
          auto rand = std::bind(dis, gen);
  
          decltype(solver.advectee(ix::th)) prtrb(solver.advectee(ix::th).shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          prtrb = where(index * dz >= 1000., 0., prtrb); // no perturbation above 1km
          solver.advectee(ix::th) += prtrb;
        }
        {
          std::mt19937 gen(rng_seed);
          std::uniform_real_distribution<> dis(-0.025e-3, 0.025e-3);
          auto rand = std::bind(dis, gen);
  
          decltype(solver.advectee(ix::rv)) prtrb(solver.advectee(ix::rv).shape()); // array to store perturbation
          std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
          prtrb = where(index * dz >= 1000., 0., prtrb); // no perturbation above 1km
          solver.advectee(ix::rv) += prtrb;
        }
      }
  
  
      // calculate the initial environmental theta and rv profiles
      // alse set w_LS and hgt_fctrs
      void env_prof(profiles_t &profs, int nz, const user_params_t &user_params)
      {
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::const_cp::l_tri;
        using libcloudphxx::common::theta_std::p_1000;

        // read the soundings
        // containers for soundings
        std::vector<real_t> pres_s, temp_s, RH_s, z_s;
        for(std::string line : LasherTrapp2001_sounding_file)
        {
          real_t pres, temp, RH, z;
          sscanf(line.c_str(), "%*f %f %f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f", &pres, &temp, &RH, &z);
          pres_s.push_back(pres * 100); 
          temp_s.push_back(temp + 273.16);  // TODO: use libcloud's T_0 
          RH_s.push_back(RH / 100); 
          z_s.push_back(z); 
        }

        real_t dz = (Z / si::metres) / (nz-1); 

        // interpolate soundings to centers of cells 
        std::vector<real_t> pres_si(nz), temp_si(nz), RH_si(nz);
        real_t offset = z_s.at(0); // consider the lowest sounding level as ground level
        pres_si[0] = pres_s[0];
        temp_si[0] = temp_s[0];
        RH_si[0] = RH_s[0];
        int cell_no = 1;
        real_t z = cell_no * dz;
        for(int i=1; i<pres_s.size(); ++i)
        {
          real_t z_up = z_s.at(i) - offset;
          real_t z_down = z_s.at(i-1) - offset;
          while(z_down <= z && z < z_up)
          {
            real_t lin_fact = (z - z_down) / (z_up - z_down);
            pres_si[cell_no] = pres_s[i-1] + lin_fact * (pres_s[i] - pres_s[i-1]);
            temp_si[cell_no] = temp_s[i-1] + lin_fact * (temp_s[i] - temp_s[i-1]);
            RH_si[cell_no] = RH_s[i-1] + lin_fact * (RH_s[i] - RH_s[i-1]);
            ++cell_no;
            z = cell_no*dz;
            if(cell_no == nz) break;
          }
          if(cell_no == nz) break;
        }
        if(cell_no != nz) throw std::runtime_error("The initial sounding is not high enough");

        // calc derived profsiles
        std::vector<real_t> th_std(nz), th_dry(nz), rv(nz);
        for(int i=0; i<nz; ++i)
        {
          th_std[i] = temp_si[i] * pow(p_1000<real_t>() / si::pascals / pres_si[i], R_d<real_t>() / c_pd<real_t>());  
          rv[i] = RH_T_p_to_rv(RH_si[i], temp_si[i] * si::kelvins, pres_si[i] * si::pascals); 
          th_dry[i] = theta_dry::std2dry<real_t>(th_std[i] * si::kelvins, quantity<si::dimensionless, real_t>(rv[i])) / si::kelvins;
        }

        // create 1D blitz arrays to wrap the derived profsiles, store the for use in intcond_hlpr
        th_dry_env.resize(nz);
        th_std_env.resize(nz);
        p_env.resize(nz);
        rv_env.resize(nz);
        th_dry_env = arr_1D_t(th_dry.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        th_std_env = arr_1D_t(th_std.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        p_env = arr_1D_t(pres_si.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        rv_env     = arr_1D_t(rv.data(), blitz::shape(nz), blitz::neverDeleteData).copy();

        // TODO: calc hydrostatic env profsiles like in dycoms? w kodzie od S. L-T tego jednak nie ma...
        profs.p_e = p_env;
        profs.rv_e = rv_env;
        profs.rl_e = 0;
        profs.th_e = th_std_env; // temp to calc rhod
  
        // compute reference state theta and rhod
        blitz::firstIndex k;
        // calculate average stability
        blitz::Range notopbot(1, nz-2);
        arr_1D_t st(nz);
        st=0;
        st(notopbot) = (profs.th_e(notopbot+1) - profs.th_e(notopbot-1)) / profs.th_e(notopbot);
        real_t st_avg = blitz::sum(st) / (nz-2) / (2.*dz);
        // reference theta
        profs.th_ref = profs.th_e(0) * exp(st_avg * k * dz);
        // virtual temp at surface
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::theta_std::p_1000;
  
        real_t T_surf = profs.th_e(0) *  pow(p_0 / p_1000<real_t>(),  R_d_over_c_pd<real_t>());
        real_t T_virt_surf = T_surf * (1. + 0.608 * profs.rv_e(0));
        real_t rho_surf = (p_0 / si::pascals) / T_virt_surf / 287. ; // TODO: R_d instead of 287
        real_t cs = 9.81 / (c_pd<real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / T_surf;
        // rhod profsile
        profs.rhod = rho_surf * exp(- st_avg * k * dz) * pow(
                     1. - cs * (1 - exp(- st_avg * k * dz)), (1. / R_d_over_c_pd<real_t>()) - 1);

        profs.th_e = th_dry_env; // actual env profsile of theta_dry
  
        // calc divergence directly
        real_t z_0 = z_rlx_vctr / si::metres;
        profs.hgt_fctr_vctr = exp(- k * dz / z_0) / z_0;
        // for scalars
        z_0 = user_params.z_rlx_sclr;
        profs.hgt_fctr_sclr = exp(- k * dz / z_0) / z_0;
      }

      // functions that set surface fluxes per timestep
      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens, 
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy)
      {
        if(timestep == 0) 
          surf_flux_sens = .1; // [K * m/s]
        else if(int((3600. / dt) + 0.5) == timestep)
        {
          if(surf_flux_sens.rank() == 2) // TODO: make it a compile-time decision
            surf_flux_sens = .3 * exp( - ( pow(blitz::tensor::i * dx - 5000., 2) +  pow(blitz::tensor::j * dy - 5000., 2) ) / (1700. * 1700.) );
          else if(surf_flux_sens.rank() == 1)
            surf_flux_sens = .3 * exp( - ( pow(blitz::tensor::i * dx - 5000., 2)  ) / (1700. * 1700.) );
        }
      }
      
      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy)
      {
        if(timestep == 0)
          surf_flux_lat = .4e-4; // [1/s]
        else if(int((3600. / dt) + 0.5) == timestep)
        {
          if(surf_flux_lat.rank() == 2) // TODO: make it a compile-time decision
            surf_flux_lat = 1.2e-4 * exp( - ( pow(blitz::tensor::i * dx - 5000., 2) +  pow(blitz::tensor::j * dy - 5000., 2) ) / (1700. * 1700.) );
          else if(surf_flux_lat.rank() == 1)
            surf_flux_lat = 1.2e-4 * exp( - ( pow(blitz::tensor::i * dx - 5000., 2)  ) / (1700. * 1700.) );
        }
      }

      // ctor
      LasherTrapp2001Common()
      {
        //aerosol bimodal lognormal dist. - DYCOMS
        this->mean_rd1 = real_t(.011e-6) * si::metres,
        this->mean_rd2 = real_t(.06e-6) * si::metres;
        this->sdev_rd1 = real_t(1.2),
        this->sdev_rd2 = real_t(1.7);
        this->n1_stp = real_t(5*125e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(5*65e6) / si::cubic_metres;  // 65 || 16
        this->div_LS = real_t(0.);
        this->ForceParameters.surf_latent_flux_in_watts_per_square_meter = false; // it's given as mean(rv w) [kg/kg m/s]
        this->ForceParameters.surf_sensible_flux_in_watts_per_square_meter = false; // it's given as mean(theta) w [ K m/s]
        this->ForceParameters.u_fric = 0.28;
      }
    };
    
    template<class case_ct_params_t, int n_dims>
    class LasherTrapp2001;

    template<class case_ct_params_t>
    class LasherTrapp2001<case_ct_params_t, 2> : public LasherTrapp2001Common<case_ct_params_t, 2>
    {
      using parent_t = LasherTrapp2001Common<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      void intcond(typename parent_t::concurr_any_t &solver,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);
        this->make_cyclic(solver.advectee(ix::th));
      }
    };

    template<class case_ct_params_t>
    class LasherTrapp2001<case_ct_params_t, 3> : public LasherTrapp2001Common<case_ct_params_t, 3>
    {
      using parent_t = LasherTrapp2001Common<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Y / si::metres) / (nps[1]-1);
        params.dk = (Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      void intcond(typename parent_t::concurr_any_t &solver,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);
        this->make_cyclic(solver.advectee(ix::th));
  
        int nz = solver.advectee().extent(ix::w);
        real_t dz = (Z / si::metres) / (nz-1); 
  
        solver.advectee(ix::v)= 0;
        solver.vab_relaxed_state(1) = solver.advectee(ix::v);
      }
    };
  };
};

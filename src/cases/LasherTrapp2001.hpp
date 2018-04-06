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
    const quantity<si::length, real_t> z_rlx_vctr = 1 * si::metres;

    // env profiles of th and rv from the sounding
    arr_1D_t th_dry_env;
    arr_1D_t th_std_env;
    arr_1D_t rv_env;


    // RH T and p to rv
    quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
    }

    template<class concurr_t>
    class LasherTrapp2001 : public CasesCommon<concurr_t>
    {

      protected:
  
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
        params.dt = user_params.dt;
        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = false;
        params.friction = true;
        params.radiation = false;
      }
  
  
      template <class index_t>
      void intcond_hlpr(concurr_t &solver, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        // we assume here that env_prof was called already, so that *_env profiles are initialized
        using ix = typename concurr_t::solver_t::ix;
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
      void env_prof(arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &th_ref, arr_1D_t &pre_ref, arr_1D_t &rhod, arr_1D_t &w_LS, arr_1D_t &hgt_fctr_vctr, arr_1D_t &hgt_fctr_sclr, int nz, const user_params_t &user_params)
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

        using ix = typename concurr_t::solver_t::ix;
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

        // calc derived profiles
        std::vector<real_t> th_std(nz), th_dry(nz), rv(nz);
        for(int i=0; i<nz; ++i)
        {
          th_std[i] = temp_si[i] * pow(p_1000<real_t>() / si::pascals / pres_si[i], R_d<real_t>() / c_pd<real_t>());  
          rv[i] = RH_T_p_to_rv(RH_si[i], temp_si[i] * si::kelvins, pres_si[i] * si::pascals); 
          th_dry[i] = theta_dry::std2dry<real_t>(th_std[i] * si::kelvins, quantity<si::dimensionless, real_t>(rv[i])) / si::kelvins;
        }

        // create 1D blitz arrays to wrap the derived profiles, store the for use in intcond_hlpr
        th_dry_env.resize(nz);
        th_std_env.resize(nz);
        rv_env.resize(nz);
        th_dry_env = arr_1D_t(th_dry.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        th_std_env = arr_1D_t(th_std.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        rv_env     = arr_1D_t(rv.data(), blitz::shape(nz), blitz::neverDeleteData).copy();

        rv_e = rv_env;
        th_e = th_std_env; // temp to calc rhod
  
        // compute reference state theta and rhod
        blitz::firstIndex k;
        // calculate average stability
        blitz::Range notopbot(1, nz-2);
        arr_1D_t st(nz);
        st=0;
        st(notopbot) = (th_e(notopbot+1) - th_e(notopbot-1)) / th_e(notopbot);
        real_t st_avg = blitz::sum(st) / (nz-2) / (2.*dz);
        // reference theta
        th_ref = th_e(0) * exp(st_avg * k * dz);
        // virtual temp at surface
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::theta_std::p_1000;
  
        real_t T_surf = th_e(0) *  pow(p_0 / p_1000<real_t>(),  R_d_over_c_pd<real_t>());
        real_t T_virt_surf = T_surf * (1. + 0.608 * rv_e(0));
        real_t rho_surf = (p_0 / si::pascals) / T_virt_surf / 287. ; // TODO: R_d instead of 287
        real_t cs = 9.81 / (c_pd<real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / T_surf;
        // rhod profile
        rhod = rho_surf * exp(- st_avg * k * dz) * pow(
                 1. - cs * (1 - exp(- st_avg * k * dz)), (1. / R_d_over_c_pd<real_t>()) - 1);

        th_e = th_dry_env; // actual env profile of theta_dry
  
        // surface sources relaxation factors
        // for vectors
        real_t z_0 = z_rlx_vctr / si::metres;
        hgt_fctr_vctr = exp(- (k-0.5) * dz / z_0); // z=0 at k=1/2
        hgt_fctr_vctr(0) = 1;
        // for scalars
        z_0 = user_params.z_rlx_sclr;
        hgt_fctr_sclr = exp(- (k-0.5) * dz / z_0);
        hgt_fctr_sclr(0) = 1;
      }

      // ctor
      LasherTrapp2001()
      {
        //aerosol bimodal lognormal dist. - DYCOMS
        this->mean_rd1 = real_t(.011e-6) * si::metres,
        this->mean_rd2 = real_t(.06e-6) * si::metres;
        this->sdev_rd1 = real_t(1.2),
        this->sdev_rd2 = real_t(1.7);
        this->n1_stp = real_t(125e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(65e6) / si::cubic_metres;  // 65 || 16
        this->div_LS = real_t(0.);
        this->ForceParameters.surf_latent_flux_in_watts_per_square_meter = false; // it's given as a change in q_v [1/s]
        this->ForceParameters.u_fric = 0.28;
      }
    };

    template<class concurr_t>
    class LasherTrapp2001_2d : public LasherTrapp2001<concurr_t>
    {
      void setopts(typename concurr_t::solver_t::rt_params_t &params, int nx, int nz, const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nx-1); 
        params.dj = (Z / si::metres) / (nz-1);
        params.dz = params.dj;
      }

      void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);
        using ix = typename concurr_t::solver_t::ix;
        this->make_cyclic(solver.advectee(ix::th));
      }
    };

    template<class concurr_t>
    class LasherTrapp2001_3d : public LasherTrapp2001<concurr_t>
    {
      void setopts(typename concurr_t::solver_t::rt_params_t &params, int nx, int ny, int nz, const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nx-1); 
        params.dj = (Y / si::metres) / (ny-1);
        params.dk = (Z / si::metres) / (nz-1);
        params.dz = params.dk;
      }

      void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);
        using ix = typename concurr_t::solver_t::ix;
        this->make_cyclic(solver.advectee(ix::th));
  
        int nz = solver.advectee().extent(ix::w);
        real_t dz = (Z / si::metres) / (nz-1); 
  
        solver.advectee(ix::v)= 0;
        solver.vab_relaxed_state(1) = solver.advectee(ix::v);
      }
    };
  };
};

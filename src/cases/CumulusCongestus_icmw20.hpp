// Cumulus Congestus case from International Cloud Modeling Workshop 2020
// However, some changes were made (?) Actual setup as in DOI: 10.5194/gmd-17-759-2024 

#pragma once
#include <random>
#include <fstream>
#include "CumulusCongestusCommon.hpp"
#include "detail/LasherTrapp2001_sounding/x7221545.adjdec2.hpp"

namespace cases 
{
  namespace CumulusCongestus
  {
    template<class case_ct_params_t, int n_dims>
    class CumulusCongestusCommon_icmw20 : public CumulusCongestusCommon<case_ct_params_t, n_dims>
    {
      protected:
      using parent_t = CumulusCongestusCommon<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;
  
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &concurr, arr_1D_t &rhod, int rng_seed, index_t index) 
      {
        parent_t::intcond_hlpr(concurr, rhod, rng_seed, index, 0.1, 0.025e-3, (this->Z / si::metres) - 1000);
      }

      // calculate the initial environmental theta and rv profiles
      // alse set w_LS and hgt_fctrs
      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
      {
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::const_cp::l_tri;
        using libcloudphxx::common::theta_std::p_1000;

        parent_t::set_profs(profs, nz, user_params);

        // read the soundings
        // containers for soundings
        std::vector<real_t> pres_s, temp_s, RH_s, z_s;
        for(std::string line : LasherTrapp2001_sounding_file)
        {
          float pres, temp, RH, z;
          sscanf(line.c_str(), "%*f %f %f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f", &pres, &temp, &RH, &z);
          pres_s.push_back(pres * 100); 
          temp_s.push_back(temp + 273.16);  // TODO: use libcloud's T_0 
          RH_s.push_back(RH / 100); 
          z_s.push_back(z); 
        }

        real_t dz = (this->Z / si::metres) / (nz-1); 

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
        if(cell_no != nz) throw std::runtime_error("UWLCM: The initial sounding is not high enough");

        // calc derived profsiles
        std::vector<real_t> th_std(nz), th_dry(nz), rv(nz);
        for(int i=0; i<nz; ++i)
        {
          th_std[i] = temp_si[i] * pow(p_1000<real_t>() / si::pascals / pres_si[i], R_d<real_t>() / c_pd<real_t>());  
          rv[i] = parent_t::RH_T_p_to_rv(RH_si[i], temp_si[i] * si::kelvins, pres_si[i] * si::pascals); 
          th_dry[i] = theta_dry::std2dry<real_t>(th_std[i] * si::kelvins, quantity<si::dimensionless, real_t>(rv[i])) / si::kelvins;
        }

        // create 1D blitz arrays to wrap the derived profsiles, store the for use in intcond_hlpr
        this->th_dry_env.resize(nz);
        this->th_std_env.resize(nz);
        this->p_env.resize(nz);
        this->rv_env.resize(nz);
        this->th_dry_env = arr_1D_t(th_dry.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        this->th_std_env = arr_1D_t(th_std.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        this->p_env = arr_1D_t(pres_si.data(), blitz::shape(nz), blitz::neverDeleteData).copy();
        this->rv_env     = arr_1D_t(rv.data(), blitz::shape(nz), blitz::neverDeleteData).copy();

        // TODO: calc hydrostatic env profsiles like in dycoms? w kodzie od S. L-T tego jednak nie ma...
        profs.p_e =  this->p_env;
        profs.rv_e = this->rv_env;
        profs.rl_e = 0;
        profs.th_e = this->th_std_env; // temp to calc rhod
  
        // calc reference profiles
        this->ref_prof(profs, nz);

        profs.w_LS = 0.; // no subsidence
        profs.th_LS = 0.; // no large-scale horizontal advection
        profs.rv_LS = 0.; 
      }

      // functions that set surface fluxes per timestep
      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        parent_t::update_surf_flux_sens_hlpr(surf_flux_sens, th_ground, U_ground, U_ground_z, timestep, dt, dx, dy, 0.1, 1700, .3);
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_sens,
                                       blitz::Array<real_t, n_dims> th_ground,   
                                       blitz::Array<real_t, n_dims> U_ground,   
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        parent_t::update_surf_flux_lat_hlpr(surf_flux_sens, th_ground, U_ground, U_ground_z, timestep, dt, dx, dy, 4e-5, 1700, 1.2e-4);
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind
                               const real_t &U_ground_z,
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy, const real_t &uv_mean) override
      {
        surf_flux_uv = where(U_ground < 1e-4, 
            - 0.0784 * (uv_ground + uv_mean) / real_t(1e-4) * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters), // 0.0784 m^2 / s^2 is the square of friction velocity = 0.28 m / s
            - 0.0784 * (uv_ground + uv_mean) / U_ground * -1  * (this->rhod_0 / si::kilograms * si::cubic_meters)
          );
      }

      // ctor
      CumulusCongestusCommon_icmw20()
      {
        this->p_0 = real_t(101800) * si::pascals;
        //aerosol bimodal lognormal dist. - as in RICO with 11x conc following the ICMW2020 setup
        this->mean_rd1 = real_t(.03e-6) * si::metres,
        this->mean_rd2 = real_t(.14e-6) * si::metres;
        this->sdev_rd1 = real_t(1.28),
        this->sdev_rd2 = real_t(1.75);
        this->n1_stp = real_t(11*90e6) / si::cubic_metres, 
        this->n2_stp = real_t(11*15e6) / si::cubic_metres;
        this->z_rlx = real_t(1e2) * si::metres;
      }
    };
    
    template<class case_ct_params_t, int n_dims>
    class CumulusCongestus_icmw20;

    template<class case_ct_params_t>
    class CumulusCongestus_icmw20<case_ct_params_t, 2> : public CumulusCongestusCommon_icmw20<case_ct_params_t, 2>
    {
      using parent_t = CumulusCongestusCommon_icmw20<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);
      }

      public:
      CumulusCongestus_icmw20(const real_t _X=-1, const real_t _Y=-1, const real_t _Z=-1)
      {
        this->X = _X < 0 ? 12e3 * si::meters : _X * si::meters;
        this->Z = _Z < 0 ? 10e3 * si::meters : _Z * si::meters;
      }
    };

    template<class case_ct_params_t>
    class CumulusCongestus_icmw20<case_ct_params_t, 3> : public CumulusCongestusCommon_icmw20<case_ct_params_t, 3>
    {
      using parent_t = CumulusCongestusCommon_icmw20<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (this->X / si::metres) / (nps[0]-1); 
        params.dj = (this->Y / si::metres) / (nps[1]-1);
        params.dk = (this->Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      void intcond(typename parent_t::concurr_any_t &concurr,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(concurr, rhod, rng_seed, k);
  
        int nz = concurr.advectee_global().extent(ix::w);
        real_t dz = (this->Z / si::metres) / (nz-1); 
  
        concurr.advectee(ix::v)= 0;
        concurr.vab_relaxed_state(1) = concurr.advectee(ix::v);
      }

      public:
      CumulusCongestus_icmw20(const real_t _X=-1, const real_t _Y=-1, const real_t _Z=-1)
      {
        this->X = _X < 0 ? 10e3 * si::meters : _X * si::meters;
        this->Y = _Y < 0 ? 10e3 * si::meters : _Y * si::meters;
        this->Z = _Z < 0 ? 10e3 * si::meters : _Z * si::meters;
      }
    };
  };
};

  


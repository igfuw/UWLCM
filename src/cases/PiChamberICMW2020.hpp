// LES of Pi Chamber for the International Cloud Modeling Workshop 2020
// https://iccp2020.tropmet.res.in/Cloud-Modeling-Workshop-2020
// Setup based on Thomas et al. 2019, 
// but boundary ae conditions similar to the top and bottom conditions used in Grabowski 2019

#pragma once
#include "Anelastic.hpp"

namespace cases 
{
  namespace PiChamber2020
  {
    namespace const_cp = libcloudphxx::common::const_cp;

    // RH T and p to rv assuming RH = r_v / r_vs
    inline quantity<si::dimensionless, real_t> RH_T_p_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
    {
      return  RH * const_cp::r_vs<real_t>(T, p);
    }

    const quantity<si::length, real_t> 
     Z    ( 1 * si::metres), 
     X    ( 2 * si::metres), 
     Y    ( 2 * si::metres); 
    
    const quantity<si::pressure, real_t> p_0(100000 * si::pascals); // total pressure, const in the whole domain

    const quantity<si::dimensionless, real_t>
      RH_top(1),    // RH at top wall
      RH_bot(1),    // RH at bottom wall
      RH_side(0.82); // RH at side walls

    const real_t abs_dist = 0.03; // distance from walls in which velocity absorber is applied to mimick momentum flux...
    const quantity<si::time, real_t> abs_char_time = real_t(0.1) * si::seconds; // characteristic time defining maximal absorber strength. du/dt = 1 / abs_char_time * (u - u_relax)

    // initial temperature at height z
    inline quantity<si::temperature, real_t> T(const real_t &z)
    {
      return quantity<si::temperature, real_t>((299. - 19. * z / 1.) * si::kelvins);
    }

    // initial vapor at height z
    inline quantity<si::dimensionless, real_t> r_t(const real_t &z)
    {
//      return 0.0216 - 0.0154 * z / 1.; // default from the ICMW2020 case, gives different RH in UWLCM
      const real_t rv_top = RH_T_p_to_rv(RH_top, 280 * si::kelvins, p_0);
      const real_t rv_bot = RH_T_p_to_rv(RH_bot, 299 * si::kelvins, p_0);
      return rv_bot - (rv_bot - rv_top) * z / 1.;
    }

    template<class case_ct_params_t, int n_dims>
    class PiChamberICMW2020Common : public Anelastic<case_ct_params_t, n_dims>
    {

      protected:
      using parent_t = Anelastic<case_ct_params_t, n_dims>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // initial standard potential temp at height z, given constant pressure = 100000 Pa
      struct th_std_fctr
      {
        real_t operator()(const real_t &z) const
        {
          return T(z) / si::kelvins;
        }
        BZ_DECLARE_FUNCTOR(th_std_fctr);
      };

      // water mixing ratio at height z
      struct r_t_fctr
      {
        quantity<si::dimensionless, real_t> operator()(const real_t &z) const
        {
          return r_t(z);
        }
        BZ_DECLARE_FUNCTOR(r_t_fctr);
      };

      void setopts_hlpr(rt_params_t &params, const user_params_t &user_params)
      {
//        params.outdir = user_params.outdir;
//        params.outfreq = user_params.outfreq;
//        params.spinup = user_params.spinup;
//        params.w_src = user_params.w_src;
        // no explicit sources from boundaries (following Wojtek)
        // but what about transverse velocities at boundaries? add relaxation?
        params.uv_src = false;
        params.th_src = false;
        params.rv_src = false;
        params.rc_src = false;
        params.rr_src = false;
//        params.dt = user_params.dt;
//        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = false;
        params.vel_subsidence = false;
        params.friction = false;
        params.coriolis = false;
        params.radiation = false;
        params.no_ccn_at_init = true;
        params.open_side_walls = true;
    //    params.n_iters=1;

        this->setopts_sgs(params);
      }
    
      template <class index_t>
      void intcond_hlpr(typename parent_t::concurr_any_t &solver,
                        arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, int rng_seed, index_t index)
      {
        int nz = solver.advectee().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (Z / si::metres) / (nz-1); 
        int nx = solver.advectee().extent(0);
        real_t dx = (X / si::metres) / (nx-1); 
    
        solver.advectee(ix::u) = 0;
        solver.advectee(ix::w) = 0;  
    
        // density profile
        solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
    
        // initial potential temperature
        solver.advectee(ix::th) = th_std_fctr{}(index * dz);
    
        // initial water vapor mixing ratio
        solver.advectee(ix::rv) = r_t_fctr{}(index * dz);

        // --- absorbers ---

        // relaxed states
        solver.vab_relaxed_state(0) = 0;
        solver.vab_relaxed_state(ix::w) = 0;

        // coefficients. TODO: right now, values of coefficients in corners are not correct and depend on the order of operations below...
        const real_t max_vab_coeff = 1. / (abs_char_time / si::seconds);
        solver.vab_coefficient() = 0;
        // top wall
        solver.vab_coefficient() = where(index * dz >= Z / si::metres - abs_dist,  max_vab_coeff * pow(sin(3.1419 / 2. * (index * dz - (Z / si::metres - abs_dist))/ abs_dist), 2), solver.vab_coefficient());
        // bottom wall
        solver.vab_coefficient() = where(index * dz <= abs_dist,  max_vab_coeff * pow(sin(3.1419 / 2. * (abs_dist - index * dz )/ abs_dist), 2), solver.vab_coefficient());
        // left wall
        if(solver.vab_coefficient().lbound(0) == 0) // NOTE: with MPI, this condition would be true even for rank =0 at mpi rank > 0 (?)
          solver.vab_coefficient() = where(blitz::tensor::i * dx <= abs_dist,  max_vab_coeff * pow(sin(3.1419 / 2. * (abs_dist - blitz::tensor::i * dx )/ abs_dist), 2), solver.vab_coefficient());
        // right wall
        if(solver.vab_coefficient().ubound(0) == nx-1) // NOTE: with MPI, this condition would be true even for rank =0 at mpi rank > 0 (?)
          solver.vab_coefficient() = where(blitz::tensor::i * dx >= X / si::metres - abs_dist,  max_vab_coeff * pow(sin(3.1419 / 2. * (blitz::tensor::i * dx - (X / si::metres - abs_dist))/ abs_dist), 2), solver.vab_coefficient());
      }

    
      public:
      // calculate the initial environmental theta and rv profiles
      void set_profs(detail::profiles_t &profs, int nz, const user_params_t &user_params)
      // pre_ref - total pressure
      // th_e - dry potential temp
      // th_ref - dry potential temp refrence profsile
      // rhod - dry density profsile
      {

        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::const_cp::ls_tri;
        using libcloudphxx::common::theta_std::p_1000;
        using setup::real_t;

        parent_t::set_profs(profs, nz, user_params);
        // mix len not reduced neat bottom surface.
        // TODO: use smaller mix lean near surfaces (side and top too) ?
        profs.mix_len(0)=profs.mix_len(2);
        profs.mix_len(1)=profs.mix_len(2);

        real_t dz = (Z / si::metres) / (nz-1);
        blitz::firstIndex index;

        // env profiles
        profs.p_e = p_0 / si::pascals;
        profs.th_e = th_std_fctr{}(index * dz);
        profs.rv_e = r_t_fctr{}(index * dz);
        profs.rl_e = 0.;

        // turn supersaturation into water in the env profile
        // some constants first (from Wojtek's code)
        real_t tt0 = 273.17;
        real_t rv = 461; // specific gas constant for vapor
        real_t ee0 = 611.;
        real_t a = R_d<real_t>() / rv / si::joules * si::kelvins * si::kilograms; // aka epsilon
        real_t b = ls_tri<real_t>() / si::joules * si::kilograms / rv / tt0;
        real_t c = ls_tri<real_t>() / c_pd<real_t>() / si::kelvins;
        real_t d = ls_tri<real_t>() / si::joules * si::kilograms / rv;
        real_t f = R_d_over_c_pd<real_t>();

        for(int k=1; k<nz; ++k)
        {
          real_t thetme = 1; // pow(p_1000<real_t>() / si::pascals / profs.p_e(k), f); // 1/Exner
          real_t thi = 1. / (T(k * dz) / si::kelvins); // 1/theta_std, with theta_std = T
          real_t y = b * thetme * tt0 * thi;
          real_t ees = ee0 * exp(b-y); // saturation vapor pressure (Tetens equation or what?)
          real_t qvs = a * ees / (profs.p_e(k) - ees);  // saturation vapor mixing ratio = R_d / R_v * ees / p_d
          // calculate linearized condensation rate
          real_t cf1 = thetme*thetme*thi*thi;  // T^{-2}
          cf1 *= c * d * profs.p_e(k) / (profs.p_e(k) - ees); // = ls_tri^2 / (C_pd * R_v * T^2) * p/p_d
          real_t delta = (r_t(k*dz) - qvs) / (1 + qvs * cf1); // how much supersaturated is the air (divided by sth)
          if(delta < 0.) delta = 0.;
          profs.rv_e(k) = r_t(k*dz) - delta;
          profs.rl_e(k) = delta;
          profs.th_e(k) = T(k*dz) / si::kelvins + c * thetme * delta;
        }

        parent_t::ref_prof(profs, nz);

        // subsidence rate
        profs.w_LS = 0;
        profs.th_LS = 0.; // no large-scale horizontal advection
        profs.rv_LS = 0.;
      }


      void update_surf_flux_sens(blitz::Array<real_t, n_dims> surf_flux_sens,
                                 blitz::Array<real_t, n_dims> th_ground,    // value of th on the ground
                                 blitz::Array<real_t, n_dims> U_ground,     // magnitude of horizontal ground wind
                                 const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                                 const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
          surf_flux_sens = 0;
      }

      void update_surf_flux_lat(blitz::Array<real_t, n_dims> surf_flux_lat,
                                       blitz::Array<real_t, n_dims> rt_ground,
                                       blitz::Array<real_t, n_dims> U_ground,
                                       const real_t &U_ground_z,
                                       const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) override
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
          surf_flux_lat = 0;
      }

      // one function for updating u or v
      // the n_dims arrays have vertical extent of 1 - ground calculations only in here
      void update_surf_flux_uv(blitz::Array<real_t, n_dims>  surf_flux_uv, // output array
                               blitz::Array<real_t, n_dims>  uv_ground,    // value of u or v on the ground
                               blitz::Array<real_t, n_dims>  U_ground,     // magnitude of horizontal ground wind
                               const real_t &U_ground_z,                   // altituted at which U_ground is diagnosed
                               const int &timestep, const real_t &dt, const real_t &dx, const real_t &dy) 
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
          surf_flux_uv = 0;
      }

      // ctor
      PiChamberICMW2020Common()
      {
        this->p_0 = p_0;
       // this->kappa = 1.28; // NaCl aerosol
        this->Z = Z;
        this->n1_stp = real_t(0) / si::cubic_metres;
        this->n2_stp = real_t(0) / si::cubic_metres;
        this->ForceParameters.coriolis_parameter = 0.; 

      }
    };

    // 2d/3d children
    template<class case_ct_params_t, int n_dims>
    class PiChamberICMW2020;

    template<class case_ct_params_t>
    class PiChamberICMW2020<case_ct_params_t, 2> : public PiChamberICMW2020Common<case_ct_params_t, 2>
    {
      using parent_t = PiChamberICMW2020Common<case_ct_params_t, 2>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // function expecting a libmpdata solver parameters struct as argument
      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Z / si::metres) / (nps[1]-1);
        params.dz = params.dj;
      }

      // function expecting a libmpdata++ solver as argument
      void intcond(typename parent_t::concurr_any_t &solver,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(solver, rhod, th_e, rv_e, rl_e, rng_seed, k);
      }

      public:
      PiChamberICMW2020()
      {
        this->X = X;
      }
    };

    template<class case_ct_params_t>
    class PiChamberICMW2020<case_ct_params_t, 3> : public PiChamberICMW2020Common<case_ct_params_t, 3>
    {
      using parent_t = PiChamberICMW2020Common<case_ct_params_t, 3>;
      using ix = typename case_ct_params_t::ix;
      using rt_params_t = typename case_ct_params_t::rt_params_t;

      // function expecting a libmpdata solver parameters struct as argument
      void setopts(rt_params_t &params, const int nps[], const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nps[0]-1); 
        params.dj = (Y / si::metres) / (nps[1]-1);
        params.dk = (Z / si::metres) / (nps[2]-1);
        params.dz = params.dk;
      }

      // function expecting a libmpdata++ solver as argument
      void intcond(typename parent_t::concurr_any_t &solver,
                   arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &rl_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(solver, rhod, th_e, rv_e, rl_e, rng_seed, k);
    
        solver.advectee(ix::v) = 0;
        solver.vab_relaxed_state(1) = 0;

        int ny = solver.advectee().extent(1);
        real_t dy = (Y / si::metres) / (ny-1); 

        // absorption coefficients
        const real_t max_vab_coeff = 1. / (abs_char_time / si::seconds);
        // back wall
        solver.vab_coefficient() = where(blitz::tensor::j * dy >= Y / si::metres - abs_dist,  max_vab_coeff * pow(sin(3.1419 / 2. * (blitz::tensor::j * dy - (Y / si::metres - abs_dist))/ abs_dist), 2), solver.vab_coefficient());
        // front wall
        solver.vab_coefficient() = where(blitz::tensor::j * dy <= abs_dist,  max_vab_coeff * pow(sin(3.1419 / 2. * (abs_dist - blitz::tensor::j * dy )/ abs_dist), 2), solver.vab_coefficient());
      }

      public:
      PiChamberICMW2020()
      {
        this->X = X;
        this->Y = Y;
      }
    };
  };
};

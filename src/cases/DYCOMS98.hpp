#pragma once
#include <random>
#include <type_traits>
#include <iomanip>
#include "CasesCommon.hpp"
#include <libmpdata++/solvers/mpdata_rhs_vip_prs_sgs.hpp>

namespace setup 
{
  namespace dycoms
  {
    namespace hydrostatic = libcloudphxx::common::hydrostatic;
    namespace theta_std = libcloudphxx::common::theta_std;
    namespace theta_dry = libcloudphxx::common::theta_dry;
    namespace lognormal = libcloudphxx::common::lognormal;
  
    const quantity<si::pressure, real_t> 
      p_0 = 101780 * si::pascals;
    const quantity<si::length, real_t> 
      z_0  = 0    * si::metres,
      Z    = 1500 * si::metres, // DYCOMS: 1500
      X    = 6400 * si::metres, // DYCOMS: 6400
      Y    = 6400 * si::metres; // DYCOMS: 6400
    const real_t z_abs = 1250;
    const real_t z_i = 795; //initial inversion height
    const quantity<si::length, real_t> z_rlx_vctr = 25 * si::metres;
  
    // liquid water potential temperature at height z
    quantity<si::temperature, real_t> th_l(const real_t &z)
    {
      quantity<si::temperature, real_t> ret;
      ret = z < z_i ?
        288.3 * si::kelvins : 
        (295. + pow(z - z_i, 1./3)) * si::kelvins;
      return ret;
    }
  
    // water mixing ratio at height z
    struct r_t
    {
      quantity<si::dimensionless, real_t> operator()(const real_t &z) const
      {
        const quantity<si::dimensionless, real_t> q_t = z < z_i ?
          9.45e-3 : 
          (5. - 3. * (1. - exp((z_i - z)/500.))) * 1e-3;
        return q_t;
      }
      BZ_DECLARE_FUNCTOR(r_t);
    };
  
    // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
    struct th_dry_fctr
    {
      real_t operator()(const real_t &z) const
      {
        return theta_dry::std2dry<real_t>(th_l(z), r_t()(z)) / si::kelvins;
      }
      BZ_DECLARE_FUNCTOR(th_dry_fctr);
    };
  
    // westerly wind
    struct u
    {
      real_t operator()(const real_t &z) const
      {
        return 3. + 4.3 * z / 1000.; 
      }
      BZ_DECLARE_FUNCTOR(u);
    };
  
    // southerly wind
    struct v
    {
      real_t operator()(const real_t &z) const
      {
        return -9. + 5.6 * z / 1000.; 
      }
      BZ_DECLARE_FUNCTOR(v);
    };
  
    // large-scale vertical wind
    struct w_LS_fctr
    {
      real_t operator()(const real_t &z) const
      {
        return - D * z; 
      }
      BZ_DECLARE_FUNCTOR(w_LS_fctr);
    };
  
    // density profile as a function of altitude
    // hydrostatic and assuming constant theta (not used now)
    struct rhod_fctr
    {
      real_t operator()(real_t z) const
      {
        quantity<si::pressure, real_t> p = hydrostatic::p(
  	z * si::metres, th_dry_fctr()(0.) * si::kelvins, r_t()(0.), z_0, p_0
        );
        
        quantity<si::mass_density, real_t> rhod = theta_std::rhod(
  	p, th_dry_fctr()(0.) * si::kelvins, r_t()(0.)
        );
  
        return rhod / si::kilograms * si::cubic_metres;
      }
  
      // to make the rhod() functor accept Blitz arrays as arguments
      BZ_DECLARE_FUNCTOR(rhod_fctr);
    };

    template<class concurr_t>
    class Dycoms98 : public CasesCommon<concurr_t>
    {
      protected:
      
      template<class ct_params_t = typename concurr_t::solver_t::ct_params_t_>
      void setopts_sgs(typename concurr_t::solver_t::rt_params_t &params,
                       typename std::enable_if<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles>::type* = 0)
      {}

      template<class ct_params_t = typename concurr_t::solver_t::ct_params_t_>
      void setopts_sgs(typename concurr_t::solver_t::rt_params_t &params,
                       typename std::enable_if<ct_params_t::sgs_scheme == libmpdataxx::solvers::smg>::type* = 0)
      {
        params.c_m = 0.0856;
        params.smg_c = 0.165;
        params.prandtl_num = 0.42;
        params.cdrag = 0.;
      }
  
      template <class T, class U>
      void setopts_hlpr(T &params, const U &user_params)
      {
        params.outdir = user_params.outdir;
        params.outfreq = user_params.outfreq;
        params.spinup = user_params.spinup;
        params.w_src = user_params.w_src;
        params.uv_src = user_params.uv_src;
        params.th_src = user_params.th_src;
        params.rv_src = user_params.rv_src;
        params.rc_src = user_params.rc_src;
        params.rr_src = user_params.rr_src;
        params.dt = user_params.dt;
        params.nt = user_params.nt;
        params.buoyancy_wet = true;
        params.subsidence = true;
        params.friction = true;
        params.radiation = true;

        setopts_sgs(params);
      }
  
      template <class index_t>
      void intcond_hlpr(concurr_t &solver, arr_1D_t &rhod, int rng_seed, index_t index)
      {
        using ix = typename concurr_t::solver_t::ix;
        int nz = solver.advectee().extent(ix::w);  // ix::w is the index of vertical domension both in 2D and 3D
        real_t dz = (Z / si::metres) / (nz-1); 
  
        solver.advectee(ix::rv) = r_t()(index * dz); 
        solver.advectee(ix::u)= u()(index * dz);
        solver.advectee(ix::w) = 0;  
       
        // absorbers
        solver.vab_coefficient() = where(index * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (index * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);
        solver.vab_relaxed_state(0) = solver.advectee(ix::u);
        solver.vab_relaxed_state(ix::w) = 0; // vertical relaxed state
  
        // density profile
        solver.g_factor() = rhod(index); // copy the 1D profile into 2D/3D array
  
        // initial potential temperature
        solver.advectee(ix::th) = th_dry_fctr()(index * dz); 
        // randomly prtrb tht
        std::mt19937 gen(rng_seed);
        std::uniform_real_distribution<> dis(-0.1, 0.1);
        auto rand = std::bind(dis, gen);
  
        decltype(solver.advectee(ix::th)) prtrb(solver.advectee(ix::th).shape()); // array to store perturbation
        std::generate(prtrb.begin(), prtrb.end(), rand); // fill it, TODO: is it officialy stl compatible?
        solver.advectee(ix::th) += prtrb;
      }
  
      // calculate the initial environmental theta and rv profiles
      // like in Wojtek's BabyEulag
      // alse set w_LS and hgt_fctrs
      // TODO: move hgt_fctrs from cases to main code
      void env_prof(arr_1D_t &th_e, arr_1D_t &p_e, arr_1D_t &rv_e, arr_1D_t &th_ref, arr_1D_t &pre_ref, arr_1D_t &rhod, arr_1D_t &w_LS, arr_1D_t &hgt_fctr_vctr, arr_1D_t &hgt_fctr_sclr, arr_1D_t& mix_len, int nz, const user_params_t &user_params)
      {
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::const_cp::l_tri;
        using libcloudphxx::common::theta_std::p_1000;
  
        // temp profile
        arr_1D_t T(nz);
        real_t dz = (Z / si::metres) / (nz-1);
  
        r_t rt;
        // input sounding at z=0, for moist air, no liquid water
        T(0) = th_l(0.) / si::kelvins *  pow(p_0 / p_1000<real_t>(),  R_d_over_c_pd<real_t>());
        p_e(0) = p_0 / si::pascals;
        th_e(0) = th_l(0.) / si::kelvins;
        rv_e(0) = rt(0.);
  
        real_t tt0 = 273.17;
        real_t rv = 461; // specific gas constant for vapor
        real_t ee0 = 611.;
        real_t a = R_d<real_t>() / rv / si::joules * si::kelvins * si::kilograms; // aka epsilon
        real_t b = l_tri<real_t>() / si::joules * si::kilograms / rv / tt0;
        real_t c = l_tri<real_t>() / c_pd<real_t>() / si::kelvins;
        real_t d = l_tri<real_t>() / si::joules * si::kilograms / rv;
        real_t f = R_d_over_c_pd<real_t>(); 

real_t lwp_env = 0;
  
        for(int k=1; k<nz; ++k)
        {
          real_t bottom = R_d<real_t>() / si::joules * si::kelvins * si::kilograms * T(k-1) * (1 + 0.61 * rv_e(k-1)); // (p / rho) of moist air at k-1
          real_t rho1 = p_e(k-1) / bottom; // rho at k-1
          p_e(k) = p_e(k-1) - rho1 * 9.81 * dz; // estimate of pre at k (dp = -g * rho * dz)
          real_t thetme = pow(p_1000<real_t>() / si::pascals / p_e(k), f); // 1/Exner
          real_t thi = 1. / (th_l(k * dz) / si::kelvins); // 1/theta_std
          real_t y = b * thetme * tt0 * thi; 
          real_t ees = ee0 * exp(b-y); // saturation vapor pressure (Tetens equation or what?)
          real_t qvs = a * ees / (p_e(k) - ees);  // saturation vapor mixing ratio = R_d / R_v * ees / p_d
// calculate linearized condensation rate
          real_t cf1 = thetme*thetme*thi*thi;  // T^{-2}
          cf1 *= c * d * p_e(k) / (p_e(k) - ees); // = l_tri^2 / (C_pd * R_v * T^2) * p/p_d
          real_t delta = (rt(k*dz) - qvs) / (1 + qvs * cf1); // how much supersaturated is the air (divided by sth)
          if(delta < 0.) delta = 0.;
          rv_e(k) = rt(k*dz) - delta;
          th_e(k) = th_l(k*dz) / si::kelvins + c * thetme * delta;
          T(k) = th_e(k) * pow(p_e(k) / (p_1000<real_t>() / si::pascals),  f);

          bottom = R_d<real_t>() / si::joules * si::kelvins * si::kilograms * T(k) * (1 + 0.61 * rv_e(k)); // (p / rho) of moist air at k-1
          rho1 = p_e(k) / bottom; // rho at k-1
          lwp_env  += delta * rho1;
std::cout << k << " env_prof temp: " << T(k)  << " env prof delta: " << delta << std::endl;
        }
        lwp_env = lwp_env * 5  * 1e3;
std::cout << "lwp env: " << lwp_env << std::endl;
  
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
      //  th_ref = th_e(0) * pow(1 + rv_e(0) / a, f) // calc dry theta at z=0 
      //           * exp(st_avg * k * dz);
        // virtual temp at surface
        using libcloudphxx::common::moist_air::R_d_over_c_pd;
        using libcloudphxx::common::moist_air::c_pd;
        using libcloudphxx::common::moist_air::R_d;
        using libcloudphxx::common::theta_std::p_1000;
  
        real_t T_surf = th_e(0) *  pow(p_0 / p_1000<real_t>(),  R_d_over_c_pd<real_t>());

        real_t T_virt_surf = T_surf * (1. + 0.608 * rv_e(0));
        real_t rho_surf = (p_0 / si::pascals) / T_virt_surf / 287. ; // TODO: R_d instead of 287, its the total, not dry density!
//        rho_surf /= (1 + rv_e(0)); // turn it into dry air density! TODO: is this correct? TODO2: approp change in the paper

     //   real_t rho_surf = (p_0 / si::pascals) / T_surf / (1. + 29. / 18. * rv_e(0)) / 287. ; // dry air density at the surface TODO: R_d instead of 287

        // real_t cs = 9.81 / (c_pd<real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / T_surf; // this is from Wojtek
         real_t cs = 9.81 / (c_pd<real_t>() / si::joules * si::kilograms * si::kelvins) / st_avg / th_e(0);  // this is correct? or total, not dry th_e(0) should be here?
        // rhod profile
        rhod = rho_surf * exp(- st_avg * k * dz) * pow(
                 1. - cs * (1 - exp(- st_avg * k * dz)), (1. / R_d_over_c_pd<real_t>()) - 1);


        // theta_std env prof to theta_dry_e
        for(int k=1; k<nz; ++k)
          th_e(k) = theta_dry::std2dry<real_t>(th_e(k) * si::kelvins, quantity<si::dimensionless, real_t>(rv_e(k))) / si::kelvins;
  
        // subsidence rate
        w_LS = w_LS_fctr()(k * dz);
  
        // surface sources relaxation factors
        // for vectors
        /*
        real_t z_0 = z_rlx_vctr / si::metres;
        hgt_fctr_vctr = exp(- (k-0.5) * dz / z_0); // z=0 at k=1/2
        hgt_fctr_vctr(0) = 1;
        // for scalars
        z_0 = user_params.z_rlx_sclr;
        hgt_fctr_sclr = exp(- (k-0.5) * dz / z_0);
        hgt_fctr_sclr(0) = 1;
        */
 
        // calc divergence directly
        real_t z_0 = z_rlx_vctr / si::metres;
        hgt_fctr_vctr = exp(- k * dz / z_0) / z_0;
        // for scalars
        z_0 = user_params.z_rlx_sclr;
        hgt_fctr_sclr = exp(- k * dz / z_0) / z_0;

        auto delta = dz;
        mix_len = min(max(k, 1) * dz * 0.845, delta);
      }

      void update_surf_flux_sens(typename concurr_t::solver_t::arr_t &surf_flux_sens, int timestep, real_t dt)
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
          surf_flux_sens = 16.; // [W/m^2]
      }

      void update_surf_flux_lat(typename concurr_t::solver_t::arr_t &surf_flux_lat, int timestep, real_t dt)
      {
        if(timestep == 0) // TODO: what if this function is not called at t=0? force such call
          surf_flux_lat = 93.; // [W/m^2]
      }

      // ctor
      Dycoms98()
      {
        //aerosol bimodal lognormal dist. - DYCOMS
        this->mean_rd1 = real_t(.011e-6) * si::metres,
        this->mean_rd2 = real_t(.06e-6) * si::metres;
        this->sdev_rd1 = real_t(1.2),
        this->sdev_rd2 = real_t(1.7);
        this->n1_stp = real_t(125e6) / si::cubic_metres, // 125 || 31
        this->n2_stp = real_t(65e6) / si::cubic_metres;  // 65 || 16
        this->div_LS = real_t(3.75e-6); // [1/s] large-scale wind divergence used to calc subsidence of SDs, TODO: use boost.units to enforce 1/s
      }
    };

    template<class concurr_t>
    class Dycoms98_2d : public Dycoms98<concurr_t>
    {
      void setopts(typename concurr_t::solver_t::rt_params_t &params, int nx, int nz, const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        params.di = (X / si::metres) / (nx-1); 
        params.dj = (Z / si::metres) / (nz-1);
        params.dz = params.dj;
      }

      void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::secondIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);
        using ix = typename concurr_t::solver_t::ix;
        this->make_cyclic(solver.advectee(ix::th));
      }
    };

    template<class concurr_t>
    class Dycoms98_3d : public Dycoms98<concurr_t>
    {
      void setopts(typename concurr_t::solver_t::rt_params_t &params, int nx, int ny, int nz, const user_params_t &user_params)
      {
        this->setopts_hlpr(params, user_params);
        this->setopts_sgs(params);
        params.di = (X / si::metres) / (nx-1); 
        params.dj = (Y / si::metres) / (ny-1);
        params.dk = (Z / si::metres) / (nz-1);
        params.dz = params.dk;
      }

      void intcond(concurr_t &solver, arr_1D_t &rhod, arr_1D_t &th_e, arr_1D_t &rv_e, arr_1D_t &p_e, int rng_seed)
      {
        blitz::thirdIndex k;
        this->intcond_hlpr(solver, rhod, rng_seed, k);
        using ix = typename concurr_t::solver_t::ix;
        this->make_cyclic(solver.advectee(ix::th));
  
        int nz = solver.advectee().extent(ix::w);
        real_t dz = (Z / si::metres) / (nz-1); 
  
        solver.advectee(ix::v)= v()(k * dz);
        solver.vab_relaxed_state(1) = solver.advectee(ix::v);
      }
    };
  };
};

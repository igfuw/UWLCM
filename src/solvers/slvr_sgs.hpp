#pragma once

#include "slvr_common.hpp"
#include "../detail/blitz_hlpr_fctrs.hpp"
#include "../formulae/stress_formulae.hpp"

// common isotropic/anisotropic Smagorinsky stuff
template <class ct_params_t>
class slvr_sgs_smg_common : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;

  public:
  using real_t = typename ct_params_t::real_t;
  using ix = typename ct_params_t::ix;

  protected:

  real_t prandtl_num, karman_c; // karman_c used only in aniso...

  typename parent_t::arr_t &rcdsn_num, &tdef_sq, &tke, &sgs_th_flux, &sgs_rv_flux;
  arrvec_t<typename parent_t::arr_t> &tmp_grad, &sgs_momenta_fluxes;
  
  void calc_rcdsn_num()
  {
    using libmpdataxx::arakawa_c::h;

    const auto g = (libcloudphxx::common::earth::g<setup::real_t>() / si::metres_per_second_squared);

    const auto dz = this->params.dz;
    const auto& tht = this->state(ix::th);
    const auto& rv = this->state(ix::rv);
    // depending on microphysics we either have rc already (blk_m1) or have to diagnose it (lgrngn)
    const auto& rc = this->get_rc(rcdsn_num); // use rcdsn_num as temp storage for rc
   
    // libcloudph stuff
    const auto l_tri = libcloudphxx::common::const_cp::l_tri<setup::real_t>() * si::kilograms / si::joules;
    const auto eps = libcloudphxx::common::moist_air::eps<setup::real_t>();// / si::dimensionless;
    const auto c_pd = libcloudphxx::common::moist_air::c_pd<setup::real_t>() * si::kilograms * si::kelvins / si::joules;
    const auto R_d = libcloudphxx::common::moist_air::R_d<setup::real_t>() * si::kilograms  * si::kelvins/ si::joules;
    using libcloudphxx::common::theta_std::exner;

    // some constant coefficients
    const auto cf1 = (1 - eps) / eps;
    const auto cf2 = l_tri / R_d;
    const auto cf3 = l_tri / c_pd;
    const auto cf4 = eps * cf2 * cf3;

    // TODO: loops are bad, very bad !
    for (int k = this->vert_rng.first(); k <= this->vert_rng.last() - 1; ++k)
    {
      const auto th_ref_kph = 0.5 * ((*this->params.th_ref)(k + 1) + (*this->params.th_ref)(k));
      const auto dthtdz_kph = (this->hrzntl_slice(tht, k + 1) - this->hrzntl_slice(tht, k)) / dz;
      const auto rv_kph = 0.5 * (this->hrzntl_slice(rv, k + 1) + this->hrzntl_slice(rv, k));
      const auto drvdz_kph = (this->hrzntl_slice(rv, k + 1) - this->hrzntl_slice(rv, k)) / dz;
      
      const auto N2unsat = g * (dthtdz_kph / th_ref_kph + cf1 / (1 + cf1 * rv_kph) * drvdz_kph);
     
      const auto T_kp1 = this->hrzntl_slice(tht, k + 1) * exner((*this->params.p_e)(k + 1) * si::pascals);
      const auto T_k = this->hrzntl_slice(tht, k) * exner((*this->params.p_e)(k) * si::pascals);
      const auto T_kph = 0.5 * (T_kp1 + T_k);
      const auto drwdz_kph = ( this->hrzntl_slice(rv, k + 1) + this->hrzntl_slice(rc, k + 1) 
                             - this->hrzntl_slice(rv, k) - this->hrzntl_slice(rc, k)         ) / dz;

      const auto gamma = (1 + cf2 * rv_kph / T_kph) / (1 + cf4 * rv_kph / (T_kph * T_kph));
      
      const auto N2sat = g * (gamma * (dthtdz_kph / th_ref_kph + cf3 * drvdz_kph / T_kph) - drwdz_kph);
      
      const auto rc_kph = 0.5 * (this->hrzntl_slice(rc, k + 1) + this->hrzntl_slice(rc, k));

      tmp_grad[ct_params_t::n_dims - 1](this->hrzntl_slice(k + h))
      =
      blitz::where(rc_kph > 1e-6, N2sat, N2unsat);
    }

    // boundary conditions
    tmp_grad[ct_params_t::n_dims - 1](this->hrzntl_slice(0 - h)) = tmp_grad[ct_params_t::n_dims - 1](this->hrzntl_slice(0 + h));
    auto lk = this->vert_rng.last();
    tmp_grad[ct_params_t::n_dims - 1](this->hrzntl_slice(lk + h)) = tmp_grad[ct_params_t::n_dims - 1](this->hrzntl_slice(lk - h));
    
    this->vert_aver_cmpct(tmp_grad[ct_params_t::n_dims - 1], rcdsn_num);
    rcdsn_num(this->ijk) /= max(1e-15, tdef_sq(this->ijk)); // TODO: is 1e-15 sensible epsilon here ?
  }
  
  template <int nd = ct_params_t::n_dims> 
  void calc_sgs_momenta_fluxes(typename std::enable_if<nd == 2>::type* = 0)
  {
    this->sgs_momenta_fluxes[0](this->ijk) = ( this->tau[2](this->i - h, this->j - h)
                               + this->tau[2](this->i + h, this->j - h)
                               + this->tau[2](this->i + h, this->j + h)
                               + this->tau[2](this->i - h, this->j + h)
                               ) / 4;
  }
  
  template <int nd = ct_params_t::n_dims> 
  void calc_sgs_momenta_fluxes(typename std::enable_if<nd == 3>::type* = 0)
  {
    this->sgs_momenta_fluxes[0](this->ijk) = ( this->tau[4](this->i - h, this->j, this->k - h)
                               + this->tau[4](this->i + h, this->j, this->k - h)
                               + this->tau[4](this->i + h, this->j, this->k + h)
                               + this->tau[4](this->i - h, this->j, this->k + h)
                               ) / 4;
    this->sgs_momenta_fluxes[1](this->ijk) = ( this->tau[5](this->i, this->j - h, this->k - h)
                               + this->tau[5](this->i, this->j + h, this->k - h)
                               + this->tau[5](this->i, this->j + h, this->k + h)
                               + this->tau[5](this->i, this->j - h, this->k + h)
                               ) / 4;
  }

  void calc_sgs_flux(int s)
  {
    if (s != ix::th && s != ix::rv) return;

    real_t conv_fctr = 1.;

    if (s == ix::th)
    {
      auto conv_fctr_sens = (libcloudphxx::common::moist_air::c_pd<real_t>() * si::kilograms * si::kelvins / si::joules);
      conv_fctr = conv_fctr_sens;
      this->vert_aver_cmpct(tmp_grad[ct_params_t::n_dims - 1], sgs_th_flux, conv_fctr);
    }
    else if (s == ix::rv)// || s == ix::rc)
    {
      auto conv_fctr_lat = (libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules);
      conv_fctr = conv_fctr_lat;
      this->vert_aver_cmpct(tmp_grad[ct_params_t::n_dims - 1], sgs_rv_flux, conv_fctr);
    }
  }
  
  void calc_drag_cmpct() override
  {
    if(params.cdrag > 0)           // kinematic momentum flux = - cdrag |U| u
    {
      formulae::stress::calc_drag_cmpct<ct_params_t::n_dims, ct_params_t::opts>(this->tau_srfc,
                                                                            this->vips(),
                                                                            *this->mem->G,
                                                                            this->cdrag,
                                                                            this->ijk,
                                                                            this->ijkm,
                                                                            this->params.ForceParameters.uv_mean);
    }
    else if(params.fricvelsq > 0)  // kinematic momentum flux = - fricvelsq u / |U|
    {
      // have to use modified ijkm due to shared-memory parallelisation, otherwise overlapping ranges
      // would lead to trobule as tau_srfc in calc_drag_cmpct_fricvel is used both as tmp storage and output 
      formulae::stress::calc_drag_cmpct_fricvel<ct_params_t::n_dims, ct_params_t::opts>(this->tau_srfc,
                                                                                this->vips(),
                                                                                *this->mem->G,
                                                                                params.fricvelsq,
                                                                                this->ijk,
                                                                                this->ijkm_sep,
                                                                                this->params.ForceParameters.uv_mean);
    }
  }
  
  void sgs_scalar_forces(const std::vector<int> &sclr_indices) override
  {
    for (const auto s : sclr_indices)
    {
      auto& field = this->state(s);

      this->xchng_pres(field, this->ijk);

      formulae::nabla::calc_grad_cmpct<parent_t::n_dims>(tmp_grad, field, this->ijk, this->ijkm, this->dijk);
      for(int d = 0; d < parent_t::n_dims; ++d)
        nancheck(tmp_grad[d](this->ijk), "tmp_grad in sgs_scalar_forces after calc_grad_cmpct");

      // document why
      this->mem->barrier();

      // ijk_vec is used, because MPI requires that thread rank 0 calculates next vector to the left of the process' domain
      formulae::stress::multiply_vctr_cmpct<ct_params_t::n_dims, ct_params_t::opts>(tmp_grad,
                                                                                    1.0 / prandtl_num,
                                                                                    this->k_m,
                                                                                    *this->mem->G,
                                                                                    this->ijk_vec);
      for(int d = 0; d < parent_t::n_dims; ++d)
        nancheck(tmp_grad[d](this->ijk), "tmp_grad in sgs_scalar_forces after multiply_vctr_cmpct");

      if (s == ix::th)
      {
        this->xchng_sgs_vctr(tmp_grad, this->surf_flux_sens, this->ijk);
      }
      else if (s == ix::rv)
      {
        this->xchng_sgs_vctr(tmp_grad , this->surf_flux_lat , this->ijk);
      }
      else
      {
        this->xchng_sgs_vctr(tmp_grad , this->surf_flux_zero, this->ijk);
      }

      for(int d = 0; d < parent_t::n_dims; ++d)
        nancheck(tmp_grad[d](this->ijk), "tmp_grad in sgs_scalar_forces after xchng_sgs_vctr");

      if (this->timestep == 0 || ((this->timestep + 1) % static_cast<int>(this->outfreq) < this->outwindow)) // timstep is increased after ante_step, i.e after update_rhs(at=0) that calls sgs_scalar_forces
        calc_sgs_flux(s);
    
      if (s == ix::th) // dth/dt = dT/dt / exner
      {
        this->tmp1(this->ijk) = formulae::stress::flux_div_cmpct<parent_t::n_dims, ct_params_t::opts>(
                                      tmp_grad,
                                      *this->mem->G,
                                      this->ijk,
                                      this->dijk
                                    );
        this->tmp1(this->ijk).reindex(this->zero) /= calc_exner()((*params.p_e)(this->vert_idx));
        this->rhs.at(s)(this->ijk) += this->tmp1(this->ijk); 
      }
      else
      {
        this->rhs.at(s)(this->ijk) += formulae::stress::flux_div_cmpct<parent_t::n_dims, ct_params_t::opts>(
                                      tmp_grad,
                                      *this->mem->G,
                                      this->ijk,
                                      this->dijk
                                    );
      }
    }
  }

  virtual void calc_sgs_diag_fields() = 0;

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at);

    // explicit application of subgrid forcings
    if(at == 0)
    {
      if (this->timestep == 0 || ((this->timestep + 1) % static_cast<int>(this->outfreq) < this->outwindow)) // timstep is increased after ante_step, i.e after update_rhs(at=0)
        this->calc_sgs_diag_fields();

      sgs_scalar_forces({ix::th, ix::rv});
      nancheck(rhs.at(ix::th)(this->ijk), "RHS of th after sgs_scalar_forces");
      nancheck(rhs.at(ix::rv)(this->ijk), "RHS of rv after sgs_scalar_forces");
    }
  }

  virtual void record_k_m() = 0;

  void diag() override
  {
    assert(this->rank == 0);
    parent_t::diag();
//    std::cout << "test u: " << min(this->state(ix::u)(this->domain)) << ' ' << max(this->state(ix::u)(this->domain)) << std::endl;
//    std::cout << "test w: " << min(this->state(ix::w)(this->domain)) << ' ' << max(this->state(ix::w)(this->domain)) << std::endl;
//    std::cout << "test tht: " << min(this->state(ix::th)(this->domain)) << ' ' << max(this->state(ix::th)(this->domain)) << std::endl;
//    std::cout << "test rv: " << min(this->state(ix::rv)(this->domain)) << ' ' << max(this->state(ix::rv)(this->domain)) << std::endl;
//
//    //std::cout << "test tht1: " << grad_tht[1](0, -1) << std::endl;
//    //std::cout << "test tht2: " << grad_tht[1](0, 0) << std::endl;
//    //std::cout << "test tht3: " << grad_tht[1](0, 1) << std::endl;
//    //std::cout << "test tht4: " << sgs_momenta_fluxes[0](0, 0) << std::endl;
//    //std::cout << "test tht5: " << sgs_momenta_fluxes[0](0, 1) << std::endl;
//    //std::cout << "test tht6: " << this->k_m(0, 0) << std::endl;
//    //std::cout << "test tht7: " << this->k_m(0, 1) << std::endl;
//    std::cout << "recording sgs" << std::endl;
    this->record_k_m();
    this->record_aux_dsc("tke", tke);
    this->record_aux_dsc("sgs_u_flux", sgs_momenta_fluxes[0]);
    if (ct_params_t::n_dims > 2)
    {
      this->record_aux_dsc("sgs_v_flux", sgs_momenta_fluxes[1]);
    }
    this->record_aux_dsc("p", this->Phi);
    this->record_aux_dsc("sgs_th_flux", sgs_th_flux);
    this->record_aux_dsc("sgs_rv_flux", sgs_rv_flux);
  } 

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 

    if(this->rank==0)
    {
      this->record_aux_const("fricvelsq", "sgs", params.fricvelsq);  
      this->record_aux_const("prandtl_num", "sgs", prandtl_num);  
      this->record_aux_const("karman_c", "sgs", karman_c);  
    }

    if(!params.cdrag && !params.fricvelsq) // no drag - fill tau with 0's
      formulae::stress::zero_drag_cmpct<ct_params_t::n_dims>(this->tau_srfc,
                                                             this->ijk,
                                                             this->ijkm);
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t prandtl_num, karman_c;
    real_t fricvelsq = 0; // square of friction velocity [m^2 / s^2]
  };

  // per-thread copy of params
  rt_params_t params;

  // ctor
  slvr_sgs_smg_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p),
    prandtl_num(p.prandtl_num),
    karman_c(p.karman_c),
    rcdsn_num(args.mem->tmp[__FILE__][0][0]),
    tdef_sq(args.mem->tmp[__FILE__][0][1]),
    tke(args.mem->tmp[__FILE__][0][2]),
    tmp_grad(args.mem->tmp[__FILE__][1]),
    sgs_momenta_fluxes(args.mem->tmp[__FILE__][2]),
    sgs_th_flux(args.mem->tmp[__FILE__][3][0]),
    sgs_rv_flux(args.mem->tmp[__FILE__][3][1])
  {
    if(params.fricvelsq > 0 && params.cdrag > 0)
      throw std::runtime_error("UWLCM: in SGS simulation either cdrag or fricvelsq need to be positive, not both");
  }

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 3); // rcdsn_num, tdef_sq, tke
    parent_t::alloc_tmp_vctr(mem, __FILE__); // tmp_grad
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // sgs_momenta_fluxes
    parent_t::alloc_tmp_sclr(mem, __FILE__, 2); // sgs_th/rv_flux
  }
};


// isotropic SMG version
template <class ct_params_t>
class slvr_sgs_smg_iso : public slvr_sgs_smg_common<ct_params_t>
{
  using parent_t = slvr_sgs_smg_common<ct_params_t>;
  using parent_t::parent_t; // inheriting constructors

  public:
  using real_t = typename ct_params_t::real_t;
  using ix = typename ct_params_t::ix;

  void multiply_sgs_visc() override
  {
    static_assert(static_cast<libmpdataxx::solvers::stress_diff_t>(ct_params_t::stress_diff) == libmpdataxx::solvers::compact,
                  "UWLCM smagorinsky model requires compact stress differencing");

    this->tdef_sq(this->ijk) = formulae::stress::calc_tdef_sq_cmpct<ct_params_t::n_dims>(this->tau, this->ijk);
    this->calc_rcdsn_num();
  
    this->k_m(this->ijk).reindex(this->zero) = where(
                                 this->rcdsn_num(this->ijk).reindex(this->zero) / this->prandtl_num < 1,
                                   (*this->params.mix_len_iso_sq)(this->vert_idx)
                                   * sqrt(this->tdef_sq(this->ijk).reindex(this->zero)
                                   * (1 - this->rcdsn_num(this->ijk).reindex(this->zero) / this->prandtl_num)),
                                   0
                                );
    this->k_m(this->hrzntl_slice(0)) = this->k_m(this->hrzntl_slice(1));
    this->xchng_sclr(this->k_m, this->ijk, 1);
    
   
    // calculate dissipation rate

    // TODO: c_eps should be an adjustable parameter for different cases
    real_t c_eps = 0.845;
    this->diss_rate(this->ijk).reindex(this->zero) = c_eps *
      pow3(this->k_m(this->ijk).reindex(this->zero) / (this->c_m * (*this->params.mix_len_iso_sq)(this->vert_idx)))
      * (*this->params.mix_len_iso_sq)(this->vert_idx);

    formulae::stress::multiply_tnsr_cmpct<ct_params_t::n_dims, ct_params_t::opts>(this->tau, 1.0, this->k_m, *this->mem->G, this->ijkm_sep);

    this->xchng_sgs_tnsr_offdiag(this->tau, this->tau_srfc, this->ijk, this->ijkm);
    
    //this->mem->barrier();
    //if (this->rank == 0)
    //{
    //  std::cout << "tdef_sq: " << min(tdef_sq(this->domain)) << " " << max(tdef_sq(this->domain)) << std::endl;
    //  std::cout << "rcdsn: " << min(rcdsn_num(this->domain)) << " " << max(rcdsn_num(this->domain)) << std::endl;
    //  std::cout << "k_m:   " << min(this->k_m(this->domain)) << " " << max(this->k_m(this->domain)) << std::endl;
    //  std::cout << "tau0:   " << min(this->tau[0](this->domain)) << " " << max(this->tau[0](this->domain)) << std::endl;
    //  std::cout << "tau1:   " << min(this->tau[1](this->domain)) << " " << max(this->tau[1](this->domain)) << std::endl;
    //  std::cout << "tau2:   " << min(this->tau[2](this->domain)) << " " << max(this->tau[2](this->domain)) << std::endl;
    //  //if (this->timestep % static_cast<int>(this->outfreq) == 0)
    //  //{
    //  //  std::cout << "k_m profile" << std::endl;
    //  //  for (int k = 0; k < 301; ++k)
    //  //  {
    //  //    std::cout << k << ' ' << sum(this->k_m(rng_t(0, 128), rng_t(0, 128), k)) / (129. * 129.) << std::endl;
    //  //  }
    //  //}
    //}
    //this->mem->barrier();
  }
  
  void calc_sgs_diag_fields() override
  {
    this->tke(this->ijk).reindex(this->zero) = pow2(this->k_m(this->ijk).reindex(this->zero)
                                                  / (this->c_m)) / (*this->params.mix_len_iso_sq)(this->vert_idx);
    this->calc_sgs_momenta_fluxes();
  }

  void record_k_m() override
  {
    this->record_aux_dsc("k_m", this->k_m);
  }

  // per-thread copy of params
  // rt_params_t params;
};


// anisotropic SMG version
template <class ct_params_t>
class slvr_sgs_smg_ani : public slvr_sgs_smg_common<ct_params_t>
{
  using parent_t = slvr_sgs_smg_common<ct_params_t>;
  using parent_t::parent_t; // inheriting constructors

  public:
  using real_t = typename ct_params_t::real_t;
  using ix = typename ct_params_t::ix;

  protected:

  void multiply_sgs_visc() override
  {
    static_assert(static_cast<libmpdataxx::solvers::stress_diff_t>(ct_params_t::stress_diff) == libmpdataxx::solvers::compact,
                  "UWLCM smagorinsky model requires compact stress differencing");

    this->tdef_sq(this->ijk) = formulae::stress::calc_tdef_sq_cmpct<ct_params_t::n_dims>(this->tau, this->ijk);
    this->calc_rcdsn_num();

    // Simon and Chow 2021, eqs. 14a, 14b, but using mixing length along, not normal!
    // l^-2 = (kz)^-2 + (Cs dl)^-2 = A
    // l = A^(-1/2)
    // l^2 = 1/A
    // horizontal assuming dx=dy in 3D
    this->k_m[0](this->ijk).reindex(this->zero) = where(
                                 this->rcdsn_num(this->ijk).reindex(this->zero) / this->prandtl_num < 1,
                                   (*this->params.mix_len_hori_sq)(this->vert_idx)
                                   * sqrt(this->tdef_sq(this->ijk).reindex(this->zero)
                                   * (1 - this->rcdsn_num(this->ijk).reindex(this->zero) / this->prandtl_num)),
                                   0
                                );
    // min/max as in Simon and Chow 2021
    // TODO: do it all in the single loop above
    this->k_m[0](this->ijk) = where(this->k_m[0](this->ijk) > 0.1 * pow(this->params.di, 2) / this->dt, 0.1 * pow(this->params.di, 2) / this->dt, this->k_m[0](this->ijk));
    this->k_m[0](this->ijk) = where(this->k_m[0](this->ijk) < 1e-6 * pow(this->params.di, 2),  1e-6 * pow(this->params.di, 2), this->k_m[0](this->ijk));

    this->k_m[0](this->hrzntl_slice(0)) = this->k_m[0](this->hrzntl_slice(1));
    this->xchng_sclr(this->k_m[0], this->ijk, 1);

    // vertical
    this->k_m[1](this->ijk).reindex(this->zero) = where(
                                 this->rcdsn_num(this->ijk).reindex(this->zero) / this->prandtl_num < 1,
                                   (*this->params.mix_len_vert_sq)(this->vert_idx)
                                   * sqrt(this->tdef_sq(this->ijk).reindex(this->zero)
                                   * (1 - this->rcdsn_num(this->ijk).reindex(this->zero) / this->prandtl_num)),
                                   0
                                );
    // min/max as in Simon and Chow 2021
    // TODO: do it all in the single loop above
    this->k_m[1](this->ijk) = where(this->k_m[1](this->ijk) > 0.1 * pow(this->params.dz, 2) / this->dt, 0.1 * pow(this->params.dz, 2) / this->dt, this->k_m[1](this->ijk));
    this->k_m[1](this->ijk) = where(this->k_m[1](this->ijk) < 1e-6 * pow(this->params.dz, 2),  1e-6 * pow(this->params.dz, 2), this->k_m[1](this->ijk));

    this->k_m[1](this->hrzntl_slice(0)) = this->k_m[1](this->hrzntl_slice(1));
    this->xchng_sclr(this->k_m[1], this->ijk, 1);

//    std::cerr << "k_m[0]: " << this->k_m[0] << std::endl;
//    std::cerr << "k_m[1]: " << this->k_m[1] << std::endl;
   
    // calculate dissipation rate

    // TODO: c_eps should be an adjustable parameter for different cases
    // NOTE: diss rate is calculated only using horizontal turbulence, which is assumed to dominate (as dz is expected to be smaller)
    real_t c_eps = 0.845;
    this->diss_rate(this->ijk).reindex(this->zero) = c_eps *
      pow3(this->k_m[0](this->ijk).reindex(this->zero) / (this->c_m * (*this->params.mix_len_hori_sq)(this->vert_idx)))
      * (*this->params.mix_len_hori_sq)(this->vert_idx);

    formulae::stress::multiply_tnsr_cmpct<ct_params_t::n_dims, ct_params_t::opts>(this->tau, 1.0, this->k_m, *this->mem->G, this->ijkm_sep);

    this->xchng_sgs_tnsr_offdiag(this->tau, this->tau_srfc, this->ijk, this->ijkm);    
  }


  void calc_sgs_diag_fields() override
  {
    // NOTE: calculated only using horizontal turbulence, which is assumed to dominate (as dz is expected to be smaller)
    this->tke(this->ijk).reindex(this->zero) = pow2(this->k_m[0](this->ijk).reindex(this->zero)
                                                  / (this->c_m)) / (*this->params.mix_len_hori_sq)(this->vert_idx);
    this->calc_sgs_momenta_fluxes();
  }

  void record_k_m() override
  {
    this->record_aux_dsc("k_m[0]", this->k_m[0]);
    this->record_aux_dsc("k_m[1]", this->k_m[1]);
  }
};


// primary template to select appropriate SGS solver depending on ct_params_t::sgs_scheme
template <class ct_params_t, class enableif = void>
class slvr_sgs
{};

// ILES version - slvr_common already has all needed functionality
template <class ct_params_t>
class slvr_sgs<
  ct_params_t,
  typename std::enable_if<ct_params_t::sgs_scheme == libmpdataxx::solvers::iles>::type
> : public slvr_common<ct_params_t>
{
  using parent_t = slvr_common<ct_params_t>;
  using parent_t::parent_t; // inheriting constructors
};

// isotropic SMG version
template <class ct_params_t>
class slvr_sgs<
  ct_params_t,
  typename std::enable_if<ct_params_t::sgs_scheme == libmpdataxx::solvers::smg>::type
> : public slvr_sgs_smg_iso<ct_params_t>
{
  using parent_t = slvr_sgs_smg_iso<ct_params_t>;
  using parent_t::parent_t; // inheriting constructors
};

// anisotropic SMG version
template <class ct_params_t>
class slvr_sgs<
  ct_params_t,
  typename std::enable_if<ct_params_t::sgs_scheme == libmpdataxx::solvers::smgani>::type
> : public slvr_sgs_smg_ani<ct_params_t>
{
  using parent_t = slvr_sgs_smg_ani<ct_params_t>;
  using parent_t::parent_t; // inheriting constructors
};


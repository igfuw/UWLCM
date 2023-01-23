#include <chrono>

template <class solver_t>
class exec_timer : public solver_t
{
  using parent_t = solver_t;

  private:
  const unsigned long nt;

  setup::clock::time_point tbeg_step, tend_step,
                           tbeg_aux , tend_aux ,
                           tbeg_loop, tend_loop;

  setup::timer tloop, trecord_all, thas, thads, thps, thps_has, thas_hads, thads_hps, thmas, thmps;

  void get_step_time(setup::timer &timer)
  {
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tend_step = setup::clock::now();
      timer += std::chrono::duration_cast<std::chrono::milliseconds>( tend_step - tbeg_step );
      tbeg_step = setup::clock::now();
    }
    this->mem->barrier();
  }

  void start_aux_clock()
  {
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tbeg_aux = setup::clock::now();
    }
    this->mem->barrier();
  }

  void stop_aux_clock(setup::timer &timer)
  {
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tend_aux = setup::clock::now();
      timer += std::chrono::duration_cast<std::chrono::milliseconds>( tend_aux - tbeg_aux );
    }
    this->mem->barrier();
  }

  protected:

  void hook_ante_loop(int nt) override
  {
    parent_t::hook_ante_loop(nt);
    this->mem->barrier();
    if (this->rank == 0) 
    {   
      tbeg_loop = setup::clock::now();
      trecord_all = setup::timer::zero(); // reset to 0, because we only want record all done in loop, not the one in ante_loop 
    }
    this->mem->barrier();
  }

  void hook_ante_step() override
  {
    get_step_time(thps_has);
    parent_t::hook_ante_step();
    get_step_time(thas);
  }

  void hook_ante_delayed_step() override
  {
    get_step_time(thas_hads);
    parent_t::hook_ante_delayed_step();
    get_step_time(thads);
  }

  void hook_post_step() override
  {
    get_step_time(thads_hps);
    parent_t::hook_post_step();
    get_step_time(thps);

    // there's no hook_post_loop, so we imitate it here to write out computation times, TODO: move to destructor?
    if(this->timestep == nt) // timestep incremented before post_step
    {
      if (this->rank == 0)
      {
        tend_loop = setup::clock::now();
        tloop = std::chrono::duration_cast<std::chrono::milliseconds>( tend_loop - tbeg_loop );

        // calculate CPU/GPU times and concurrency, valid only for async runs and not taking into account diagnostics in record_all
        setup::timer  tsync_in = parent_t::tsync,
                      tgpu = parent_t::tasync_wait_in_record_all + parent_t::tsync_wait + parent_t::tasync_wait + tsync_in, // time of pure GPU calculations (= wait time of CPU)
                      tcpugpu = tsync_in + parent_t::tasync_gpu + parent_t::tsync_gpu - tgpu, // time of concurrent CPU and GPU calculations (= total time of GPU calculations - tgpu)
                      tcpu = tloop - tgpu - tcpugpu;
         
        std::cout <<  "wall time in milliseconds: " << std::endl
          << "loop:                            " << tloop.count() << std::endl
          << "  hook_ante_step:                  " << thas.count() << " ("<< setup::real_t(thas.count())/tloop.count()*100 <<"%)" << std::endl
          << "    hook_mixed_rhs_ante_step:        " << thmas.count() << " ("<< setup::real_t(thmas.count())/tloop.count()*100 <<"%)" << std::endl
          << "      async_wait:                    " << parent_t::tasync_wait.count() << " ("<< setup::real_t(parent_t::tasync_wait.count())/tloop.count()*100 <<"%)" << std::endl
          << "      sync:                            " << parent_t::tsync.count() << " ("<< setup::real_t(parent_t::tsync.count())/tloop.count()*100 <<"%)" << std::endl
          << "  step:                            " << thas_hads.count() << " ("<< setup::real_t(thas_hads.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_ante_delayed_step:          " << thads.count() << " ("<< setup::real_t(thads.count())/tloop.count()*100 <<"%)" << std::endl
          << "    sync_wait:                       " << parent_t::tsync_wait.count() << " ("<< setup::real_t(parent_t::tsync_wait.count())/tloop.count()*100 <<"%)" << std::endl
          << "    async:                           " << parent_t::tasync.count() << " ("<< setup::real_t(parent_t::tasync.count())/tloop.count()*100 <<"%)" << std::endl
          << "  delayed step:                    " << thads_hps.count() << " ("<< setup::real_t(thads_hps.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_post_step:                  " << thps.count() << " ("<< setup::real_t(thps.count())/tloop.count()*100 <<"%)" << std::endl
          << "    hook_mixed_rhs_post_step:        " << thmps.count() << " ("<< setup::real_t(thmps.count())/tloop.count()*100 <<"%)" << std::endl
          << "    record_all (in loop):            " << trecord_all.count() << " ("<< setup::real_t(trecord_all.count())/tloop.count()*100 <<"%)" << std::endl
          << "      async_wait in record_all:      " << parent_t::tasync_wait_in_record_all.count() << " ("<< setup::real_t(parent_t::tasync_wait_in_record_all.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_post_step->hook_ante_step:  " << thps_has.count() << " ("<< setup::real_t(thps_has.count())/tloop.count()*100 <<"%)" << std::endl;

          std::cout << std::endl
          << "CPU/GPU concurrency stats, only make sense for async lgrngn runs" << std::endl
          << "and does not take into account GPU time in record_all, so most accurate without diag:" << std::endl
          << "  pure CPU calculations: " << tcpu.count() << " ("<< setup::real_t(tcpu.count())/tloop.count()*100 <<"%)" << std::endl
          << "  pure GPU calculations: " << tgpu.count() << " ("<< setup::real_t(tgpu.count())/tloop.count()*100 <<"%)" << std::endl
          << "  concurrent CPU&GPU:    " << tcpugpu.count() << " ("<< setup::real_t(tcpugpu.count())/tloop.count()*100 <<"%)" << std::endl
          << "  tsync_gpu:  " << parent_t::tsync_gpu.count() << " ("<< setup::real_t(parent_t::tsync_gpu.count())/tloop.count()*100 <<"%)" << std::endl
          << "  tasync_gpu: " << parent_t::tasync_gpu.count() << " ("<< setup::real_t(parent_t::tasync_gpu.count())/tloop.count()*100 <<"%)" << std::endl;
      }
    }
  }

  void record_all() override
  {
    assert(this->rank == 0);

    tbeg_aux = setup::clock::now();
    parent_t::record_all();
    tend_aux = setup::clock::now();
    trecord_all = std::chrono::duration_cast<std::chrono::milliseconds>( tend_aux - tbeg_aux );
  }

  void hook_mixed_rhs_ante_step() override
  {
    start_aux_clock();
    parent_t::hook_mixed_rhs_ante_step();
    stop_aux_clock(thmas);
  }

  void hook_mixed_rhs_post_step() override
  {
    start_aux_clock();
    parent_t::hook_mixed_rhs_post_step();
    stop_aux_clock(thmps);
  }

  public:

  // ctor
  exec_timer( 
    typename parent_t::ctor_args_t args, 
    const typename parent_t::rt_params_t &p
  ) : 
    parent_t(args, p),
    nt(p.user_params.nt)
  {}  
};

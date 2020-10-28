#include <chrono>

template <class solver_t>
class exec_timer : public solver_t
{
  using parent_t = solver_t;

  private:
  const unsigned long nt;

  typename parent_t::clock::time_point tbeg_step, tend_step,
                                       tbeg_aux , tend_aux ,
                                       tbeg_loop, tend_loop;

  typename parent_t::timer tloop, trecord_all, thas, thads, thps, thps_has, thas_hads, thads_hps, thmas, thmps;

  void get_step_time(typename parent_t::timer &timer)
  {
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tend_step = parent_t::clock::now();
      timer += std::chrono::duration_cast<std::chrono::milliseconds>( tend_step - tbeg_step );
      tbeg_step = parent_t::clock::now();
    }
    this->mem->barrier();
  }

  void start_aux_clock()
  {
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tbeg_aux = parent_t::clock::now();
    }
    this->mem->barrier();
  }

  void stop_aux_clock(typename parent_t::timer &timer)
  {
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tend_aux = parent_t::clock::now();
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
      tbeg_loop = parent_t::clock::now();
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
        tend_loop = parent_t::clock::now();
        tloop = std::chrono::duration_cast<std::chrono::milliseconds>( tend_loop - tbeg_loop );

        std::cout <<  "wall time in milliseconds: " << std::endl
          << "loop:                            " << tloop.count() << std::endl
          << "  hook_ante_step:                  " << thas.count() << " ("<< setup::real_t(thas.count())/tloop.count()*100 <<"%)" << std::endl
          << "    hook_mixed_rhs_ante_step:        " << thmas.count() << " ("<< setup::real_t(thmas.count())/tloop.count()*100 <<"%)" << std::endl
          << "    async_wait:                      " << parent_t::tasync_wait.count() << " ("<< setup::real_t(parent_t::tasync_wait.count())/tloop.count()*100 <<"%)" << std::endl
          << "    sync:                            " << parent_t::tsync.count() << " ("<< setup::real_t(parent_t::tsync.count())/tloop.count()*100 <<"%)" << std::endl
          << "  step:                            " << thas_hads.count() << " ("<< setup::real_t(thas_hads.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_ante_delayed_step:          " << thads.count() << " ("<< setup::real_t(thads.count())/tloop.count()*100 <<"%)" << std::endl
          << "    sync_wait:                       " << parent_t::tsync_wait.count() << " ("<< setup::real_t(parent_t::tsync_wait.count())/tloop.count()*100 <<"%)" << std::endl
          << "    async:                           " << parent_t::tasync.count() << " ("<< setup::real_t(parent_t::tasync.count())/tloop.count()*100 <<"%)" << std::endl
          << "  delayed step:                    " << thads_hps.count() << " ("<< setup::real_t(thads_hps.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_post_step:                  " << thps.count() << " ("<< setup::real_t(thps.count())/tloop.count()*100 <<"%)" << std::endl
          << "    hook_mixed_rhs_post_step:        " << thmps.count() << " ("<< setup::real_t(thmps.count())/tloop.count()*100 <<"%)" << std::endl
          << "    record_all:                      " << trecord_all.count() << " ("<< setup::real_t(trecord_all.count())/tloop.count()*100 <<"%)" << std::endl
          << "      async_wait in record_all:      " << parent_t::tasync_wait_in_record_all.count() << " ("<< setup::real_t(parent_t::tasync_wait_in_record_all.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_post_step->hook_ante_step:  " << thps_has.count() << " ("<< setup::real_t(thps_has.count())/tloop.count()*100 <<"%)" << std::endl;
      }
    }
  }

  void record_all() override
  {
    assert(this->rank == 0);

    tbeg_aux = parent_t::clock::now();
    parent_t::record_all();
    tend_aux = parent_t::clock::now();
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
    nt(p.nt)
  {}  
};

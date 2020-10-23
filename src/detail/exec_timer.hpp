#include <chrono>

template <class solver_t>
class exec_timer : public solver_t
{
  using parent_t = solver_t;

  private:
  typename parent_t::clock::time_point tbeg,      tend,
                                       tbeg_loop, tend_loop;

  typename parent_t::timer tloop, trecord_all, thas, thads, thps, thps_has, thas_hads, thads_hps;

  void time(std::chrono::milliseconds &timer)
  {
    this->mem->barrier();
    if (this->rank == 0) 
    {
      tend = parent_t::clock::now();
      timer += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      tbeg = parent_t::clock::now();
    }
    this->mem->barrier();
  }


  protected:

  void hook_ante_loop()
  {
    parent_t::hook_ante_loop();
    this->mem->barrier();
    if (this->rank == 0) 
      tbeg_loop = parent_t::clock::now();
    this->mem->barrier();
  }

  void hook_ante_step()
  {
    timer(thps_has);
    parent_t::hook_ante_step();
    timer(thas);
  }

  void hook_ante_delayed_step()
  {
    timer(thas_hads);
    parent_t::hook_ante_delayed_step();
    timer(thads);
  }

  void hook_post_step()
  {
    timer(thads_hps);
    parent_t::hook_post_step();
    timer(thps);

    // there's no hook_post_loop, so we imitate it here to write out computation times, TODO: move to destructor?
    if(this->timestep == this->params.nt) // timestep incremented before post_step
    {
      if (this->rank == 0)
      {
        tend_loop = parent_t::clock::now();
        tloop = std::chrono::duration_cast<std::chrono::milliseconds>( tend_loop - tbeg_loop );

        std::cout <<  "wall time in milliseconds: " << std::endl
          << "loop:                            " << tloop.count() << std::endl
          << "  hook_ante_step:                  " << thas.count() << " ("<< setup::real_t(thas.count())/tloop.count()*100 <<"%)" << std::endl
          << "    async_wait:                      " << parent_t::tasync_wait.count() << " ("<< setup::real_t(parent_t::tasync_wait.count())/tloop.count()*100 <<"%)" << std::endl
          << "    update_rhs0:                     " << parent_t::tupdate_rhs0.count() << " ("<< setup::real_t(parent_t::tupdate_rhs0.count())/tloop.count()*100 <<"%)" << std::endl
          << "    apply_rhs0:                      " << parent_t::tapply_rhs0.count() << " ("<< setup::real_t(parent_t::tapply_rhs0.count())/tloop.count()*100 <<"%)" << std::endl
          << "    sync:                            " << parent_t::tsync.count() << " ("<< setup::real_t(parent_t::tsync.count())/tloop.count()*100 <<"%)" << std::endl
          << "  step:                            " << thas_hads.count() << " ("<< setup::real_t(thas_hads.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_ante_delayed_step:          " << thads.count() << " ("<< setup::real_t(thads.count())/tloop.count()*100 <<"%)" << std::endl
          << "    sync_wait:                       " << parent_t::tsync_wait.count() << " ("<< setup::real_t(parent_t::tsync_wait.count())/tloop.count()*100 <<"%)" << std::endl
          << "    async:                           " << parent_t::tasync.count() << " ("<< setup::real_t(parent_t::tasync.count())/tloop.count()*100 <<"%)" << std::endl
          << "  delayed step:                    " << thads_hps.count() << " ("<< setup::real_t(thads_hps.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_post_step:                  " << thps.count() << " ("<< setup::real_t(thps.count())/tloop.count()*100 <<"%)" << std::endl
          << "    update_rhs1:                     " << parent_t::tupdate_rhs1.count() << " ("<< setup::real_t(parent_t::tupdate_rhs1.count())/tloop.count()*100 <<"%)" << std::endl
          << "    apply_rhs1:                      " << parent_t::tapply_rhs1.count() << " ("<< setup::real_t(parent_t::tapply_rhs1.count())/tloop.count()*100 <<"%)" << std::endl
          << "    record_all:                      " << trecord_all.count() << " ("<< setup::real_t(trecord_all.count())/tloop.count()*100 <<"%)" << std::endl
          << "      async_wait in record_all:      " << parent_t::tasync_wait_in_record_all.count() << " ("<< setup::real_t(parent_t::tasync_wait_in_record_all.count())/tloop.count()*100 <<"%)" << std::endl
          << "  hook_post_step->hook_ante_step:  " << thps.count() << " ("<< setup::real_t(thps.count())/tloop.count()*100 <<"%)" << std::endl
          << std::endl
          << "  update_rhs parts from update_rhs0 + update_rhs1: " << std::endl
          << "    update_rhs in slvr_common: " << parent_t::tupdate_rhs_slvr_common.count() << " ("<< setup::real_t(parent_t::tupdate_rhs_slvr_common.count())/tloop.count()*100 <<"%)" << std::endl
          << "    update_rhs in slvr_sgs:    " << parent_t::tupdate_rhs_slvr_sgs.count() << " ("<< setup::real_t(parent_t::tupdate_rhs_slvr_sgs.count())/tloop.count()*100 <<"%)" << std::endl;
      }
    }
  }

  void record_all()
  {
    assert(this->rank == 0);

    tbeg = parent_t::clock::now();
    parent_t::record_all();
    tend = parent_t::clock::now();
    trecord_all = std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
  }

};

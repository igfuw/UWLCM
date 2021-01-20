// function that calculates execution time of any other member function called via ptr
#pragma once

#if defined(UWLCM_TIMING)
  template<class clock, class timer, class F, class ptr, typename... Args>
  timer func_time(F func, ptr p, Args&&... args){
      timer t1=clock::now();
      p->func(std::forward<Args>(args)...);
      return std::chrono::duration_cast<timer>(clock::now()-t1);
  }
#else
  template<class clock, class timer, class F, class ptr, typename... Args>
  timer func_time(F func, ptr p, Args&&... args){
      p->func(std::forward<Args>(args)...);
      return timer();
  }
#endif

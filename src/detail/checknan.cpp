#pragma once

#ifdef NDEBUG
#define nancheck(arr, name) ((void)0)
#else
#define nancheck(arr, name) {nancheck_hlprs::nancheck_hlpr(arr, name);}
#endif

#ifdef NDEBUG
#define nancheck2(arrcheck, arrout, name) ((void)0)
#else
#define nancheck2(arrcheck, arrout, name) {nancheck_hlprs::nancheck2_hlpr(arrcheck, arrout, name);}
#endif

#ifdef NDEBUG
#define negcheck(arr, name) ((void)0)
#else
#define negcheck(arr, name) {nancheck_hlprs::negcheck_hlpr(arr, name);}
#endif

// actually not to zero, but to 1e-10 (we need rv>0 in libcloud and cond substepping numerical errors colud lead to rv<0 if we would set it here to 0)
// we don't make it a critical section, because it is also used in production runs
#define negtozero(arr, name) {if(min(arr) < 0.) {\
                               std::cout << "A negative number detected in: " << name << std::endl;\
                               std::cout << arr;\
                               std::cout << "CHEATING: turning negative values to small positive values" << std::endl;\
                               arr = where(arr <= 0., 1e-10, arr);\
                               }}

#ifndef NDEBUG
namespace nancheck_hlprs
{
  template<class arr_t>
  void nancheck_hlpr(const arr_t &arr, std::string name)
  {
    #pragma omp critical
    {
      if(!std::isfinite(sum(arr))) 
      {
        std::cout << "A not-finite number detected in: " << name << std::endl;
        std::cout << arr;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void nancheck2_hlpr(const arr_t &arrcheck, const arr_t &arrout, std::string name)
  {
    #pragma omp critical
    {
      if(!std::isfinite(sum(arrcheck))) 
      {
        std::cout << "A not-finite number detected in: " << name << std::endl;
        std::cout << arrcheck;
        std::cout << arrout;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void negcheck_hlpr(const arr_t &arr, std::string name)
  {
    #pragma omp critical
    {
      if(min(arr) < 0.) 
      {
        std::cout << "A negative number detected in: " << name << std::endl;
        std::cout << arr;
        assert(0);
      }
    }
  }
};
#endif

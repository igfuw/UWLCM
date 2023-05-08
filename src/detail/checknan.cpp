#pragma once

#define NEGTOZERO_SET_VAR 1e-4

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

#ifdef NDEBUG
#define negcheck2(arrcheck, arrout, name) ((void)0)
#else
#define negcheck2(arrcheck, arrout, name) {nancheck_hlprs::negcheck2_hlpr(arrcheck, arrout, name);}
#endif

#ifdef NDEBUG
// actually not to zero, but to NEGTOZERO_SET_VAR (we need rv>0 in libcloud and cond substepping numerical errors colud lead to rv<0 if we would set it here to 0)
#define negtozero(arr, name) {arr = where(arr <= 0., NEGTOZERO_SET_VAR, arr);}
#else
#define negtozero(arr, name) {nancheck_hlprs::negtozero_hlpr(arr, name);}
#endif

#ifdef NDEBUG
// same as above, but with printing additional arrays
#define negtozero2(arr, name, outarrs, outnames) {arr = where(arr <= 0., NEGTOZERO_SET_VAR, arr);}
#else
#define negtozero2(arr, name, outarrs, outnames) {nancheck_hlprs::negtozero_hlpr(arr, name, outarrs, outnames);}
#endif

#ifndef NDEBUG
namespace nancheck_hlprs
{
  template<class arr_t>
  void nancheck_hlpr(const arr_t &arr, const std::string &name)
  {
    if(!std::isfinite(sum(arr)))
    {
      #pragma omp critical
      {
        std::cerr << "A not-finite number detected in: " << name << std::endl;
        std::cerr << arr;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void nancheck2_hlpr(const arr_t &arrcheck, const arr_t &arrout, const std::string &name)
  {
    if(!std::isfinite(sum(arrcheck)))
    {
      #pragma omp critical
      {
        std::cerr << "A not-finite number detected in: " << name << std::endl;
        std::cerr << arrcheck;
        std::cerr << arrout;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void negcheck_hlpr(const arr_t &arr, const std::string &name)
  {
    if(min(arr) < 0.)
    {
      #pragma omp critical
      {
        std::cerr << "A negative number detected in: " << name << std::endl;
        std::cerr << arr;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void negcheck2_hlpr(const arr_t &arrcheck, const arr_t &arrout,const std::string &name)
  {
    if(min(arrcheck) < 0.) 
    {
      #pragma omp critical
      {
        std::cout << "A negative number detected in: " << name << std::endl;
        std::cout << arrcheck;
        std::cout << arrout;
        assert(0);
      }
    }
  }

  template<class arr_t>
  void negtozero_hlpr(arr_t arr, const std::string &name)
  {
    auto minVal = min(arr);
    if(minVal < 0.) 
    {
      auto minIdx = minIndex(arr);
      #pragma omp critical
      {
        std::cerr << count(arr<=0.) << " non-positive numbers detected in: " << name  << ". Minval: " << minVal << " minIndex: " << minIdx << std::endl;
        std::cerr << "CHEATING: turning non-positive values to small positive values (" + std::to_string(NEGTOZERO_SET_VAR) + ")" << std::endl;
      }
      arr = where(arr <= 0., NEGTOZERO_SET_VAR, arr);
    }
  }

  // as above but with printing of additional arrays
  template<class arr_t, class arrvec_t>
  void negtozero_hlpr(arr_t arr, const std::string &name, /*shuld be const*/ arrvec_t outarrs, const std::vector<std::string> &outnames)
  {
    assert(outarrs.size() == outnames.size());

    auto minVal = min(arr);
    if(minVal < 0.) 
    {
      auto minIdx = minIndex(arr);
      #pragma omp critical
      {
        std::cerr << count(arr<=0.) << " non-positive numbers detected in: " << name  << ". Minval: " << minVal << " minIndex: " << minIdx << std::endl;
	for(int c=0; c<outarrs.size(); ++c)
	{
	  std::cerr << outnames.at(c) << "(minIndex): " << outarrs[c](minIdx) << std::endl;
	}
        std::cerr << "CHEATING: turning non-positive values to small positive values (" + std::to_string(NEGTOZERO_SET_VAR) + ")" << std::endl;
      }
      arr = where(arr <= 0., NEGTOZERO_SET_VAR, arr);
    }
  }


};
#endif

// function that calculates execution time of any other member function called via ptr
#pragma once

// async_forwarder code taken from https://kholdstare.github.io/technical/2012/12/18/perfect-forwarding-to-async-2.html (C) Alexander Kondratskiy
// it is used to pass any type of reference (lvalue or rvalue) to std::async
template <typename T>
class async_forwarder
{
    // Store value directly
    T val_;

public:
    /**
     * Move an rvalue of T into the wrapper,
     * incurring no copies.
     */
    async_forwarder(T&& t) 
     : val_(std::move(t)) { }

    // ensure no copies are made
    async_forwarder(async_forwarder const& other) = delete;

    // move constructor
    async_forwarder(async_forwarder&& other)
        : val_(std::move(other.val_)) { }

    // Move the value out.
    // Note: can only occur once!
    operator T&& ()       { return std::move(val_); }
    operator T&& () const { return std::move(val_); }
};

// This particular specialization
// is essentially std::ref
template <typename T>
class async_forwarder<T&>
{
    T& val_;

public:
    /**
     * Wrap the reference when passed an lvalue reference,
     * to fool std::async
     */
    async_forwarder(T& t) : val_(t) { }

    // ensure no copies are made
    async_forwarder(async_forwarder const& other) = delete;

    // move constructor
    async_forwarder(async_forwarder&& other)
        : val_(other.val_) { }

    // User-defined conversion that automatically
    // converts to the appropriate type
    operator T&       ()       { return val_; }
    operator T const& () const { return val_; }
};


#if defined(UWLCM_TIMING)
  template<class F, class ptr, typename... Args>
  setup::timer func_time(F func, ptr p, Args&&... args){
    auto t1=setup::clock::now();
    (p->*func)(std::forward<Args>(args)...);
    return std::chrono::duration_cast<setup::timer>(setup::clock::now()-t1);
  }
#endif

#if defined(UWLCM_TIMING)
  template<class F, class ptr, typename... Args>
  std::future<setup::timer> async_launcher(F func, ptr p, Args&&... args) // func and p are pointers, so their copies are lightweight
  {
    return std::async(
             std::launch::async,
             func_time<F, ptr, Args...>,
             func, 
             p,
             async_forwarder<Args>(std::forward<Args>(args))... // ATTENTION! args are passed by reference to async
           );
  }
#else
  template<class F, class ptr, typename... Args>
  std::future<void> async_launcher(F func, ptr p, Args&&... args) // func and p are pointers, so their copies are lightweight
  {
    return std::async(
             std::launch::async,
             func, 
             p,
             async_forwarder<Args>(std::forward<Args>(args))... // ATTENTION! args are passed by reference to async
           );
  }
#endif

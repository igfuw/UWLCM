#pragma once

#ifdef NDEBUG
#define nancheck(arr, name) ((void)0)
#else
#define nancheck(arr, name) {if(!std::isfinite(sum(arr))) {\
                               std::cout << "A not-finite number detected in: " << name << std::endl;\
                               std::cout << arr;\
                               abort();}}
#endif

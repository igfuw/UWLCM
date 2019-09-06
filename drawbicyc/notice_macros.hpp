// error reporting

#include <iostream>

#define error_macro(msg) \
{ \
  std::cerr << "error: " << msg << std::endl; \
  throw std::exception(); \
}

#define notice_macro(msg) \
{ \
  std::cerr << " info: " << msg << std::endl; \
}

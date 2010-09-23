#pragma once
#ifndef STRING_UTILS_HPP
#define STRING_UTILS_HPP

#include <sstream>
#include <string>

/// collection of various string utilities
class StringUtils
{
 public:
  /* *****************
  * member functions *
  *******************/

  ///casts C++ data types in std::string
  template<typename T>
  static inline std::string stringify(const T& x)
  {
    std::ostringstream o;
    o << x;
    return o.str();
  }
};
#endif //  #ifndef STRING_UTILS_HPP

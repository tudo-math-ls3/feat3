#pragma once
#ifndef KERNEL_UTIL_INHERITANCE_HPP
#define KERNEL_UTIL_INHERITANCE_HPP 1

/**
* \file inheritance.hpp
*
* \brief Static compile-time inheritance analysis.
*
* \author Dirk Ribbrock
*/

namespace FEAST
{
  /// \cond
  template <typename T, typename U>
  class Conversion
  {
    typedef char small;
    class Big {char dummy[2]; };
    static small Test (const U& );
    static Big Test(...);
    static T MakeT();

  public:
    enum {exists = sizeof(Test(MakeT())) == sizeof(small) };
    enum {sameType = false};
  };

  template <typename T>
  class Conversion<T, T>
  {
  public:
    enum {exists = 1, sameType = 1};
  };
  /// \endcond
}

/// checks if U inherits from T
#define SUPERSUBCLASS(T, U) \
  (::FEAST::Conversion<const U*, const T*>::exists && \
   !::FEAST::Conversion<const T*, const void*>::sameType)
#endif // KERNEL_UTIL_INHERITANCE_HPP

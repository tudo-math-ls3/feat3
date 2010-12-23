#pragma once
#ifndef UTIL_INHERITANCE_HPP
/// Header guard
#define UTIL_INHERITANCE_HPP 1

/**
 * \file inheritance.hpp
 *
 * \brief Static compile-time inheritance analysis.
 *
 * \author Dirk Ribbrock
 */

namespace FEAST
{
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
}

/**
 * \brief Check that U inherits from T.
 */
#define SUPERSUBCLASS(T, U) \
  (::FEAST::Conversion<const U*, const T*>::exists && \
   !::FEAST::Conversion<const T*, const void*>::sameType)
#endif // UTIL_INHERITANCE_HPP

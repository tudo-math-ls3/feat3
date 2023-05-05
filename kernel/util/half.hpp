// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_HALF_HPP
#define KERNEL_UTIL_HALF_HPP 1

// if cuda supports fp16 arithmetics, include the cuda sdk fp16 header for cpu and gpu computing in datatype __half
#ifdef FEAT_HAVE_HALFMATH
#include <cuda_fp16.h>
#endif // FEAT_HAVE_HALFMATH

namespace FEAT
{
#if defined(FEAT_HAVE_HALFMATH) || defined(DOXYGEN)
  /**
   * \brief Half data type.
   *
   * Half data type derived from nvidias cuda __half data type.
   *
   * \todo add overloads/specializations in math.hpp for this type
   *
   * \author Dirk Ribbrock
   */
  class Half : public __half
  {
  public:
    // use constructors from __half parent class
    using __half::__half;

    Half(const int & other) : __half(double(other))
    {
    }

    Half(const unsigned int & other) : __half(double(other))
    {
    }

    Half(const long & other) : __half(double(other))
    {
    }

    Half(const unsigned long & other) : __half(double(other))
    {
    }

    Half& operator+=(const Half& other)
    {
      float dother(other);
      float dthis(*this);
      dthis += dother;
      *this = dthis;
      return *this;
    }

    Half& operator*=(const Half& other)
    {
      float dother(other);
      float dthis(*this);
      dthis *= dother;
      *this = dthis;
      return *this;
    }

    Half& operator/=(const Half& other)
    {
      float dother(other);
      float dthis(*this);
      dthis /= dother;
      *this = dthis;
      return *this;
    }

    Half& operator-=(const Half& other)
    {
      float dother(other);
      float dthis(*this);
      dthis -= dother;
      *this = dthis;
      return *this;
    }

    Half operator-()
    {
      float dthis(*this);
      dthis = -dthis;
      return Half(dthis);
    }
  };
#endif // FEAT_HAVE_HALFMATH

} // namespace FEAT

#endif // KERNEL_UTIL_HALF_HPP

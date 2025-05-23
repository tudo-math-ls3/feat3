// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SCALE_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SCALE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {

      template <typename DT_>
      void Scale::value_generic(DT_ * r, const DT_ * const x, const DT_ s, const Index size)
      {
        if (x == r)
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < size ; ++i)
          {
            r[i] *= s;
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < size ; ++i)
          {
            r[i] = x[i] * s;
          }
        }
      }

      template <typename ValueType_>
      void Scale::value_blocked_generic(ValueType_ * r, const ValueType_ * const x, const ValueType_ s, const Index size)
      {
        if (x == r)
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < size ; ++i)
          {
            for(int j = 0; j < ValueType_::n; ++j)
            {
              r[i][j] *= s[j];
            }
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < size ; ++i)
          {
            for(int j = 0; j < ValueType_::n; ++j)
            {
              r[i][j] = x[i][j] * s[j];
            }
          }
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_SCALE_GENERIC_HPP

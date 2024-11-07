// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      DT_ DotProduct::value_generic(const DT_ * const x, const DT_ * const y, const Index size)
      {
        DT_ r(0);

        if(x == y)
        {
          FEAT_PRAGMA_OMP(parallel for reduction(+:r))
          for (Index i = 0 ; i < size ; ++i)
          {
            r += x[i] * x[i];
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for reduction(+:r))
          for (Index i = 0 ; i < size ; ++i)
          {
            r += x[i] * y[i];
          }
        }

        return r;
      }

      template <typename ValueType_>
      ValueType_ DotProduct::value_blocked_generic(const ValueType_ * const x, const ValueType_ * const y, const Index size)
      {
        ValueType_ r(0);

        if(x == y)
        {

          for (Index i(0) ; i < size ; ++i)
          {
            for(int j(0); j < ValueType_::n; ++j) {
              r[j] += x[i][j] * x[i][j];
            }
          }
        }
        else
        {

          for (Index i(0) ; i < size ; ++i)
          {
            for(int j(0); j < ValueType_::n; ++j) {
              r[j] += x[i][j] * y[i][j];
            }
          }
        }

        return r;
      }

      template <typename DT_>
      DT_ TripleDotProduct::value_generic(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
      {
        DT_ r(0);

        if (x == y)
        {
          FEAT_PRAGMA_OMP(parallel for reduction(+:r))
          for (Index i = 0 ; i < size ; ++i)
            r += x[i] * x[i] * z[i];
        }
        else if (x == z)
        {
          FEAT_PRAGMA_OMP(parallel for reduction(+:r))
          for (Index i = 0 ; i < size ; ++i)
            r += x[i] * x[i] * y[i];
        }
        else if (y == z)
        {
          FEAT_PRAGMA_OMP(parallel for reduction(+:r))
          for (Index i = 0 ; i < size ; ++i)
            r += x[i] * y[i] * y[i];
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for reduction(+:r))
          for (Index i = 0 ; i < size ; ++i)
            r += x[i] * y[i] * z[i];
        }

        return r;
      }

      template <typename ValueType_>
      ValueType_ TripleDotProduct::value_blocked_generic(const ValueType_ * const x, const ValueType_ * const y, const ValueType_ * const z, const Index size)
      {
        ValueType_ r(0);

        if (x == y)
        {

          for(Index i(0); i < size; ++i)
          {
            for(int j(0); j < ValueType_::n; ++j)
            {
              r[j] += x[i][j] * x[i][j] * z[i][j];
            }
          }
        }
        else if (x == z)
        {

          for (Index i(0) ; i < size ; ++i)
          {
            for(int j(0); j < ValueType_::n; ++j)
            {
              r[j] += x[i][j] * x[i][j] * y[i][j];
            }
          }
        }
        else if (y == z)
        {

          for (Index i(0) ; i < size ; ++i)
          {
            for(int j(0); j < ValueType_::n; ++j)
            {
              r[j] += x[i][j] * y[i][j] * y[i][j];
            }
          }
        }
        else
        {

          for(Index i(0); i < size; ++i)
          {
            for(int j(0); j < ValueType_::n; ++j)
            {
              r[j] += x[i][j] * y[i][j] * z[i][j];
            }
          }
        }

        return r;
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP

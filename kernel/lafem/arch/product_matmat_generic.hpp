// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

#include <iostream>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      void ProductMatMat::dense_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index columns, const Index inner)
      {
        if (Math::abs(beta) < Math::eps<DT_>())
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < rows ; ++i)
          {
            for (Index j = 0 ; j < columns ; ++j)
            {
              DT_ sum(0.);
              Index xindex(i * inner);
              Index yindex(j);
              for (Index xcol(0) ; xcol < inner ; ++xcol)
              {
                sum  = sum + x[xindex + xcol] * y[yindex + xcol * columns];
              }
              r[i * columns + j] = alpha * sum;
            }
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < rows ; ++i)
          {
            for (Index j = 0 ; j < columns ; ++j)
            {
              DT_ sum(0.);
              Index xindex(i * inner);
              Index yindex(j);
              for (Index xcol(0) ; xcol < inner ; ++xcol)
              {
                sum  = sum + x[xindex + xcol] * y[yindex + xcol * columns];
              }
              r[i * columns + j] = beta * z[i * columns + j] + alpha * sum;
            }
          }
        }
      }

      template <typename DT_, typename IT_>
      void ProductMatMat::dsd_generic(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index /*used_elements*/,
                                         const DT_ * const y, const Index rows,  const Index columns, const Index /*inner*/)
      {
        if (Math::abs(beta) < Math::eps<DT_>())
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < rows ; ++i)
          {
            for (Index j = 0 ; j < columns ; ++j)
            {
              DT_ sum(0.);
              Index xindex = row_ptr[i];
              Index yindex(j);
              for (Index tmp = xindex ; tmp < row_ptr[i+1] ; ++tmp)
              {
                sum  = sum + val[tmp] * y[yindex + col_ind[tmp] * columns];
              }
              r[i * columns + j] = alpha * sum;
            }
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < rows ; ++i)
          {
            for (Index j = 0 ; j < columns ; ++j)
            {
              DT_ sum(0.);
              Index xindex = row_ptr[i];
              Index yindex(j);
              for (Index tmp = xindex ; tmp < row_ptr[i+1] ; ++tmp)
              {
                sum  = sum + val[tmp] * y[yindex + col_ind[tmp] * columns];
              }
              r[i * columns + j] = beta * r[i * columns + j] + alpha * sum;
            }
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_GENERIC_HPP

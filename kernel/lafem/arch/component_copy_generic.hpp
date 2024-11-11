// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_COMPONENT_COPY_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_COMPONENT_COPY_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_COMPONENT_COPY_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      void ComponentCopy::value_generic(DT_ * r, const DT_ * const x, const int stride, const int block, const Index size)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < size; ++i)
        {
          r[i*Index(stride) + Index(block)] = x[i];
        }
      }

      template <typename DT_>
      void ComponentCopy::value_to_generic(const DT_ * const r, DT_ * x, const int stride, const int block, const Index size)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < size; ++i)
        {
          x[i] = r[i*Index(stride) + Index(block)];
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_COMPONENT_COPY_GENERIC_HPP

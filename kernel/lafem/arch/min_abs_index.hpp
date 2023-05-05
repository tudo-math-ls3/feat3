// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_HPP
#define KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct MinAbsIndex
      {
        template <typename DT_>
        static Index value(const DT_ * const x, const Index size)
        {
          return value_generic(x, size);
        }

        static Index value(const float * const x, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, size)
        }

        static Index value(const double * const x, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, size)
        }

        template <typename DT_>
        static Index value_generic(const DT_ * const x, const Index size);

        static Index value_mkl(const float * const x, const Index size);
        static Index value_mkl(const double * const x, const Index size);

        template <typename DT_>
        static Index value_cuda(const DT_ * const x, const Index size);
      };

#ifdef FEAT_EICKT
      extern template Index MinAbsIndex::value_generic(const float * const, const Index);
      extern template Index MinAbsIndex::value_generic(const double * const, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/min_abs_index_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_HPP

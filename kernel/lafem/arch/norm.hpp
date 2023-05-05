// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_NORM_HPP
#define KERNEL_LAFEM_ARCH_NORM_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct Norm2
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const Index size)
        {
          return value_generic(x, size);
        }

#ifdef FEAT_HAVE_HALFMATH
        static Half value(const Half * const x, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_generic, value_generic, x, size)
        }
#endif

        static float value(const float * const x, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, size)
        }

        static double value(const double * const x, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, size)
        }

        template <typename DT_>
        static DT_ value_generic(const DT_ * const x, const Index size);

        static float value_mkl(const float * const x, const Index size);
        static double value_mkl(const double * const x, const Index size);

        template <typename DT_>
        static DT_ value_cuda(const DT_ * const x, const Index size);
      };

#ifdef FEAT_EICKT
      extern template float Norm2::value_generic(const float * const, const Index);
      extern template double Norm2::value_generic(const double * const, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/norm_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_NORM2_HPP

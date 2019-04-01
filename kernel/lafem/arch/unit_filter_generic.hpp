// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void UnitFilter<Mem::Main>::filter_rhs_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
      {
        for(Index i(0); i < ue ; ++i)
        {
          v[sv_indices[i]] = sv_elements[i];
        }
      }

      template <typename DT_, typename IT_>
      void UnitFilter<Mem::Main>::filter_def_generic(DT_ * v, const IT_ * const sv_indices, const Index ue)
      {
        for(Index i(0); i < ue ; ++i)
        {
          v[sv_indices[i]] = DT_(0);
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP

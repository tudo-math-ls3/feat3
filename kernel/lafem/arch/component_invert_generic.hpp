// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_COMPONENT_INVERT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_COMPONENT_INVERT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_COMPONENT_INVERT_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template<typename DT_>
      void ComponentInvert<Mem::Main>::value_generic(DT_* r, const DT_* const x, const DT_ s, const Index size)
      {
        if (r == x)
        {
          for(Index i(0); i < size; ++i)
          {
            r[i] = s / r[i];
          }
        }
        else
        {
          for(Index i(0); i < size; ++i)
          {
            r[i] = s / x[i];
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_COMPONENT_INVERT_HPP

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_NORM2_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  DT_ Norm2Dense::exec_generic_impl(const DT_ * const x, const Index size)
  {
    DT_ r(0);

    FEAT_PRAGMA_OMP(parallel for reduction(+:r))
    for (Index i = 0 ; i < size ; ++i)
    {
      r += x[i] * x[i];
    }

    return Math::sqrt(r);
  }
} // namespace FEAT::LAFEM::Arch

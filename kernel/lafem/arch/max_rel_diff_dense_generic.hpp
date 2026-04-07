// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MAX_REL_DIFF_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  DT_ MaxRelDiffDense::exec_generic_impl(const DT_ * const x, const DT_ * const y, const Index size)
  {
    DT_ r(0);
    bool has_nan = false;
    bool has_inf = false;
    for (Index i(0) ; i < size ; ++i)
    {
      has_nan = has_nan || Math::isnan(x[i]) || Math::isnan(y[i]);
      has_inf = has_inf || Math::isinf(x[i]) || Math::isinf(y[i]);
      r = Math::max(r, Math::abs(x[i] - y[i]) / (Math::abs(x[i]) + Math::abs(y[i]) + DT_(1)));
    }

    // if we found a NaN or an infinity, we'll return a NaN here
    return (has_nan || has_inf ? Math::nan<DT_>() : r);
  }
} // namespace FEAT::LAFEM::Arch

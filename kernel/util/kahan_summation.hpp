// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_KAHAN_SUMMATION_HPP
#define KERNEL_UTIL_KAHAN_SUMMATION_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>

namespace FEAT
{
  /**
   * This struct holds the intermediate summation result and the calculation error from the previous sum operation
   * See https://en.wikipedia.org/wiki/Kahan_summation_algorithm and
   * https://stackoverflow.com/questions/10330002/sum-of-small-double-numbers-c/10330857#10330857
   * for details.
  */
  struct KahanAccumulation
  {
    double sum;
    double correction;

    KahanAccumulation() :
      sum(0.0),
      correction(0.0)
    {
    }
  };

  /**
   * This struct executes the actual Kahan Summation step.
   * See https://en.wikipedia.org/wiki/Kahan_summation_algorithm and
   * https://stackoverflow.com/questions/10330002/sum-of-small-double-numbers-c/10330857#10330857
   * for details.
   */
  KahanAccumulation KahanSum(KahanAccumulation accumulation, double value);

} // namespace FEAT

#endif // KERNEL_UTIL_KAHAN_SUMMATION_HPP

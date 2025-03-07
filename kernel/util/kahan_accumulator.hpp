// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

namespace FEAT
{
  /**
   * \brief Kahan Summation accumulator class template
   *
   * This class implements the Kahan summation algorithm for a given floating point data type,
   * see e.g. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
   *
   * \author Peter Zajac
   */
  template<typename DT_>
  class KahanAccumulator
  {
  public:
    /// the current sum value
    DT_ value;
    /// the current correction
    DT_ correction;

    /// default constructor
    KahanAccumulator() :
      value(DT_(0)),
      correction(DT_(0))
    {
    }

    /// constructor
    explicit KahanAccumulator(DT_ val) :
      value(val),
      correction(DT_(0))
    {
    }

    /// clears the accumulator, i.e. sets all values to zero
    void clear()
    {
      value = correction = DT_(0);
    }

    /// updates the summed accumulator with a new summand
    KahanAccumulator& operator+=(DT_ val)
    {
      DT_ y = val - this->correction;
      DT_ t = this->value + y;
      this->correction = (t - this->value) - y;
      this->value = t;
      return *this;
    }

    /// returns the value of the sum
    operator DT_() const
    {
      return value;
    }
  }; // class KahanAccumulator<...>
} // namespace FEAT

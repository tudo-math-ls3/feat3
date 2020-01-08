// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_META_MATH_HPP
#define KERNEL_UTIL_META_MATH_HPP 1

#include <kernel/util/assertion.hpp>

namespace FEAT
{
  /**
   * \brief Template Meta-Program Math namespace
   *
   * This namespace encapsulated mathematical functions, which are written as template meta-programs
   * so that they may be evaluated at compile-time.
   */
  namespace MetaMath
  {
    /**
     * \brief Factorial template meta-program
     *
     * This meta-program calculates the coefficient
     * \f[ f(n,m) := \prod_{k=m}^n k \f]
     * at compile-time. The factorial <em>n!</em> is given by <em>f(n,0)</em>.
     *
     * \author Peter Zajac
     */
    template<int n_, int m_ = 0>
    struct Factorial
    {
      // Ensure that the parameter m_ is non-negative. The case n_ = 0 is handled by specialisation below.
      static_assert(m_ >= 0, "parameter m_ must be non-negative");
      // Ensure that n_ >= m_. The case m_ = n_ is handled by partial specialisation below.
      static_assert(n_ > m_, "parameter m_ must not be greater than parameter n_");

      /// value of the factorial
      static constexpr int value = n_ * Factorial<n_ - 1, m_>::value;
    }; // struct Factorial<n_,m_>

    /// \cond internal
    // partial specialisation for n_ = m_
    template<int n_>
    struct Factorial<n_, n_>
    {
      static_assert(n_ >= 0, "parameters n_ = m_ must be non-negative");
      static constexpr int value = n_;
    }; // struct Factorial<n_,n_>

    // explicit specialisation for n_ = m_ = 0; the partial specialisation above for n_ = m_ would return 0
    template<>
    struct Factorial<0, 0>
    {
      static constexpr int value = 1;
    }; // struct Factorial<0,0>
    /// \endcond

    /**
     * \brief Binomial template meta-program
     *
     * This meta-program calculates the <em>binomial coefficient</em>
     * \f[ {n \choose k} := \frac{n!}{k!(n-k)!} \f]
     * at compile-time.
     *
     * \author Peter Zajac
     */
    template<int n_, int k_>
    struct Binomial
    {
      // note: the valid cases k_ = 0 and k_ = n_ are specialised below
      static_assert(k_ > 0, "parameter k_ must be non-negative");
      static_assert(n_ > k_, "parameter k_ must not be greater than parameter n_");

      // (n,k) = (n-1,k-1) + (n-1, k)
      static constexpr int value = Binomial<n_-1,k_-1>::value + Binomial<n_-1,k_>::value;
    }; // struct Binomial<n_,k_>

    /// \cond internal
    // partial specialisation for k_ = n_
    template<int n_>
    struct Binomial<n_,n_>
    {
      static_assert(n_ > 0, "parameter n_ must be non-negative");
      static constexpr int value = 1;
    };

    // partial specialisation for k_ = 0
    template<int n_>
    struct Binomial<n_, 0>
    {
      static_assert(n_ > 0, "parameter n_ must be non-negative");
      static constexpr int value = 1;
    };

    // explicit specialisation for k_ = n_ = 0; this is only needed for disambiguation
    template<>
    struct Binomial<0, 0>
    {
      static constexpr int value = 1;
    };
    /// \endcond
  } // namespace MetaMath
} // namespace FEAT

#endif // KERNEL_UTIL_META_MATH_HPP

#pragma once
#ifndef KERNEL_UTIL_FACTORIAL_HPP
#define KERNEL_UTIL_FACTORIAL_HPP 1

#include <kernel/base_header.hpp>

namespace FEAST
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
    enum
    {
      /// value of the factorial
      value = n_ * Factorial<n_ - 1, m_>::value
    };
  }; // struct Factorial<n_,m_>

  /// \cond internal
  // partial specialisation for n_ = m_
  template<int n_>
  struct Factorial<n_, n_>
  {
    static_assert(n_ >= 0, "parameters n_ = m_ must be non-negative");
    enum
    {
      value = n_
    };
  }; // struct Factorial<n_,n_>

  // explicit specialisation for n_ = m_ = 0; the partial specialisation above for n_ = m_ would return 0
  template<>
  struct Factorial<0, 0>
  {
    enum
    {
      value = 1
    };
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
    enum
    {
      // (n,k) = (n-1,k-1) + (n-1, k)
      value = Binomial<n_-1,k_-1>::value + Binomial<n_-1,k_>::value
    };
  }; // struct Binomial<n_,k_>

  /// \cond internal
  // partial specialisation for k_ = n_
  template<int n_>
  struct Binomial<n_,n_>
  {
    static_assert(n_ > 0, "parameter n_ must be non-negative");
    enum
    {
      value = 1
    };
  };

  // partial specialisation for k_ = 0
  template<int n_>
  struct Binomial<n_, 0>
  {
    static_assert(n_ > 0, "parameter n_ must be non-negative");
    enum
    {
      value = 1
    };
  };

  // explicit specialisation for k_ = n_ = 0; this is only needed for disambiguation
  template<>
  struct Binomial<0, 0>
  {
    enum
    {
      value = 1
    };
  };
  /// \endcond

} // namespace FEAST

#endif // KERNEL_UTIL_FACTORIAL_HPP

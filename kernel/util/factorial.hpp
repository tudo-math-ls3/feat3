#pragma once
#ifndef KERNEL_UTIL_FACTORIAL_HPP
#define KERNEL_UTIL_FACTORIAL_HPP 1

#include <kernel/util/assertion.hpp>

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

    /// dummy enumeration
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
   * \brief Calculates the (partial) factorial.
   *
   * This function calculates the coefficient
   * \f[ f(n,m) := \prod_{k=m}^n k \f]
   * The factorial <em>n!</em> is given by <em>f(n,0)</em>.
   *
   * \param[in] n, m
   * The coefficients for the (partial) factorial. It must hold \e m <= \e n <= 20.
   *
   * \warning
   * The parameter \e n must not be greater than 20, since for \e n > 20 the resulting factorial \e n! is greater
   * than 2^64, thus leading to integer overflow!
   *
   * \returns
   * The partial factorial from \e m to \e n.
   */
  inline unsigned long long factorial(
    unsigned long long n,
    unsigned long long m = 0ull)
  {
    CONTEXT("factorial()");

    // n > 20 ==> n! > 2^64 ==> integer overflow!
    ASSERT(n <= 20, "parameter n must be less than 21");
    ASSERT(m <= n, "parameter m must not be greater than n");

    if(n <= 1ull)
      return 1ull;

    // calculate the factorial
    unsigned long long k = 1ull;
    for(m = std::max(1ull, m); m <= n; ++m)
    {
      k *= m;
    }
    return k;
  }

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

    /// dummy enumeration
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

  /**
   * \brief Calculates the binomial coefficient.
   *
   * This function calculates the binomial coefficient
   * \f[ {n \choose k} := \frac{n!}{k!(n-k)!} \f]
   *
   * \param[in] n,k
   * The parameters for the binomial coefficient.
   *
   * \returns
   * The binomial coefficient <em>n choose k</em>.
   */
  inline unsigned long long binomial(
    unsigned long long n,
    unsigned long long k)
  {
    if(k > n)
    {
      return 0; // by definition
    }
    else if((k == 0ull) || (k == n))
    {
      return 1; // by definition
    }

    // exploit symmetry: (n \choose k) = (n \choose n-k)
    k = std::min(k, n-k);

    // use multiplicative formula: (n \choose k+1) = (n - k) * (n \choose k) / (k + 1)
    unsigned long long m = n;
    for(unsigned long long i(1); i < k; ++i)
    {
      m *= (n - i);
      m /= (i + 1);
    }

    return m;
  }
} // namespace FEAST

#endif // KERNEL_UTIL_FACTORIAL_HPP

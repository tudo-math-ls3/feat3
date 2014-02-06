#pragma once
#ifndef KERNEL_UITL_MATH_HPP
#define KERNEL_UITL_MATH_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <cmath>
#include <limits>

namespace FEAST
{
  /**
   * \brief FEAST Math namespace
   *
   * This namespace encapsulates common mathematical functions used throughout the FEAST kernel.
   * This namespace contains all math functions from the standard <c>cmath</c> header.
   */
  namespace Math
  {
    // include C++ overloads of C89 math functions
    using std::abs;
    using std::acos;
    using std::asin;
    using std::atan;
    using std::atan2;
    using std::ceil;
    using std::cos;
    using std::cosh;
    using std::exp;
    using std::floor;
    using std::fmod;
    using std::log;
    using std::log10;
    using std::modf;
    using std::pow;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    using std::tan;
    using std::tanh;

    /**
     * \brief Returns the square of a value.
     *
     * \param[in] x
     * The value to be squared.
     *
     * Returns x*x
     */
    template<typename T_>
    inline T_ sqr(T_ x)
    {
      return x * x;
    }

    /**
     * \brief Returns the cube of a value.
     *
     * \param[in] x
     * The value to be cubed.
     *
     * Returns x*x*x
     */
    template<typename T_>
    inline T_ cub(T_ x)
    {
      return x * x * x;
    }

    /**
     * \brief Returns the minimum of two values.
     *
     * \param[in] a,b
     * The two values whose minimum is to be returned.
     *
     * \returns The minimum of \p a and \p b.
     */
    template<typename T_>
    inline T_ min(T_ a, T_ b)
    {
      return (a < b ? a : b);
    }

    /**
     * \brief Returns the maximum of two values.
     *
     * \param[in] a,b
     * The two values whose maximum is to be returned.
     *
     * \returns The maximum of \p a and \p b.
     */
    template<typename T_>
    inline T_ max(T_ a, T_ b)
    {
      return (a < b ? b : a);
    }

    /**
     * \brief Clamps a value to a range.
     *
     * \param[in] x
     * The value to be clamped.
     *
     * \param[in] a,b
     * The clamp boundaries.
     *
     * \returns
     * max(a, min(x, b))
     */
    template<typename T_>
    inline T_ clamp(T_ x, T_ a, T_ b)
    {
      return max(a, min(x, b));
    }

    /**
     * \brief Calculates the (partial) factorial.
     *
     * This function calculates the coefficient
     * \f[ f(n,m) := \prod_{k=m}^n k \f]
     * The factorial <em>n!</em> is given by <em>f(n,0)</em>.
     *
     * \attention
     * Please keep in mind that the factorial becomes very huge very fast, i.e. for
     * integral data types, one should honor the following boundaries to avoid overflow:
     *  -  8-bit: n <= 5
     *  - 16-bit: n <= 7 (signed); n <= 8 (unsigned)
     *  - 32-bit: n <= 12
     *  - 64-bit: n <= 20
     *
     * \param[in] n, m
     * The coefficients for the (partial) factorial.
     *
     * \returns
     * The partial factorial from \e m to \e n.
     *
     * \author Peter Zajac
     */
    template<typename T_>
    inline T_ factorial(T_ n, T_ m = T_(0))
    {
      // calculate the factorial
      T_ k(T_(1));
      for(m = max(T_(1), m); m <= n; ++m)
      {
        k *= m;
      }
      return k;
    }

    /**
     * \brief Calculates the binomial coefficient.
     *
     * This function calculates the binomial coefficient
     * \f[ {n \choose k} := \frac{n!}{k!(n-k)!} \f]
     *
     * \attention
     * This function works only for integral data types;
     * it will give incorrect results for floating point data types!
     *
     * \param[in] n,k
     * The parameters for the binomial coefficient.
     *
     * \returns
     * The binomial coefficient <em>n choose k</em>.
     *
     * \author Peter Zajac
     */
    template<typename T_>
    inline T_ binomial(T_ n, T_ k)
    {
      if(k > n)
      {
        return T_(0); // by definition
      }
      else if((k <= T_(0)) || (k == n))
      {
        return T_(1); // by definition
      }

      // exploit symmetry: (n \choose k) = (n \choose n-k)
      k = min(k, n-k);

      // use multiplicative formula: (n \choose k+1) = (n - k) * (n \choose k) / (k + 1)
      T_ m = n;
      for(T_ i(1); i < k; ++i)
      {
        m *= (n - i);
        m /= (i + 1);
      }

      return m;
    }

    /**
     * \brief Math Limits class template
     *
     * This class is an extended version of the <c>std::numeric_limits</c> class template, from which
     * it derives.
     */
    template<typename T_>
    class Limits :
      public std::numeric_limits<T_>
    {
    public:
      /// Returns the mathematical constant pi = 3.14...
      static T_ pi()
      {
        return T_(2) * acos(T_(0));
      }
    }; // class Limits<...>

    /**
     * \brief Returns the machine precision constant for a floating-point data type.
     */
    template<typename DataType_>
    static inline DataType_ eps()
    {
      return Limits<DataType_>::epsilon();
    }
  } // namespace Math
} // namespace FEAST

#endif // KERNEL_UITL_MATH_HPP

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_MATH_HPP
#define KERNEL_UTIL_MATH_HPP 1

// includes, FEAT
#include <kernel/util/type_traits.hpp>

// includes, system
#include <cmath>
#include <limits>

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
extern "C"
{
#    include <quadmath.h>
}
#  define CAT_(x,y) x##y
#  define CAT(x,y) CAT_(x,y)
#  define WRAP_QUAD_MATH1(func) \
    inline __float128 func(__float128 x) {return ::CAT(func,q)(x);}
#  define WRAP_QUAD_MATH2(func) \
    inline __float128 func(__float128 x, __float128 y) {return ::CAT(func,q)(x, y);}
#  define WRAP_QUAD_MATH2PTR(func) \
    inline __float128 func(__float128 x, __float128* y) {return ::CAT(func,q)(x, y);}
#else
#  define WRAP_QUAD_MATH1(func)
#  define WRAP_QUAD_MATH2(func)
#  define WRAP_QUAD_MATH2PTR(func)
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__

    // single argument function wrapper
#define WRAP_STD_MATH1(func) \
    inline float func(float x) {return std::func(x);} \
    inline double func(double x) {return std::func(x);} \
    inline long double func(long double x) {return std::func(x);}

    // double argument function wrapper
#define WRAP_STD_MATH2(func) \
    inline float func(float x, float y) {return std::func(x,y);} \
    inline double func(double x, double y) {return std::func(x,y);} \
    inline long double func(long double x, long double y) {return std::func(x,y);}

    // double argument function wrapper, second argument is pointer
#define WRAP_STD_MATH2PTR(func) \
    inline float func(float x, float* y) {return std::func(x,y);} \
    inline double func(double x, double* y) {return std::func(x,y);} \
    inline long double func(long double x, long double* y) {return std::func(x,y);}

    // single argument function wrapper, bool return type
#define WRAP_STD_MATH1BRET(func) \
    inline bool func(float x) {return std::func(x);} \
    inline bool func(double x) {return std::func(x);} \
    inline bool func(long double x) {return std::func(x);}

namespace FEAT
{
  /**
   * \brief FEAT Math namespace
   *
   * This namespace encapsulates common mathematical functions used throughout the FEAT kernel.
   * This namespace contains all math functions from the standard <c>cmath</c> header.
   */
  namespace Math
  {
    // include C++ overloads of C89 math functions
    WRAP_STD_MATH1(ceil)
    WRAP_STD_MATH1(floor)
    WRAP_STD_MATH2(fmod)
    WRAP_STD_MATH2PTR(modf)

    // wrap quadmath functions
    WRAP_QUAD_MATH1(ceil)
    WRAP_QUAD_MATH1(floor)
    WRAP_QUAD_MATH2(fmod)
    WRAP_QUAD_MATH2PTR(modf)

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
     * \brief Updates the minimum and maximum.
     *
     * This functions performs:
     * - a := min(x, a)
     * - b := max(x, b)
     *
     * \param[in] x
     * The new value for the minimum/maximum update.
     *
     * \param[in,out] a
     * The minimum that is to be updated.
     *
     * \param[in,out] b
     * The maximum that is to be updated.
     */
    template<typename T_>
    inline void minimax(T_ x, T_& a, T_& b)
    {
      if(x < a)
        a = x;
      if(b < x)
        b = x;
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
     * \brief Computes the integral base-10 logarithm of an integer, i.e. its number of non-zero decimal digits.
     *
     * \param[in] x
     * The number whose digit count is to be computed.
     *
     * \returns
     * The number of non-zero decimal digits of \p x.
     */
    template<typename T_>
    inline T_ ilog10(T_ x)
    {
      static_assert(Type::Traits<T_>::is_int, "ilog10 can only be applied to integral types");
      T_ i(0);
      while(x != T_(0))
      {
        ++i;
        x /= T_(10);
      }
      return i;
    }

    /**
     * \brief Returns the sign of a value.
     *
     * \param[in] x
     * The value whose sign is to be returned.
     */
    template<typename T_>
    inline T_ signum(T_ x)
    {
      return (x < T_(0) ? T_(-1) : (x > T_(0) ? T_(1) : T_(0)));
    }

    /**
     * \brief Returns the status of the sign bit.
     *
     * \param[in] x
     * The value whose signbit status is to be returned.
     */
    template<typename T_>
    inline bool signbit(T_ x)
    {
      return x < T_(0);
    }

    /**
     * \brief Returns the absolute value.
     *
     * \param[in] x The value to get the absolute value from.
     *
     * \returns abs(x)
     */
    template<typename T_>
    inline T_ abs(T_ x)
    {
      return (x < T_(0.) ? -x : x);
    }

    // wrap std::abs
    WRAP_STD_MATH1(abs)
#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    inline __float128 abs(__float128 x) {return ::fabsq(x);}
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__

    /**
     * \brief Returns the square-root of a value.
     *
     * \param[in] x The value to calculate the square-root from.
     *
     * \returns sqrt(x)
     *
     * \compilerhack disable "warning C4723: potential divide by 0" for MSVC compiler
     */
#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(push)
#pragma warning(disable:4723)
#endif
    template<typename T_>
    inline T_ sqrt(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "sqrt can only be applied to floating point types");

      if(x <= T_(0))
        return T_(0);

      // use Newton iteration: y_{k+1} := y_k/2 * (3 - (y_k^2)/x)
      // we choose y_0 = min(1,x); this ensures that the sequence y_k is monotonically increasing
      // if y_{k+1} is not greater than y_k, we return y_k
      const T_ z = T_(1) / x;
      T_ y(Math::min(T_(1), x));
      T_ yn(y);
      do
      {
        y = yn;
        yn = T_(0.5)*y * (T_(3) - (y*y*z));
      } while(yn > y);

      return y;
    }
#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(pop)
#endif

    // wrap std::sqrt
    WRAP_STD_MATH1(sqrt)
    WRAP_QUAD_MATH1(sqrt)

    /**
     * \brief Returns the sine of a value.
     *
     * \param[in] x The value to calculate the sine from.
     *
     * \returns sin(x)
     */
    template<typename T_>
    inline T_ sin(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "sin can only be applied to floating point types");

      // use the exponential sum formula:
      //           infty        x^(2*n+1)
      // sin(x) :=  sum (-1)^n -----------
      //            n=0         (2*n+1)!
      T_ y(x), yl(x+T_(1)), z(x);
      T_ fn(1.0);
      int n(1);
      do
      {
        // update 1/(2*n+1)!
        fn /= T_(++n);
        fn /= T_(++n);
        yl = y;
        y += T_(1 - int(n&2)) * (z *= x*x) * T_(fn);
      } while(yl != y);

      return y;
    }

    // wrap std::sin
    WRAP_STD_MATH1(sin)
    WRAP_QUAD_MATH1(sin)

    /**
     * \brief Returns the cosine of a value.
     *
     * \param[in] x The value to calculate the cosine from.
     *
     * \returns cos(x)
     */
    template<typename T_>
    inline T_ cos(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "cos can only be applied to floating point types");

      // use the exponential sum formula:
      //           infty        x^(2*n)
      // cos(x) :=  sum (-1)^n ---------
      //            n=0         (2*n)!
      T_ y(T_(1)), yl(T_(0)), z(T_(1));
      T_ fn(1.0);
      int n(0);
      do
      {
        // update 1/(2*n)!
        fn /= T_(++n);
        fn /= T_(++n);
        yl = y;
        y += T_(1 - int(n&2)) * (z *= x*x) * T_(fn);
      } while(yl != y);

      return y;
    }

    // wrap std::cos
    WRAP_STD_MATH1(cos)
    WRAP_QUAD_MATH1(cos)

    /**
     * \brief Returns the tangent of a value.
     *
     * \param[in] x The value to calculate the tangent from.
     *
     * \returns tan(x)
     */
    template<typename T_>
    inline T_ tan(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "tan can only be applied to floating point types");

      return sin(x) / cos(x);
    }

    // wrap std::tan
    WRAP_STD_MATH1(tan)
    WRAP_QUAD_MATH1(tan)

    /**
     * \brief Returns the hyperbolic sine of a value.
     *
     * \param[in] x The value to calculate the hyperbolic sine from.
     *
     * \returns sinh(x)
     */
    template<typename T_>
    inline T_ sinh(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "sinh can only be applied to floating point types");

      // use the exponential sum formula:
      //            infty  x^(2*n+1)
      // sinh(x) :=  sum  -----------
      //             n=0   (2*n+1)!
      T_ y(x), yl(x+T_(1)), z(x);
      T_ fn(1.0);
      int n(1);
      do
      {
        // update 1/(2*n+1)!
        fn /= T_(++n);
        fn /= T_(++n);
        yl = y;
        y += (z *= x*x) * T_(fn);
      } while(yl != y);

      return y;
    }

    // wrap std::sinh
    WRAP_STD_MATH1(sinh)
    WRAP_QUAD_MATH1(sinh)

    /**
     * \brief Returns the hyperbolic cosine of a value.
     *
     * \param[in] x The value to calculate the hyperbolic cosine from.
     *
     * \returns cosh(x)
     */
    template<typename T_>
    inline T_ cosh(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "cosh can only be applied to floating point types");

      // use the exponential sum formula:
      //            infty   x^(2*n)
      // cosh(x) :=  sum   ---------
      //             n=0     (2*n)!
      T_ y(T_(1)), yl(T_(0)), z(T_(1));
      T_ fn(1.0);
      int n(0);
      do
      {
        // update 1/(2*n)!
        fn /= T_(++n);
        fn /= T_(++n);
        yl = y;
        y += (z *= x*x) * T_(fn);
      } while(yl != y);

      return y;
    }

    // wrap std::cosh
    WRAP_STD_MATH1(cosh)
    WRAP_QUAD_MATH1(cosh)

    /**
     * \brief Returns the hyperbolic tangent of a value.
     *
     * \param[in] x The value to calculate the hyperbolic tangent from.
     *
     * \returns tanh(x)
     */
    template<typename T_>
    inline T_ tanh(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "tanh can only be applied to floating point types");

      return sinh(x) / cosh(x);
    }

    // wrap std::tanh
    WRAP_STD_MATH1(tanh)
    WRAP_QUAD_MATH1(tanh)

    /**
     * \brief Returns the exponential of a value.
     *
     * \param[in] x The value to calculate the exponential from.
     *
     * \returns exp(x)
     */
    template<typename T_>
    inline T_ exp(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "exp can only be applied to floating point types");

      T_ y(T_(1)), yl(T_(0)), z(y);
      int n(0);
      do
      {
        yl = y;
        y += ((z *= x) /= T_(++n));
        // Note about the stopping criterion:
        // For x > 0, the sequence y_k must be strictly increasing.
        // For x < 0, the sequence y_k must be alternating.
        // And encode this into the most beautiful crypto-expression C++ has to offer ^_^
      } while(x > T_(0) ? yl < y : n & 1 ? y < yl : yl < y);
      return yl;
    }

    // wrap std::exp
    WRAP_STD_MATH1(exp)
    WRAP_QUAD_MATH1(exp)

    /**
     * \brief Returns the natural logarithm of a value.
     *
     * \param[in] x The value to calculate the natural logarithm from.
     *
     * \returns log(x)
     */
    template<typename T_>
    inline T_ log(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "log can only be applied to floating point types");

      // use Newton iteration: y_{k+1} = y_k + 2*(x - exp(y_k))/(x + exp(y_k))
      T_ y(T_(0)), yl(T_(0));
      do
      {
        yl = y;
        T_ ey(Math::exp(y));
        y += T_(2) * (x - ey) / (x + ey);
        // Note about the stopping criterion:
        // For x > 1, the sequence y_k must be strictly increasing.
        // For x < 1, the sequence y_k must be strictly decreasing.
        // Again, encode this into one beautiful expression
      } while(x < T_(1) ? y < yl : yl < y);
      return yl;
    }

    // wrap std::log
    WRAP_STD_MATH1(log)
    WRAP_QUAD_MATH1(log)

    /**
     * \brief Returns the logarithm to the base 10 of a value.
     *
     * \param[in] x The value to calculate the 10-logarithm from.
     *
     * \returns log10(x)
     */
    template<typename T_>
    inline T_ log10(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "log10 can only be applied to floating point types");

      return log(x) / log(T_(10));
    }

    // wrap std::log10
    WRAP_STD_MATH1(log10)
    WRAP_QUAD_MATH1(log10)

    /**
     * \brief Returns x raised to the power of y.
     *
     * \param[in] x The base value. Must be positive.
     * \param[in] y The exponent value.
     *
     * \returns x^y
     */
    template<typename T_>
    inline T_ pow(T_ x, T_ y)
    {
      static_assert(Type::Traits<T_>::is_float, "pow can only be applied to floating point types");

      return Math::exp(y * Math::log(x));
    }

    // wrap std::pow
    WRAP_STD_MATH2(pow)
    WRAP_QUAD_MATH2(pow)

    /**
     * \brief Returns the arctangent of a value.
     *
     * \param[in] x The value to calculate the arctangent from.
     *
     * \returns atan(x)
     */
    template<typename T_>
    inline T_ atan(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "atan can only be applied to floating point types");

      // the exponential sum converges only for |x| < 1, but we can reduce any |x| >= 1 by
      // atan(x) = 2*atan( x / (1 + sqrt(1 + x^2)))
      int k(0);
      for(; Math::abs(x) >= T_(1); ++k)
        x /= (T_(1) + Math::sqrt(T_(1) + x*x));

      // use the exponential sum formula:
      //            infty        x^(2*n+1)
      // atan(x) :=  sum (-1)^n -----------
      //             n=0         (2*n+1)
      T_ y(x), yl(x+T_(1)), z(x);
      int n(1);
      do
      {
        yl = y;
        T_ t((z *= x*x) / T_(n += 2));
        y += T_(1 - int(n&2)) * t;
      } while(yl != y);

      return T_(1<<k) * y;
    }

    // wrap std::atan
    WRAP_STD_MATH1(atan)
    WRAP_QUAD_MATH1(atan)

    /**
     * \brief Returns the arctangent of y/x.
     *
     * \param[in] y The nominator of the value to calculate the arctangent from.
     * \param[in] x The denominator of the value to calculate the arctangent from.
     *
     * \returns atan2(y,x) = atan(y/x)
     */
    template<typename T_>
    inline T_ atan2(T_ y, T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "atan2 can only be applied to floating point types");

      // see http://en.wikipedia.org/wiki/Atan2#Variations_and_notation
      return T_(2) * atan((Math::sqrt(x*x + y*y) - x) / y);
    }

    // wrap std::atan2
    WRAP_STD_MATH2(atan2)
    WRAP_QUAD_MATH2(atan2)

    /**
     * \brief Returns the mathematical constant pi = 3.1415...
     */
    template<typename T_>
    inline T_ pi()
    {
      static_assert(Type::Traits<T_>::is_float, "pi can only be applied to floating point types");

      // use the Bailey-Borwein-Plouffe formula:
      //       infty   1     (   4      2      1      1  )
      // pi :=  sum  ----  * ( ---- - ---- - ---- - ---- )
      //        k=0  16^k    ( 8k+1   8k+4   8k+5   8k+6 )
      T_ y(T_(0)), yl(T_(0));
      int k(0);
      const T_ z(T_(1) / T_(16));
      T_ t(T_(1));
      do
      {
        yl = y;
        y += t * (T_(4)/T_(8*k+1) - T_(2)/T_(8*k+4) - T_(1)/T_(8*k+5) - T_(1)/T_(8*k+6));
        t *= z;
        ++k;
      } while(yl != y);

      return y;
    }

    /// \cond internal
    template<>
    inline float pi<float>()
    {
      return 4.0f * std::atan(1.0f);
    }

    template<>
    inline double pi<double>()
    {
      return 4.0 * std::atan(1.0);
    }

    template<>
    inline long double pi<long double>()
    {
      return 4.0l * std::atan(1.0l);
    }
    /// \endcond

    /**
     * \brief Returns the machine precision constant for a floating-point data type.
     */
    template<typename T_>
    inline T_ eps()
    {
      static_assert(Type::Traits<T_>::is_float, "eps can only be applied to floating point types");

      const T_ z(T_(1));
      const T_ t(T_(0.5));
      T_ y(t), yl(t);
      do
      {
        yl = y;
        y *= t;
      } while(z < T_(z+y));
      return yl;
    }

    /// \cond internal
    template<>
    inline float eps<float>()
    {
      return std::numeric_limits<float>::epsilon();
    }
    template<>
    inline double eps<double>()
    {
      return std::numeric_limits<double>::epsilon();
    }
    template<>
    inline long double eps<long double>()
    {
      return std::numeric_limits<long double>::epsilon();
    }

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    template<>
    inline __float128 eps<__float128>()
    {
      return FLT128_EPSILON;
    }
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__
    /// \endcond

    /**
     * \brief Returns the maximum positive finite (full-precision) value for a data type.
     *
     * The return value of this function coincides with <c>std::numeric_limits<T_>::max()</c>.
     */
    template<typename T_>
    inline T_ huge()
    {
      return std::numeric_limits<T_>::max();
    }

    /**
     * \brief Returns the minimum positive finite (full-precision) value for a data type.
     *
     * The return value of this function coincides with <c>std::numeric_limits<T_>::min()</c>.
     */
    template<typename T_>
    inline T_ tiny()
    {
      return std::numeric_limits<T_>::min();
    }

    /// \cond internal
#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    template<>
    inline __float128 huge<__float128>()
    {
      return FLT128_MAX;
    }

    template<>
    inline __float128 tiny<__float128>()
    {
      return FLT128_MIN;
    }
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__
    /// \endcond

    /**
     * \brief Returns a quiet Not-A-Number (NaN)
     *
     * \note The generic implementation simply returns 0/0, which should
     * result in a NaN for any IEEE-754 conforming implementation.
     *
     * \returns NaN
     */
    template<typename T_>
    inline T_ nan()
    {
      // divide 0 by 0, which hopefully yields NaN
      return T_(0) / T_(0);
    }

    /// \cond internal
    template<>
    inline float nan<float>()
    {
      return std::nanf("");
    }

    template<>
    inline double nan<double>()
    {
      return std::nan("");
    }

    template<>
    inline long double nan<long double>()
    {
      return std::nanl("");
    }

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    template<>
    inline __float128 nan<__float128>()
    {
      return ::nanq("");
    }
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_HALFMATH) && !defined(__CUDACC__)
    template<>
    inline half_float::half nan<half_float::half>()
    {
      return half_float::nanh("");
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_FLOATX) && !defined(__CUDACC__)
    /*template<int exp_bits_, int sig_bits_, typename Backend_>
    inline flx::floatx<exp_bits_, sig_bits_, Backend_> nan<flx::floatx<exp_bits_, sig_bits_, Backend_>>()
    {
      // FloatX doesn't offer its own nan implementation,
      // so create a backend type NaN and convert it to FloatX
      return flx::floatx<exp_bits_, sig_bits_, Backend_>(Math::nan<Backend_>());
    }*/
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__
    /// \endcond

    /**
     * \brief Returns the arcsine of a value.
     *
     * \param[in] x The value to calculate the arcsine from.
     *
     * \returns asin(x)
     */
    template<typename T_>
    inline T_ asin(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "asin can only be applied to floating point types");

      return signum(x) * Math::atan(Math::sqrt((x*x) / (T_(1) - x*x)));
    }

    WRAP_STD_MATH1(asin)
    WRAP_QUAD_MATH1(asin)

    /**
     * \brief Returns the arccosine of a value.
     *
     * \param[in] x The value to calculate the arccosine from.
     *
     * \returns acos(x)
     */
    template<typename T_>
    inline T_ acos(T_ x)
    {
      static_assert(Type::Traits<T_>::is_float, "acos can only be applied to floating point types");

      return T_(0.5) * pi<T_>() - Math::asin(x);
    }

    WRAP_STD_MATH1(acos)
    WRAP_QUAD_MATH1(acos)

    /**
     * \brief Checks whether a value is finite.
     *
     * A value is finite if it is neither NaN nor inifinity.
     *
     * \note There exists no generic implementation for this function.
     *
     * \param[in] x The value to be checked for finiteness.
     *
     * \returns \c true if \p x is finite, otherwise \c false.
     */
    template<typename T_>
    inline bool isfinite(T_ x);

    WRAP_STD_MATH1BRET(isfinite)

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    inline bool isfinite(__float128 x)
    {
      // https://chromium.googlesource.com/native_client/nacl-gcc/+/ng/master/libquadmath/math/finiteq.c
      return (::finiteq(x) != 0);
    }
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_HALFMATH) && !defined(__CUDACC__)
    inline bool isfinite(half_float::half x)
    {
      return half_float::isfinite(x);
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_FLOATX) && !defined(__CUDACC__)
    template<int exp_bits_, int sig_bits_, typename Backend_>
    inline bool isfinite(const flx::floatx<exp_bits_, sig_bits_, Backend_>& x)
    {
      // FloatX doesn't offer its own isfinite implementation,
      // so test its backend implementation instead
      return isfinite(static_cast<Backend_>(x));
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

    /**
     * \brief Checks whether a value is infinite.
     *
     * \note There exists no generic implementation for this function.
     *
     * \param[in] x The value to be checked for infinity.
     *
     * \returns \c true if \p x is +/-infinity, otherwise \c false.
     */
    template<typename T_>
    inline bool isinf(T_ x);

    WRAP_STD_MATH1BRET(isinf)

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    inline bool isinf(__float128 x)
    {
      // https://chromium.googlesource.com/native_client/nacl-gcc/+/ng/master/libquadmath/math/isinfq.c
      return (::isinfq(x) != 0);
    }
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_HALFMATH) && !defined(__CUDACC__)
    inline bool isinf(half_float::half x)
    {
      return half_float::isinf(x);
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_FLOATX) && !defined(__CUDACC__)
    template<int exp_bits_, int sig_bits_, typename Backend_>
    inline bool isinf(const flx::floatx<exp_bits_, sig_bits_, Backend_>& x)
    {
      // FloatX doesn't offer its own isinf implementation,
      // so test its backend implementation instead
      return isinf(static_cast<Backend_>(x));
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

    /**
     * \brief Checks whether a value is Not-A-Number.
     *
     * \note There exists no generic implementation for this function.
     *
     * \param[in] x The value to be checked for NaN.
     *
     * \returns \c true if \p x is NaN, otherwise \c false.
     */
    template<typename T_>
    inline bool isnan(T_ x);

    WRAP_STD_MATH1BRET(isnan)

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    inline bool isnan(__float128 x)
    {
      // https://chromium.googlesource.com/native_client/nacl-gcc/+/ng/master/libquadmath/math/isnanq.c
      return (::isnanq(x) != 0);
    }
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_HALFMATH) && !defined(__CUDACC__)
    inline bool isnan(half_float::half x)
    {
      return half_float::isnan(x);
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_FLOATX) && !defined(__CUDACC__)
    template<int exp_bits_, int sig_bits_, typename Backend_>
    inline bool isnan(const flx::floatx<exp_bits_, sig_bits_, Backend_>& x)
    {
      // FloatX doesn't offer its own isnan implementation,
      // so test its backend implementation instead
      return isnan(static_cast<Backend_>(x));
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

    /**
     * \brief Checks whether a value is normal.
     *
     * A value is \e normal if it is neither infinity, NaN, zero or subnormal.
     *
     * \note There exists no generic implementation for this function.
     *
     * \param[in] x The value to be checked for normalty.
     *
     * \returns \c true if \p x is normal, otherwise \c false.
     */
    template<typename T_>
    inline bool isnormal(T_ x);

    WRAP_STD_MATH1BRET(isnormal)

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    inline bool isnormal(__float128 x)
    {
      // check whether the value is finite
      if(::finiteq(x) == 0)
        return false;
      // check whether the value is not below minimal normal value
      return !(::fabsq(x) < FLT128_MIN);
    }
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_HALFMATH) && !defined(__CUDACC__)
    inline bool isnormal(half_float::half x)
    {
      return half_float::isnormal(x);
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

#if defined(FEAT_HAVE_FLOATX) && !defined(__CUDACC__)
    template<int exp_bits_, int sig_bits_, typename Backend_>
    inline bool isnormal(const flx::floatx<exp_bits_, sig_bits_, Backend_>& x)
    {
      // FloatX doesn't offer its own isnormal implementation,
      // so test its backend implementation instead
      return isnormal(static_cast<Backend_>(x));
    }
#endif // FEAT_HAVE_HALFMATH && !__CUDA_CC__

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
      static_assert(Type::Traits<T_>::is_int, "factorial can only be applied to integral types");

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
      static_assert(Type::Traits<T_>::is_int, "binomial can only be applied to integral types");

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
     * \brief Inverts a matrix and returns its determinant.
     *
     * This function inverts a dense n x n matrix by means of partially pivoted Gaussian elimination
     * and returns its determinant.
     *
     * \attention
     * This function does not check whether the input matrix is regular, therefore one should
     * always check whether the determinant returned by this function is \e normal by using
     * the Math::isnormal() function. If the returned determinant is not normal, i.e. if
     * Math::isnormal returns \c false, then the matrix inversion most failed and the output
     * matrix \p a most probably contains garbage!
     *
     * \note
     * This function returns zero if any of the input parameters are invalid.
     *
     * \param[in] n
     * The dimension of the matrix to be inverted. Must be > 0.
     *
     * \param[in] stride
     * The stride of the matrix. Must be >= n.
     *
     * \param[in] a
     * On entry, the matrix to be inverted. On exit, the inverse matrix. Must not be \c nullptr.
     *
     * \param[in] p
     * A temporary pivot array of length at least <b>n</b>. Must not be \c nullptr.
     *
     * \returns
     * The determinant of the input matrix \p a.
     *
     * <b>Implementational Details:</b>\n
     * This functions implements the 'classical' column-pivoted Gauss-Jordan
     * elimination algorithm, which computes the inverse of a matrix by
     * transforming the system
     *
     *            \f[ [~A~|~I~] \rightarrow [~I~|~A^{-1}~] \f]
     *
     * There are three major modifications in contrast to the original
     * widely known algorithm:
     *
     * First, this algorithm works in-situ, i.e. it successively replaces
     * columns of the left-hand-side input matrix by the corresponding
     * columns of the right-hand-side identity matrix during the column
     * elimination process.
     *
     * Second, this algorithm eliminates whole columns of the left-hand-side
     * matrix, whereas the original algorithm only eliminates all entries
     * below the pivot element. Due to this fact, no backward substitution
     * is necessary after the primary elimination process.
     *
     * Third, this algorithm does not physically swap rows during the
     * pivoting step, but simply \e remembers the pivoted positions in
     * a temporary pivot array \e p, which is initialized to identity.
     * So in the k-th iteration of the primary elimination loop:
     *
     * - p[k] contains the index of the row which is currently being
     *   used for elimination of column k
     * - p[j] for j < k contains the indices of all columns which have
     *   already been eliminated, i.e. the corresponding columns already
     *   contain the entries of the right-hand-side matrix
     * - p[j] for j > k contains the indices of all columns which have
     *   not yet been eliminated, i.e. the corresponding columns still
     *   hold the entries of the left-hand-side matrix and are therefore
     *   candidates for future pivoting steps.
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_>
    DT_ invert_matrix(const IT_ n, const IT_ stride, DT_ a[], IT_ p[])
    {
      // make sure that the parameters are valid
      if((n <= IT_(0)) || (stride < n) || (a == nullptr) || (p == nullptr))
        return DT_(0);

      // invert 1x1 explicitly
      if(n == IT_(1))
      {
        DT_ det = a[0];
        a[0] = DT_(1) / det;
        return det;
      }

      // initialize identity permutation
      for(IT_ i(0); i < n; ++i)
      {
        p[i] = i;
      }

      // initialize determinant to 1
      DT_ det = DT_(1);

      // primary column elimination loop
      for(IT_ k(0); k < n; ++k)
      {
        // step 1: find a pivot for the elimination of column k
        {
          // for this, we only check the rows p[j] with j >= k, as all
          // rows p[j] with j < k have already been eliminated and are
          // therefore not candidates for pivoting
          DT_ pivot = Math::abs(a[p[k]*stride + p[k]]);
          IT_ i = k;

          // loop over all unprocessed rows
          for(IT_ j(k+1); j < n; ++j)
          {
            // get our matrix value and check whether it can be a pivot
            DT_ abs_val = Math::abs(a[p[j]*stride + p[j]]);
            if(abs_val > pivot)
            {
              pivot = abs_val;
              i = j;
            }
          }

          // do we have to swap rows i and k?
          if(i > k)
          {
            // swap rows "virtually" by exchanging their permutation positions
            IT_ t = p[k];
            p[k] = p[i];
            p[i] = t;
          }
        }

        // compute pivot row offset
        const IT_ pk_off = p[k]*stride;

        // step 2: process pivot row
        {
          // update determinant by multiplying with the pivot element
          det *= a[pk_off + p[k]];

          // get our inverted pivot element
          const DT_ pivot = DT_(1) / a[pk_off + p[k]];

          // replace column entry by unit column entry
          a[pk_off + p[k]] = DT_(1);

          // divide the whole row by the inverted pivot
          for(IT_ j(0); j < n; ++j)
          {
            a[pk_off+j] *= pivot;
          }
        }

        // step 3: eliminate pivot column

        // loop over all rows of the matrix
        for(IT_ i(0); i < n; ++i)
        {
          // skip the pivot row
          if(i == p[k])
            continue;

          // compute row and pivot offsets
          const IT_ row_off = i*stride;

          // fetch elimination factor
          const DT_ factor =  a[row_off + p[k]];

          // replace by unit column entry
          a[row_off + p[k]] = DT_(0);

          // process the row
          for(IT_ j(0); j < n; ++j)
          {
            a[row_off + j] -= a[pk_off + j] * factor;
          }
        }
      }

      // return determinant
      return det;
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
    }; // class Limits<...>

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
    template<>
    class Limits<__float128>
    {
    public:
      static constexpr bool is_specialized = true;
      static /*constexpr*/ __float128 min() noexcept { return FLT128_MIN; }
      static /*constexpr*/ __float128 max() noexcept { return FLT128_MAX; }
      static /*constexpr*/ __float128 lowest() noexcept { return -FLT128_MAX; }
      static constexpr int digits = FLT128_MANT_DIG;
      static constexpr int digits10 = FLT128_DIG;
      // Note: The following formula was taken from the MSC implementation...
      static constexpr int max_digits10 = (2 + FLT128_MANT_DIG * 301 / 1000);
      static constexpr bool is_signed = true;
      static constexpr bool is_integer = false;
      static constexpr bool is_exact = false;
      static constexpr int radix = 2;
      static /*constexpr*/ __float128 epsilon() noexcept { return FLT128_EPSILON; }
      static constexpr __float128 round_error() noexcept { return __float128(0.5); }
      static constexpr int min_exponent = FLT128_MIN_EXP;
      static constexpr int min_exponent10 = FLT128_MIN_10_EXP;
      static constexpr int max_exponent = FLT128_MAX_EXP;
      static constexpr int max_exponent10 = FLT128_MAX_10_EXP;
      static constexpr bool has_infinity = true;
      static constexpr bool has_quiet_NaN = true;
      static constexpr bool has_signaling_NaN = false;
      static constexpr std::float_denorm_style has_denorm = std::denorm_absent;
      static constexpr bool has_denorm_loss = false;
      static /*constexpr*/ __float128 infinity() noexcept { return max()*max(); }
      static /*constexpr*/ __float128 quiet_NaN() noexcept { return ::nanq(nullptr); }
      static /*constexpr*/ __float128 signaling_NaN() noexcept { return ::nanq(nullptr); }
      static /*constexpr*/ __float128 denorm_min() noexcept { return FLT128_DENORM_MIN; }
      static constexpr bool is_iec559 = false;
      static constexpr bool is_bounded = true;
      static constexpr bool is_modulo = false;
      static constexpr bool traps = true;
      static constexpr bool tinyness_before = true;
      static constexpr std::float_round_style round_style = std::round_to_nearest;
    };
#endif // FEAT_HAVE_QUADMATH && !__CUDA_CC__
  } // namespace Math
} // namespace FEAT

#endif // KERNEL_UTIL_MATH_HPP

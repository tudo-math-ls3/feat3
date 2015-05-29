#pragma once
#ifndef KERNEL_LINEAR_ALGEBRA_HPP
#define KERNEL_LINEAR_ALGEBRA_HPP 1

/**
 * \file
 * \brief Basic Linear Algebra header.
 *
 * This header file contains a large set of the most common linear algebra functions.
 *
 * \author Peter Zajac
 */

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>

// includes, system
#include <complex>

namespace FEAST
{
  /**
   * \brief Linear Algebra namespace
   *
   * This namespace encapsulates a set of the most common basic linear algebra functions, which is comparable
   * to the famous BLAS library.
   */
  namespace LinAlg
  {
    /// \cond internal
    namespace Intern
    {
      /// helper function: calculates the square of a value
      template<typename T_>
      inline T_ sqr(const T_& x)
      {
        return x*x;
      }

      /// helper function: swaps two values
      template<typename TypeX_>
      inline void swap(
        TypeX_& x,
        TypeX_& y)
      {
        TypeX_ t = x;
        x = y;
        y = t;
      }
    } // namespace Intern
    /// \endcond

    /**
     * \brief Clears a vector, i.e. sets all entries of a vector to a given value.
     *
     * This function performs the following operation:
     * \f[ \forall\ 0\leq i< n :\ x(i) := \alpha \f]
     *
     * \param[in] n
     * The length of the vector \p x.
     *
     * \param[out] x
     * The vector that is to be cleared.
     *
     * \param[in] alpha
     * The value to which the entries of \p x are to be set to.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeSize_>
    inline void vec_clear(
      const TypeSize_ n,
      TypeX_ x[],
      const TypeX_ alpha = TypeX_(0))
    {
      CONTEXT("LinAlg::vec_clear()");
      ASSERT_(x != nullptr);

      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        x[i] = alpha;
      }
    }

    /**
     * \brief Copies one vector into another one.
     *
     * This function performs the following operation:
     * \f[ y := x \f]
     *
     * \param[in] n
     * The length of the vectors \p x and \p y.
     *
     * \param[out] y
     * The destination vector which recieves the copy of \p x.
     *
     * \param[in] x
     * The source vector that is to be copied to \p y.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeY_,
      typename TypeSize_>
    inline void vec_copy(
      const TypeSize_ n,
      TypeY_ y[],
      const TypeX_ x[])
    {
      CONTEXT("LinAlg::vec_copy()");
      ASSERT_(x != nullptr);
      ASSERT_(y != nullptr);

      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        y[i] = TypeY_(x[i]);
      }
    }

    /**
     * \brief Swaps to vectors.
     *
     * This function performs a deep swap of two dense vectors.
     *
     * \param[in] n
     * The length of the vectors \p x and \p y.
     *
     * \param[in,out] x,y
     * The vectors which are to be swapped.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeY_,
      typename TypeSize_>
    inline void vec_swap(
      const TypeSize_ n,
      TypeX_ x[],
      TypeY_ y[])
    {
      CONTEXT("LinAlg::vec_swap()");
      ASSERT_(x != nullptr);
      ASSERT_(y != nullptr);

      TypeX_ t;
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        t = x[i];
        x[i] = TypeX_(y[i]);
        y[i] = TypeY_(t);
      }
    }

    /**
     * \brief Performs a left-shift of a vector.
     *
     * This function performs an in-situ left-shift, i.e.:
     * \f[ x(1,\dots,n) := x(k+1,\dots,k+n) \f]
     *
     * \remarks
     * The sub-vector <c>x(n+1,...,n+k)</c> is left unmodified.
     *
     * \param[in] n
     * The number of entries that are to be shifted.
     *
     * \param[in] k
     * The offset that \p x is to be shifted by.
     *
     * \param[in,out] x
     * The vector that is to be shifted. Is assumed to have at least <c>(k+n)</c> entries.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeSize_>
    inline void vec_shift_l(
      const TypeSize_ n,
      const TypeSize_ k,
      TypeX_ x[])
    {
      CONTEXT("LinAlg::vec_shift_l()");
      ASSERT_(x != nullptr);

      if(k > 0)
      {
        for(TypeSize_ i(0) ; i < n ; ++i)
        {
          x[i] = x[k + i];
        }
      }
    }

    /**
     * \brief Performs a right-shift of a vector.
     *
     * This function performs an in-situ right-shift, i.e.:
     * \f[ x(k+1,\dots,k+n) := x(1,\dots,n) \f]
     *
     * \remarks
     * The sub-vector <c>x(1,...,k)</c> is left unmodified.
     *
     * \param[in] n
     * The number of entries that are to be shifted.
     *
     * \param[in] k
     * The offset that \p x is to be shifted by.
     *
     * \param[in,out] x
     * The vector that is to be shifted. Is assumed to have at least <c>(k+n)</c> entries.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeSize_>
    inline void vec_shift_r(
      const TypeSize_ n,
      const TypeSize_ k,
      TypeX_ x[])
    {
      CONTEXT("LinAlg::vec_shift_r()");
      ASSERT_(x != nullptr);

      if(k > 0)
      {
        for(TypeSize_ i(n) ; i > 0 ; --i)
        {
          x[k + i - 1] = x[i - 1];
        }
      }
    }

    /**
     * \brief Scales a vector.
     *
     * \f[ x := \alpha \cdot x \f]
     *
     * \param[in] n
     * The length of the vector \p x.
     *
     * \param[in,out] x
     * The vector that is to be scaled.
     *
     * \param[in] alpha
     * The value by which the vector is to be scaled.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeSize_>
    inline void vec_scale(
      const TypeSize_ n,
      TypeX_ x[],
      const TypeX_ alpha)
    {
      CONTEXT("LinAlg::vec_scale()");
      ASSERT_(x != nullptr);

      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        x[i] *= alpha;
      }
    }

    /**
     * \brief Calculates the dot-product of two vectors.
     *
     * \tparam TypeR_
     * The desired return type of the dot-product.
     * \note
     * Independent of the types of \p x and \p y, the caller must always explicitly specify the return type \c TypeR_
     * of this function as a template parameter!
     *
     * \param[in] n
     * The length of the vectors \p x and \p y.
     *
     * \param[in] x,y
     * The vectors whose dot-product is to be calculated.
     *
     * \returns
     * The dot-product of \p x and \p y. If \p n == 0, then this function returns zero.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeR_,
      typename TypeX_,
      typename TypeY_,
      typename TypeSize_>
    inline TypeR_ vec_dot(
      const TypeSize_ n,
      const TypeX_ x[],
      const TypeY_ y[])
    {
      CONTEXT("LinAlg::vec_dot()");
      ASSERT_(x != nullptr);
      ASSERT_(y != nullptr);

      TypeR_ r = TypeR_(0);
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        r += TypeR_(x[i]) * TypeR_(y[i]);
      }

      return r;
    }

    /**
     * \brief Performs an AXPY operation of two vectors.
     *
     * This function performs the so-called AXPY operation:
     * \f[ y := y + \alpha \cdot x \f]
     *
     * \param[in] n
     * The length of the vectors \p x and \p y.
     *
     * \param[in,out] y
     * The vector that recieves the result of the operation.
     *
     * \param[in] x
     * The input vector that is to be added onto \p y.
     *
     * \param[in] alpha
     * The value by which \p x is to be scaled before being added onto \p y.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeY_,
      typename TypeSize_>
    inline void vec_axpy(
      const TypeSize_ n,
      TypeY_ y[],
      const TypeX_ x[],
      const TypeY_ alpha = TypeY_(1))
    {
      CONTEXT("LinAlg::vec_axpy()");
      ASSERT_(y != nullptr);
      ASSERT_(x != nullptr);

      /// \compilerhack Intel C++ 14 loop vectorisation bug
#if defined(FEAST_COMPILER_INTEL) && (FEAST_COMPILER_INTEL < 1500)
      for(TypeSize_ i(0) ; (i+1) < (n+1) ; ++i)
#else
      for(TypeSize_ i(0) ; i < n ; ++i)
#endif
      {
        y[i] += alpha * TypeY_(x[i]);
      }
    }

    /**
     * \brief Calculates a linear combination of two vectors.
     *
     * This function performs the following operation:
     * \f[ z := \alpha\cdot x + \beta\cdot y \f]
     *
     * \param[in] n
     * The length of the vectors.
     *
     * \param[out] z
     * The vector that recieves the result of the operation.
     *
     * \param[in] x, y
     * The vectors whose linear combination is to be calculated.
     *
     * \param[in] alpha, beta
     * The scaling factors for \p x and \p y.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeY_,
      typename TypeZ_,
      typename TypeSize_>
    inline void vec_lin_comb(
      const TypeSize_ n,
      TypeZ_ z[],
      const TypeX_ x[],
      const TypeY_ y[],
      const TypeZ_ alpha,
      const TypeZ_ beta)
    {
      CONTEXT("LinAlg::vec_lin_comb()");
      ASSERT_(x != nullptr);
      ASSERT_(y != nullptr);
      ASSERT_(z != nullptr);

      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        z[i] = alpha * TypeZ_(x[i]) + beta * TypeZ_(y[i]);
      }
    }

    /**
     * \brief Calculates a component-wise product of two vectors.
     *
     * This function performs the following operation:
     * \f[ \forall\ 0\leq i< n :\ z(i) := \alpha\cdot x(i)\cdot y(i) \f]
     *
     * \param[in] n
     * The length of the vectors.
     *
     * \param[out] z
     * The vector that recieves the result of the operation.
     *
     * \param[in] x, y
     * The vectors which are to be multiplied component-wise.
     *
     * \param[in] alpha
     * The scaling factor for \p x and \p y.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeY_,
      typename TypeZ_,
      typename TypeSize_>
    inline void vec_comp_mult(
      const TypeSize_ n,
      TypeZ_ z[],
      const TypeX_ x[],
      const TypeY_ y[],
      const TypeZ_ alpha)
    {
      CONTEXT("LinAlg::vec_comp_mult()");
      ASSERT_(x != nullptr);
      ASSERT_(y != nullptr);
      ASSERT_(z != nullptr);

      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        z[i] = alpha * TypeZ_(x[i]) * TypeZ_(y[i]);
      }
    }

    /**
     * \brief Calculates the absolute-sum norm of a vector.
     *
     * This function returns:
     * \f[ \|x\|_1 := \sum_{0\leq i<n} |x(i)| \f]
     *
     * \param[in] n
     * The length of the vector \p x.
     *
     * \param[in] x
     * The vector whose norm is to be calculated.
     *
     * \returns
     * The absolute-sum norm of \p x.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeSize_>
    inline TypeX_ vec_norm_asum(
      const TypeSize_ n,
      const TypeX_ x[])
    {
      CONTEXT("LinAlg::vec_norm_asum()");
      ASSERT_(x != nullptr);

      TypeX_ r = TypeX_(0);
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        r += Math::abs(x[i]);
      }

      return r;
    }

    /**
     * \brief Calculates the absolute-sum norm of a complex vector.
     * \copydetails vec_norm_asum(const TypeSize_ n, const TypeX_ x[])
     */
    template<
      typename TypeReal_,
      typename TypeSize_>
    inline TypeReal_ vec_norm_asum(
      const TypeSize_ n,
      const std::complex<TypeReal_> x[])
    {
      CONTEXT("LinAlg::vec_norm_asum() [complex version]");
      ASSERT_(x != nullptr);

      TypeReal_ r = TypeReal_(0);
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        r += Math::abs(x[i]);
      }

      return r;
    }

    /**
     * \brief Calculates the euclid norm of a vector.
     *
     * This function returns:
     * \f[ \|x\|_2 := \Big(\sum_{0\leq i<n} |x(i)|^2\Big)^{\frac{1}{2}} \f]
     *
     * \param[in] n
     * The length of the vector \p x.
     *
     * \param[in] x
     * The vector whose norm is to be calculated.
     *
     * \returns
     * The euclid norm of \p x.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeSize_>
    inline TypeX_ vec_norm_euclid(
      const TypeSize_ n,
      const TypeX_ x[])
    {
      using Intern::sqr;
      CONTEXT("LinAlg::vec_norm_euclid()");
      ASSERT_(x != nullptr);

      TypeX_ r = TypeX_(0);
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        r += sqr(x[i]);
      }

      return Math::sqrt(r);
    }

    /**
     * \brief Calculates the euclid norm of a complex vector.
     * \copydetails vec_norm_euclid(const TypeSize_ n, const TypeX_ x[])
     */
    template<
      typename TypeReal_,
      typename TypeSize_>
    inline TypeReal_ vec_norm_euclid(
      const TypeSize_ n,
      const std::complex<TypeReal_> x[])
    {
      using Intern::sqr;
      CONTEXT("LinAlg::vec_norm_euclid() [complex version]");
      ASSERT_(x != nullptr);

      TypeReal_ r = TypeReal_(0);
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        // We don't use std::abs here, as this would involve calculating the square root of the right hand
        // expression below, which we would then need to square afterwards...
        r += sqr(x[i].real()) + sqr(x[i].imag());
      }

      return Math::sqrt(r);
    }

    /**
     * \brief Calculates the maximum norm of a vector.
     *
     * This function returns:
     * \f[ \|x\|_\infty := \max_{0\leq i< n} |x(i)| \f]
     *
     * \param[in] n
     * The length of the vector \p x.
     *
     * \param[in] x
     * The vector whose norm is to be calculated.
     *
     * \returns
     * The maximum norm of \p x.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeX_,
      typename TypeSize_>
    inline TypeX_ vec_norm_max(
      const TypeSize_ n,
      const TypeX_ x[])
    {
      CONTEXT("LinAlg::vec_norm_max()");
      ASSERT_(x != nullptr);

      TypeX_ t, r = TypeX_(0);
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        t = Math::abs(x[i]);
        r = r < t ? t : r;
      }

      return r;
    }

    /**
     * \brief Calculates the maximum norm of a complex vector.
     * \copydetails vec_norm_max(const TypeSize_ n, const TypeX_ x[])
     */
    template<
      typename TypeReal_,
      typename TypeSize_>
    inline TypeReal_ vec_norm_max(
      const TypeSize_ n,
      const std::complex<TypeReal_> x[])
    {
      using Intern::sqr;
      CONTEXT("LinAlg::vec_norm_max() [complex version]");
      ASSERT_(x != nullptr);

      // The square root is monotone, so we'll apply it at the end instead of in every loop iteration.
      TypeReal_ t, r = TypeReal_(0);
      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        t = sqr(x[i].real()) + sqr(x[i].imag());
        r = r < t ? t : r;
      }

      return Math::sqrt(r);
    }

    /**
     * \brief Clears a dense matrix.
     *
     * This function clears a dense matrix, i.e. sets all entries of a matrix to a specified value.
     *
     * \param[in] m
     * The number of rows of \p a.
     *
     * \param[in] n
     * The number of columns of \p a.
     *
     * \param[in] stride_a
     * The stride of \p a.
     *
     * \param[out] a
     * The matrix which is to be cleared.
     *
     * \param[in] alpha
     * The value to which the entries of \p a are to be set to.
     *
     * \see vec_clear
     *
     * \author Peter Zajac
     */
    template<
      typename TypeA_,
      typename TypeSize_>
    inline void mat_clear(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride_a,
      TypeA_ a[],
      const TypeA_ alpha = TypeA_(0))
    {
      CONTEXT("LinAlg::mat_clear()");
      ASSERT_(stride_a >= n);
      ASSERT_(a != nullptr);

      if(stride_a == n)
      {
        vec_clear(m * n, a, alpha);
      }
      else
      {
        for(TypeSize_ i(0) ; i < m ; ++i)
        {
          vec_clear(n, &a[i * stride_a], alpha);
        }
      }
    }

    /**
     * \brief Creates an identity matrix.
     *
     * \param[in] n
     * The number of rows/columns of \p a.
     *
     * \param[in] stride_a
     * The stride of \p a.
     *
     * \param[out] a
     * The matrix that is to be set to the identity matrix.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeA_,
      typename TypeSize_>
    inline void mat_identity(
      const TypeSize_ n,
      const TypeSize_ stride_a,
      TypeA_ a[])
    {
      CONTEXT("LinAlg::mat_identity()");
      ASSERT_(stride_a >= n);
      ASSERT_(a != nullptr);

      for(TypeSize_ i(0) ; i < n ; ++i)
      {
        TypeSize_ k = i * stride_a;
        for(TypeSize_ j(0) ; j < i ; ++j)
        {
          a[k + j] = TypeA_(0);
        }
        a[k + i] = TypeA_(1);
        for(TypeSize_ j(i + 1) ; j < n ; ++j)
        {
          a[k + j] = TypeA_(0);
        }
      }
    }

    /// \cond internal
    namespace Intern
    {
      // matrix copy helper class
      template<bool>
      struct MatCopy;

      // matrix copy helper: B := A
      template<>
      struct MatCopy<false>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeSize_>
        static inline void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_b,
          TypeB_ b[],
          const TypeSize_ stride_a,
          const TypeA_ a[])
        {
          CONTEXT("LinAlg::Intern::MatCopy<false>::apply()");
          if((stride_a == n) && (stride_b == n))
          {
            vec_copy(m * n, b, a);
          }
          else
          {
            for(TypeSize_ i(0) ; i < m ; ++i)
            {
              vec_copy(n, &b[i * stride_b], &a[i * stride_a]);
            }
          }
        }
      }; // class MatCopy<false>

      // matrix copy helper: B := A^T
      template<>
      struct MatCopy<true>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeSize_>
        static inline void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_b,
          TypeB_ b[],
          const TypeSize_ stride_a,
          const TypeA_ a[])
        {
          CONTEXT("LinAlg::Intern::MatCopy<true>::apply()");
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            for(TypeSize_ j(0) ; j < n ; ++j)
            {
              b[i * stride_b + j] = TypeB_(a[j * stride_a + i]);
            }
          }
        }
      }; // class MatCopy<true>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Copies a dense matrix into another one.
     *
     * This function peforms the following operation:
     * \f[ B := OP(A) \f]
     *
     * \tparam trans_a_
     * Specifies whether the input matrix is to be transposed or not.
     *
     * \li <c>trans_a_ = false</c>: <em>OP(A) = A</em>
     * \li <c>trans_a_ = true </c>: <em>OP(A) = A^T</em>
     *
     * \param[in] m
     * The number of rows of \e B and the number of rows of <em>OP(A)</em>, respectively.
     *
     * \param[in] n
     * The number of columns of \e B and the number of columns of <em>OP(A)</em>, respectively.
     *
     * \param[in] stride_b
     * The stride of \e B.
     *
     * \param[out] b
     * The m-by-n matrix \e B that recieves the copy of <em>OP(A)</em>.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in] a
     * The matrix \e A that is to be copied.
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      typename TypeA_,
      typename TypeB_,
      typename TypeSize_>
    inline void mat_copy(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride_b,
      TypeB_ b[],
      const TypeSize_ stride_a,
      const TypeA_ a[])
    {
      CONTEXT("LinAlg::mat_copy()");
      ASSERT_(stride_b >= n);
      ASSERT_(trans_a_ ? (stride_a >= m) : (stride_a >= n));
      ASSERT_(b != nullptr);
      ASSERT_(a != nullptr);

      Intern::MatCopy<trans_a_>::apply(m, n, stride_b, b, stride_a, a);
    }

    /// \cond internal
    namespace Intern
    {
      // matrix swap helper class
      template<bool>
      struct MatSwap;

      // matrix swap helper: B :=: A
      template<>
      struct MatSwap<false>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeSize_>
        static inline void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_a,
          TypeA_ a[],
          const TypeSize_ stride_b,
          TypeB_ b[])
        {
          CONTEXT("LinAlg::Intern::MatSwap<false>::apply()");
          if((stride_a == n) && (stride_b == n))
          {
            vec_swap(m * n, a, b);
          }
          else
          {
            for(TypeSize_ i(0) ; i < m ; ++i)
            {
              vec_swap(n, &a[i * stride_a], &b[i * stride_b]);
            }
          }
        }
      }; // class MatSwap<false>

      // matrix swap helper: B :=: A^T
      template<>
      struct MatSwap<true>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeSize_>
        static inline void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_a,
          TypeA_ a[],
          const TypeSize_ stride_b,
          TypeB_ b[])
        {
          CONTEXT("LinAlg::Intern::MatSwap<true>::apply()");
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            for(TypeSize_ j(0) ; j < n ; ++j)
            {
              TypeB_ t = b[i * stride_b + j];
              b[i * stride_b + j] = TypeB_(a[j * stride_a + i]);
              a[j * stride_a + i] = TypeA_(t);
            }
          }
        }
      }; // class MatSwap<false>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Swaps two dense matrices.
     *
     * This function performs a deep swap of two dense matrices.
     * \f[ B :=: OP(A) \f]
     *
     * \tparam trans_a_
     * Specifies whether the matrix A is to be transposed or not.
     *
     * \li <c>trans_a_ = false</c>: <em>OP(A) = A</em>
     * \li <c>trans_a_ = true </c>: <em>OP(A) = A^T</em>
     *
     * \param[in] m
     * The number of rows of \e B and the number of rows of <em>OP(A)</em>, respectively.
     *
     * \param[in] n
     * The number of columns of \e B and the number of columns of <em>OP(A)</em>, respectively.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in,out] a
     * The matrix \e A that is to be swapped with \e B.
     *
     * \param[in] stride_b
     * The stride of \e B.
     *
     * \param[in,out] b
     * The m-by-n matrix \e B that is to be swapped with \e A.
     *
     * \see vec_swap
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      typename TypeA_,
      typename TypeB_,
      typename TypeSize_>
    inline void mat_swap(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride_a,
      TypeA_ a[],
      const TypeSize_ stride_b,
      TypeB_ b[])
    {
      CONTEXT("LinAlg::mat_swap()");
      ASSERT_(trans_a_ ? stride_a >= m : stride_a >= n);
      ASSERT_(stride_b >= n);
      ASSERT_(a != nullptr);
      ASSERT_(b != nullptr);

      Intern::MatSwap<trans_a_>::apply(m, n, stride_a, a, stride_b, b);
    }

    /**
     * \brief Scales a dense matrix.
     *
     * This function performs the following operation:
     * \f[ A := \alpha\cdot A\f]
     *
     * \param[in] m
     * The number of rows of \e A.
     *
     * \param[in] n
     * The number of columns of \e A.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in,out] a
     * The m-by-n matrix \e A that is to be scaled.
     *
     * \param[in] alpha
     * The scaling factor.
     *
     * \see vec_scale
     *
     * \author Peter Zajac
     */
    template<
      typename TypeA_,
      typename TypeSize_>
    inline void mat_scale(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride_a,
      TypeA_ a[],
      const TypeA_ alpha)
    {
      CONTEXT("LinAlg::mat_scale()");
      ASSERT_(a != nullptr);
      ASSERT_(m >= 0);
      ASSERT_(n >= 0);
      ASSERT_(stride_a >= n);

      if(stride_a == n)
      {
        vec_scale(m * n, a, alpha);
      }
      else
      {
        for(TypeSize_ i(0) ; i < m ; ++i)
        {
          vec_scale(n, &a[i * stride_a], alpha);
        }
      }
    }

    /// \cond internal
    namespace Intern
    {
      // matrix axpy helper class
      template<bool>
      struct MatAxpy;

      // matrix axpy helper: B += alpha*A
      template<>
      struct MatAxpy<false>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeSize_>
        static inline void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_b,
          TypeB_ b[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeB_ alpha)
        {
          CONTEXT("LinAlg::Intern::MatAxpy<false>::apply()");
          if((stride_a == n) && (stride_b == n))
          {
            vec_axpy(m * n, b, a, alpha);
          }
          else
          {
            for(TypeSize_ i(0) ; i < m ; ++i)
            {
              vec_axpy(n, &b[i * stride_b], &a[i * stride_a], alpha);
            }
          }
        }
      }; // class MatAxpy<false>

      // matrix axpy helper: B += alpha*A^T
      template<>
      struct MatAxpy<true>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeSize_>
        static inline void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_b,
          TypeB_ b[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeB_ alpha)
        {
          CONTEXT("LinAlg::Intern::MatAxpy<true>::apply()");
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            for(TypeSize_ j(0) ; j < n ; ++j)
            {
              b[i * stride_b + j] += alpha*TypeB_(a[j * stride_a + i]);
            }
          }
        }
      }; // class MatAxpy<true>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Performs an AXPY operation of two dense matrices.
     *
     * This function performs the following operation:
     *   \f[ B := B + \alpha\cdot OP(A) \f]
     *
     * \tparam trans_a_
     * Specifies whether the input matrix is to be transposed or not.
     *
     * \li <c>trans_a_ = false</c>: <em>OP(A) = A</em>
     * \li <c>trans_a_ = true </c>: <em>OP(A) = A^T</em>
     *
     * \param[in] m
     * The number of rows of \e B and the number of rows of <em>OP(A)</em>, respectively.
     *
     * \param[in] n
     * The number of columns of \e B and the number of columns of <em>OP(A)</em>, respectively.
     *
     * \param[in] stride_b
     * The stride of \e B.
     *
     * \param[in,out] b
     * The m-by-n matrix \e B that recieves the result.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in] a
     * The matrix \e A that is to be added.
     *
     * \param[in] alpha
     * The scaling factor.
     *
     * \see vec_axpy
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      typename TypeA_,
      typename TypeB_,
      typename TypeSize_>
    inline void mat_axpy(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride_b,
      TypeB_ b[],
      const TypeSize_ stride_a,
      const TypeA_ a[],
      const TypeB_ alpha = TypeB_(1))
    {
      CONTEXT("LinAlg::mat_axpy()");
      ASSERT_(stride_b >= n);
      ASSERT_(b != nullptr);
      ASSERT_(trans_a_ ? stride_a >= m : stride_a >= n);
      ASSERT_(a != nullptr);

      Intern::MatAxpy<trans_a_>::apply(m, n, stride_b, b, stride_a, a, alpha);
    }

    /// \cond internal
    namespace Intern
    {
      // mat-mat mult helper class
      template<bool, bool>
      struct MatMatMult;

      // mat-mat mult helper: C += alpha*A*B
      template<>
      struct MatMatMult<false, false>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeC_,
          typename TypeSize_>
        inline static void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ l,
          const TypeSize_ stride_c,
          TypeC_ c[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeSize_ stride_b,
          const TypeB_ b[],
          const TypeC_ alpha)
        {
          CONTEXT("LinAlg::Intern::MatMatMult<false,false>::apply()");
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            for(TypeSize_ j(0) ; j < n ; ++j)
            {
              TypeC_ t = TypeC_(0);
              for(TypeSize_ k(0) ; k < l ; ++k)
              {
                t += TypeC_(a[i * stride_a + k]) * TypeC_(b[k * stride_b + j]);
              }
              c[i * stride_c + j] += alpha*t;
            }
          }
        }
      };

      // mat-mat mult helper: C += alpha*(A^T)*B
      template<>
      struct MatMatMult<true, false>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeC_,
          typename TypeSize_>
        inline static void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ l,
          const TypeSize_ stride_c,
          TypeC_ c[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeSize_ stride_b,
          const TypeB_ b[],
          const TypeC_ alpha)
        {
          CONTEXT("LinAlg::Intern::MatMatMult<true,false>::apply()");
          // \todo check whether this can be implemented more efficiently
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            for(TypeSize_ j(0) ; j < n ; ++j)
            {
              TypeC_ t = TypeC_(0);
              for(TypeSize_ k(0) ; k < l ; ++k)
              {
                t += TypeC_(a[k * stride_a + i]) * TypeC_(b[k * stride_b + j]);
              }
              c[i * stride_c + j] += alpha*t;
            }
          }
        }
      };

      // mat-mat mult helper: C += alpha*A*(B^T)
      template<>
      struct MatMatMult<false, true>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeC_,
          typename TypeSize_>
        inline static void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ l,
          const TypeSize_ stride_c,
          TypeC_ c[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeSize_ stride_b,
          const TypeB_ b[],
          const TypeC_ alpha)
        {
          CONTEXT("LinAlg::Intern::MatMatMult<false,true>::apply()");
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            for(TypeSize_ j(0) ; j < n ; ++j)
            {
              c[i * stride_c + j] += alpha * vec_dot<TypeC_>(l, &a[i * stride_a], &b[j * stride_b]);
            }
          }
        }
      };

      // mat-mat mult helper: C += alpha*(A^T)*(B^T)
      template<>
      struct MatMatMult<true, true>
      {
        template<
          typename TypeA_,
          typename TypeB_,
          typename TypeC_,
          typename TypeSize_>
        inline static void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ l,
          const TypeSize_ stride_c,
          TypeC_ c[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeSize_ stride_b,
          const TypeB_ b[],
          const TypeC_ alpha)
        {
          CONTEXT("LinAlg::Intern::MatMatMult<true,true>::apply()");
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            for(TypeSize_ j(0) ; j < n ; ++j)
            {
              TypeC_ t = TypeC_(0);
              for(TypeSize_ k(0) ; k < l ; ++k)
              {
                t += TypeC_(a[k * stride_a + i]) * TypeC_(b[j * stride_b + k]);
              }
              c[i * stride_c + j] += alpha*t;
            }
          }
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Calculates a matrix-product of two dense matrices.
     *
     * This function performs the following operation:
     * \f[ C := C + \alpha \cdot OP(A) \cdot OP(B) \f]
     *
     * \tparam trans_a_, trans_b_
     * Specifies whether the matrices \p A and \p B are to be transposed or not.
     *
     * \li <c>trans_a_ = false</c>: <em>OP(A) = A</em>
     * \li <c>trans_a_ = true </c>: <em>OP(A) = A^T</em>
     * \li <c>trans_b_ = false</c>: <em>OP(B) = B</em>
     * \li <c>trans_b_ = true </c>: <em>OP(B) = B^T</em>
     *
     * \param[in] m
     * The number of rows of \e C and the number or rows of <em>OP(A)</em>, respectively.
     *
     * \param[in] n
     * The number of columns of \e C and the number or columns of <em>OP(B)</em>, respectively.
     *
     * \param[in] l
     * The number of columns of <em>OP(A)</em> and the number of rows of <em>OP(B)</em>, respectively.
     *
     * \param[in] stride_c
     * The stride of \e C.
     *
     * \param[in,out] c
     * The matrix \e C that recieves the result of the operation.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in] a
     * The left multiplicant matrix \e A.
     *
     * \param[in] stride_b
     * The stride of \e B.
     *
     * \param[in] b
     * The right multiplicant matrix \e B.
     *
     * \param[in] alpha
     * The scaling factor.
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      bool trans_b_,
      typename TypeA_,
      typename TypeB_,
      typename TypeC_,
      typename TypeSize_>
    inline void mat_mat_mult(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ l,
      const TypeSize_ stride_c,
      TypeC_ c[],
      const TypeSize_ stride_a,
      const TypeA_ a[],
      const TypeSize_ stride_b,
      const TypeB_ b[],
      const TypeC_ alpha = TypeC_(1))
    {
      CONTEXT("LinAlg::mat_mat_mult()");
      ASSERT_(stride_c >= n);
      ASSERT_(c != nullptr);
      ASSERT_(trans_a_ ? stride_a >= m : stride_a >= l);
      ASSERT_(a != nullptr);
      ASSERT_(trans_b_ ? stride_b >= l : stride_b >= n);
      ASSERT_(b != nullptr);

      Intern::MatMatMult<trans_a_, trans_b_>::apply(m, n, l, stride_c, c, stride_a, a, stride_b, b, alpha);
    }

    /// \cond internal
    namespace Intern
    {
      // mat-vec mult helper class
      template<bool>
      struct MatVecMult;

      // mat-vec mult helper: y += alpha*A*x
      template<>
      struct MatVecMult<false>
      {
        template<
          typename TypeA_,
          typename TypeX_,
          typename TypeY_,
          typename TypeSize_>
        inline static void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_a,
          TypeY_ y[],
          const TypeA_ a[],
          const TypeX_ x[],
          const TypeY_ alpha = TypeY_(1))
        {
          CONTEXT("LinAlg::Intern::MatVecMult<false>::apply()");
          for(TypeSize_ i(0) ; i < m ; ++i)
          {
            y[i] += alpha * vec_dot<TypeY_>(n, x, &a[i * stride_a]);
          }
        }
      };

      // mat-vec mult helper: y += alpha*(A^T)*x
      template<>
      struct MatVecMult<true>
      {
        template<
          typename TypeA_,
          typename TypeX_,
          typename TypeY_,
          typename TypeSize_>
        inline static void apply(
          const TypeSize_ m,
          const TypeSize_ n,
          const TypeSize_ stride_a,
          TypeY_ y[],
          const TypeA_ a[],
          const TypeX_ x[],
          const TypeY_ alpha = TypeY_(1))
        {
          CONTEXT("LinAlg::Intern::MatVecMult<true>::apply()");
          for(TypeSize_ j(0) ; j < n ; ++j)
          {
            vec_axpy(m, y, &a[j * stride_a], alpha*TypeY_(x[j]));
          }
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Performs a matrix-vector-multiplication with a dense matrix.
     *
     * This function performs the following operation:
     * \f[ y := y + \alpha \cdot OP(A) \cdot x \f]
     *
     * \tparam trans_a_
     * Specifies whether the matrix \e A is to be transposed or not.
     *
     * \li <c>trans_a_ = false</c>: <em>OP(A) = A</em>
     * \li <c>trans_a_ = true </c>: <em>OP(A) = A^T</em>
     *
     * \param[in] m
     * The length of \e y and the number of rows of <em>OP(A)</em>, respectively.
     *
     * \param[in] n
     * The length of \e x and the number of columns of <em>OP(A)</em>, respectively.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in,out] y
     * The vector \e y that recieves the result of the operation.
     *
     * \param[in] a
     * The matrix \e A that is to be used.
     *
     * \param[in] x
     * The vector \e x that \e A is to be multiplied by.
     *
     * \param[in] alpha
     * The scaling factor.
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      typename TypeA_,
      typename TypeX_,
      typename TypeY_,
      typename TypeSize_>
    inline void mat_vec_mult(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride_a,
      TypeY_ y[],
      const TypeA_ a[],
      const TypeX_ x[],
      const TypeY_ alpha = TypeY_(1))
    {
      CONTEXT("LinAlg::mat_vec_mult()");
      ASSERT_(a != nullptr);
      ASSERT_(x != nullptr);
      ASSERT_(y != nullptr);
      ASSERT_(trans_a_ ? stride_a >= m : stride_a >= n);

      Intern::MatVecMult<trans_a_>::apply(m, n, stride_a, y, a, x, alpha);
    }

    /**
     * \brief Calculates the row-sum norm of a matrix.
     *
     * This function returns:
     * \f[ \|A\|_\infty := \max_{0\leq i < m} \sum_{0\leq j<n} |a(i,j)| \f]
     *
     * \param[in] m
     * The number of rows of \e A.
     *
     * \param[in] n
     * The number of columns of \e A.
     *
     * \param[in] stride
     * The stride of \e A.
     *
     * \param[in] a
     * The matrix \e A whose norm is to be calculated.
     *
     * \returns
     * The row-sum-norm of \e A.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeA_,
      typename TypeSize_>
    inline TypeA_ mat_norm_row_sum(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride,
      const TypeA_ a[])
    {
      CONTEXT("LinAlg::mat_norm_row_sum()");
      TypeA_ r = TypeA_(0);
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        r = Math::max(r, vec_norm_asum(n, &a[i * stride]));
      }
      return r;
    }

    /**
     * \brief Calculates the row-sum norm of a complex matrix.
     * \copydetails mat_norm_row_sum(const TypeSize_ m, const TypeSize_ n, const TypeSize_ stride, const TypeA_ a[])
     */
    template<
      typename TypeReal_,
      typename TypeSize_>
    inline TypeReal_ mat_norm_row_sum(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride,
      const std::complex<TypeReal_> a[])
    {
      CONTEXT("LinAlg::mat_norm_row_sum() [complex version]");
      TypeReal_ r = TypeReal_(0);
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        r = Math::max(r, vec_norm_asum(n, &a[i * stride]));
      }
      return r;
    }

    /**
     * \brief Calculates the Frobenius norm of a matrix.
     *
     * This function returns:
     * \f[ \|A\|_F := \Big(\sum_{0\leq i<m}\sum_{0\leq j<n} |a(i,j)|^2\Big)^\frac{1}{2} \f]
     *
     * \param[in] m
     * The number of rows of \e A.
     *
     * \param[in] n
     * The number of columns of \e A.
     *
     * \param[in] stride
     * The stride of \e A.
     *
     * \param[in] a
     * The matrix \e A whose norm is to be calculated.
     *
     * \returns
     * The Forbenius norm of \e A.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeA_,
      typename TypeSize_>
    inline TypeA_ mat_norm_frobenius(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride,
      const TypeA_ a[])
    {
      using Intern::sqr;
      CONTEXT("LinAlg::mat_norm_frobenius()");
      TypeA_ r = TypeA_(0);
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        TypeSize_ k = i * stride;
        for(TypeSize_ j(0) ; j < n ; ++j)
        {
          r += sqr(a[k+j]);
        }
      }
      return Math::sqrt(r);
    }

    /**
     * \brief Calculates the Forbenius norm of a complex matrix.
     * \copydetails mat_norm_frobenius(const TypeSize_ m, const TypeSize_ n, const TypeSize_ stride, const TypeA_ a[])
     */
    template<
      typename TypeReal_,
      typename TypeSize_>
    inline TypeReal_ mat_norm_frobenius(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride,
      const std::complex<TypeReal_> a[])
    {
      using Intern::sqr;
      CONTEXT("LinAlg::()mat_norm_frobenius [complex version]");
      TypeReal_ r = TypeReal_(0);
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        TypeSize_ k = i * stride;
        for(TypeSize_ j(0) ; j < n ; ++j)
        {
          r += sqr(a[k+j].real()) + sqr(a[k+j].imag());
        }
      }
      return Math::sqrt(r);
    }

    /**
     * \brief Calculates the maximum-norm of a matrix.
     *
     * This function returns:
     * \f[ \|A\| := \max_{0\leq i < m} \max_{0\leq j< n} |a(i,j)| \f]
     *
     * \param[in] m
     * The number of rows of \e A.
     *
     * \param[in] n
     * The number of columns of \e A.
     *
     * \param[in] stride
     * The stride of \e A.
     *
     * \param[in] a
     * The matrix \e A whose norm is to be calculated.
     *
     * \returns
     * The maximum-norm of \e A.
     *
     * \author Peter Zajac
     */
    template<
      typename TypeA_,
      typename TypeSize_>
    inline TypeA_ mat_norm_max(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride,
      const TypeA_ a[])
    {
      CONTEXT("LinAlg::mat_norm_max()");
      TypeA_ r = TypeA_(0);
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        r = Math::max(r, vec_norm_max(n, &a[i * stride]));
      }
      return r;
    }

    /**
     * \brief Calculates the maximum-norm of a complex matrix.
     * \copydetails mat_norm_max(const TypeSize_ m, const TypeSize_ n, const TypeSize_ stride, const TypeA_ a[])
     */
    template<
      typename TypeReal_,
      typename TypeSize_>
    inline TypeReal_ mat_norm_max(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride,
      const std::complex<TypeReal_> a[])
    {
      CONTEXT("LinAlg::mat_norm_max() [complex version]");
      TypeReal_ r = TypeReal_(0);
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        r = Math::max(r, vec_norm_max(n, &a[i * stride]));
      }
      return r;
    }

    /**
     * \brief Calculates a pivoted LU-factorisation of a dense matrix.
     *
     * This function calculates a row-pivoted LU-factorisation of a given matrix A, i.e.
     * \f[ A = P\cdot L\cdot U\f]
     * where
     * \li P is a permutation matrix
     * \li L is a unit lower triangular (trapezoidal if \c m > \c n) matrix
     * \li U is an upper triangular (trapezoidal if \c m < \c n) matrix.
     *
     * The input matrix \c A is overwritten by its \c L and \c U factors, whereas the permutation matrix \c P is
     * encoded into a pivot array \c p.
     *
     * \remark
     * Although this function is also capable of factorising <em>non-square</em> matrices, i.e. the case where
     * <c>m != n</c>, there are (currently) no functions which make use of such a factorisation.
     *
     * \param[in] m
     * The number of rows of the matrix \e A.
     *
     * \param[in] n
     * The number of columns of the matrix \e A.
     *
     * \param[in] stride_a
     * The stride of the matrix \e A.
     *
     * \param[in,out] a
     * On entry, the \p m by \p n matrix \e A that is to be factorised.\n
     * On exit, the \e L and \e U factors of the PLU factorisation.
     *
     * \param[out] p
     * The pivot array which represents the permutation matrix \e P of the PLU factorisation.
     * Must be at least of length \c m.
     *
     * \returns
     * \c true if the factorisation was successful or \c false if a zero pivot was encountered.
     *
     * \see mat_solve_vec, mat_solve_mat
     *
     * \author Peter Zajac
     */
    template<
      typename TypeA_,
      typename TypeP_,
      typename TypeSize_>
    bool mat_factorise(
      const TypeSize_ m,
      const TypeSize_ n,
      const TypeSize_ stride_a,
      TypeA_ a[],
      TypeP_ p[])
    {
      CONTEXT("LinAlg::mat_factorise()");
      ASSERT_(a != nullptr);
      ASSERT_(p != nullptr);
      ASSERT_(stride_a >= n);

      // loop over all rows of A
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        // find a pivot
        p[i] = TypeP_(i);
        TypeA_ tp = Math::abs(a[i * (stride_a + 1)]);
        for(TypeSize_ j(i + 1) ; j < n ; ++j)
        {
          TypeA_ t = Math::abs(a[j * stride_a + i]);
          if(t > tp)
          {
            // row j is our new candidate for the pivot row
            tp = t;
            p[i] = TypeP_(j);
          }
        }

        // swap rows if necessary
        if(p[i] != TypeP_(i))
        {
          vec_swap(n, &a[i * stride_a], &a[p[i] * stride_a]);
        }

        // try to invert the pivot
        /// \todo Check whether there is a more reliable approach to check for zero-pivots.
        try
        {
          // note: we store the inverse of the main diagonal element rather than the element itself to perform
          //       multiplications rather than divisions later on.
          a[i * (stride_a + 1)] = tp = TypeA_(1) / a[i * (stride_a + 1)];
        }
        catch(...)
        {
          return false;
        }

        // eliminate column i of A_ji for j > i
        for(TypeSize_ j(i + 1) ; j < n ; ++j)
        {
          // calculate elimination factor: L_ji := A_ji / U_ii
          a[j * stride_a + i] *= tp;

          // eliminate A_ji
          vec_axpy(n - i - 1, &a[j * stride_a + i + 1], &a[i * (stride_a + 1) + 1], -a[j * stride_a + i]);
        }
      }

      return true;
    }

    /// \cond internal
    namespace Intern
    {
      template<bool>
      struct MatVecSolve;

      template<>
      struct MatVecSolve<false>
      {
        template<
          typename TypeA_,
          typename TypeX_,
          typename TypeP_,
          typename TypeSize_>
        static void apply(
          const TypeSize_ n,
          TypeX_ x[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeP_ p[])
        {
          CONTEXT("LinAlg::Intern::MatVecSolve<false>::apply()");
          // solve P*(L+I) * y = b
          for(TypeSize_ i(0) ; i < n ; ++i)
          {
            // apply pivoting
            if(p[i] != TypeP_(i))
            {
              swap(x[i], x[p[i]]);
            }

            // eliminate row i of L
            size_t k = i * stride_a;
            for(TypeSize_ j(0) ; j < i ; ++j)
            {
              x[i] -= TypeX_(a[k + j]) * x[j];
            }
          }

          // solve (U+D) * x = y
          for(TypeSize_ i(n) ; i > 0 ; )
          {
            TypeSize_ k = (--i) * stride_a;
            for(TypeSize_ j(i + 1) ; j < n ; ++j)
            {
              x[i] -= TypeX_(a[k + j]) * x[j];
            }
            x[i] *= TypeX_(a[k + i]);
          }
        }
      }; // MatVecSolve<false>

      template<>
      struct MatVecSolve<true>
      {
        template<
          typename TypeA_,
          typename TypeX_,
          typename TypeP_,
          typename TypeSize_>
        static void apply(
          const TypeSize_ n,
          TypeX_ x[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeP_ p[])
        {
          CONTEXT("LinAlg::Intern::MatVecSolve<true>::apply()");
          // solve (U^T+D) * y = b
          for(TypeSize_ i(0) ; i < n ; ++i)
          {
            TypeSize_ k = i * stride_a;
            x[i] *= TypeX_(a[k + i]);
            for(TypeSize_ j(i + 1) ; j < n ; ++j)
            {
              x[j] -= TypeX_(a[k + j]) * x[i];
            }
          }

          // solve (L^T+I) * P^T * x = y
          for(TypeSize_ i(n) ; i > 0 ; )
          {
            // eliminate row i of L
            TypeSize_ k = (--i) * stride_a;
            for(TypeSize_ j(0) ; j < i ; ++j)
            {
              x[j] -= TypeX_(a[k + j]) * x[i];
            }

            // apply pivoting
            if(p[i] != TypeP_(i))
            {
              swap(x[i], x[p[i]]);
            }
          }
        }
      }; // MatVecSolve<true>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Solves a linear system with a right-hand-side vector.
     *
     * This function solves a linear system <em>Ax = b</em> where \em x and \em b are vectors of length \em n and
     * where the n-by-n system matrix \em A is given by its PLU-factorisation calculated by the #mat_factorise
     * function.
     *
     * \param[in] n
     * The number of rows/columns of the system matrix \e A.
     *
     * \param[in,out] x
     * On entry, the vector containing the right-hand-side \e b of the linear system.\n
     * On exit, the vector containing the solution \e x of the linear system.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in] a
     * The PLU factorisation of the system matrix, as returned by #mat_factorise.
     *
     * \param[in] p
     * The pivot array of the PLU factorisation, as returned by #mat_factorise.
     *
     * \see mat_solve_mat
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      typename TypeA_,
      typename TypeX_,
      typename TypeP_,
      typename TypeSize_>
    inline void mat_solve_vec(
      const TypeSize_ n,
      TypeX_ x[],
      const TypeSize_ stride_a,
      const TypeA_ a[],
      const TypeP_ p[])
    {
      CONTEXT("LinAlg::mat_solve_vec()");
      ASSERT_(stride_a >= n);
      ASSERT_(x != nullptr);
      ASSERT_(a != nullptr);
      ASSERT_(p != nullptr);

      Intern::MatVecSolve<trans_a_>::apply(n, x, stride_a, a, p);
    }

    /// \cond internal
    namespace Intern
    {
      template<bool>
      struct MatMatSolve;

      template<>
      struct MatMatSolve<false>
      {
        template<
          typename TypeA_,
          typename TypeX_,
          typename TypeP_,
          typename TypeSize_>
        static void apply(
          const TypeSize_ n,
          const TypeSize_ m,
          const TypeSize_ stride_x,
          TypeX_ x[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeP_ p[])
        {
          CONTEXT("LinAlg::Intern::MatMatSolve<false>::apply()");

          // solve P*(L+I)*Y = B
          for(TypeSize_ i(0) ; i < n ; ++i)
          {
            // apply pivoting
            if(p[i] != TypeP_(i))
            {
              vec_swap(m, &x[i * stride_x], &x[p[i] * stride_x]);
            }

            // eliminate row i of L
            for(TypeSize_ j(0) ; j < i ; ++j)
            {
              vec_axpy(m, &x[i * stride_x], &x[j * stride_x], -TypeX_(a[i * stride_a + j]));
            }
          }

          // solve (U+D)*X = Y
          for(TypeSize_ i(n) ; i > 0 ; )
          {
            --i;
            for(TypeSize_ j(i + 1) ; j < n ; ++j)
            {
              vec_axpy(m, &x[i * stride_x], &x[j * stride_x], -TypeX_(a[i * stride_a + j]));
            }
            vec_scale(m, &x[i * stride_x], TypeX_(a[i * stride_a + i]));
          }
        }
      }; // MatMatSolve<false>

      template<>
      struct MatMatSolve<true>
      {
        template<
          typename TypeA_,
          typename TypeX_,
          typename TypeP_,
          typename TypeSize_>
        static void apply(
          const TypeSize_ n,
          const TypeSize_ m,
          const TypeSize_ stride_x,
          TypeX_ x[],
          const TypeSize_ stride_a,
          const TypeA_ a[],
          const TypeP_ p[])
        {
          CONTEXT("LinAlg::Intern::MatMatSolve<true>::apply()");

          // solve (U^T+D) * Y = B
          for(TypeSize_ i(0) ; i < n ; ++i)
          {
            vec_scale(m, &x[i * stride_x], TypeX_(a[i * stride_a + i]));
            for(TypeSize_ j = i + 1 ; j < n ; ++j)
            {
              vec_axpy(m, &x[j * stride_x], &x[i * stride_x], -TypeX_(a[i * stride_a + j]));
            }
          }

          // solve (L^T+I) * P^T * Y = B
          for(TypeSize_ i(n) ; i > 0 ; )
          {
            --i;
            // eliminate row i of L
            for(TypeSize_ j(0) ; j < i ; ++j)
            {
              vec_axpy(m, &x[j * stride_x], &x[i * stride_x], -TypeX_(a[i * stride_a + j]));
            }

            // apply pivoting
            if(p[i] != TypeP_(i))
            {
              vec_swap(m, &x[i * stride_x], &x[p[i] * stride_x]);
            }
          }
        }
      }; // MatMatSolve<true>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Solves a linear system with a right-hand-side matrix.
     *
     * This function solves a linear system <em>AX = B</em> where \em X and \em B are n-by-m matrices and where
     * the n-by-n system matrix \em A is given by its PLU-factorisation calculated by the #mat_factorise function.
     *
     * \param[in] n
     * The number of rows/columns of the system matrix \e A.
     *
     * \param[in] m
     * The number of columns of the right-hand-side matrix \e X.
     *
     * \param[in] stride_x
     * The stride of \e X.
     *
     * \param[in,out] x
     * On entry, the matrix containing the right-hand-side \e B of the linear system.\n
     * On exit, the matrix containing the solution \e X of the linear system.
     *
     * \param[in] stride_a
     * The stride of \e A.
     *
     * \param[in] a
     * The PLU factorisation of the system matrix, as returned by #mat_factorise.
     *
     * \param[in] p
     * The pivot array of the PLU factorisation, as returned by #mat_factorise.
     *
     * \see mat_solve_vec
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      typename TypeA_,
      typename TypeX_,
      typename TypeP_,
      typename TypeSize_>
    inline void mat_solve_mat(
      const TypeSize_ n,
      const TypeSize_ m,
      const TypeSize_ stride_x,
      TypeX_ x[],
      const TypeSize_ stride_a,
      const TypeA_ a[],
      const TypeP_ p[])
    {
      CONTEXT("LinAlg::mat_solve_mat()");
      ASSERT_(stride_x >= m);
      ASSERT_(stride_a >= n);
      ASSERT_(x != nullptr);
      ASSERT_(a != nullptr);
      ASSERT_(p != nullptr);

      Intern::MatMatSolve<trans_a_>::apply(n, m, stride_x, x, stride_a, a, p);
    }

    /// \cond internal
    namespace Intern
    {
      // CSR mat-vec-mult helper class
      template<bool>
      struct CSRVecMult;

      // CSR mat-vec-mult: y += alpha * A * x
      template<>
      struct CSRVecMult<false>
      {
        template<
          typename TypePtr_,
          typename TypeIdx_,
          typename TypeA_,
          typename TypeX_,
          typename TypeY_,
          typename TypeSize_>
        inline static void apply(
          TypeSize_ n,
          TypeY_ y[],
          const TypePtr_ row_ptr[],
          const TypePtr_ row_end[],
          const TypeIdx_ col_idx[],
          const TypeA_ a[],
          const TypeX_ x[],
          const TypeY_ alpha)
        {
          CONTEXT("LinAlg::Intern::CSRVecMult<false>::apply()");
          for(TypeSize_ i(0) ; i < n; ++i)
          {
            TypeY_ t = TypeY_(0);
            for(TypePtr_ j(row_ptr[i]) ; j < row_end[i] ; ++j)
            {
              t += TypeY_(a[j]) * TypeY_(x[col_idx[j]]);
            }
            y[i] += alpha * t;
          }
        }
      };

      // CSR mat-vec-mult: y += alpha * A^T * x
      template<>
      struct CSRVecMult<true>
      {
        template<
          typename TypePtr_,
          typename TypeIdx_,
          typename TypeA_,
          typename TypeX_,
          typename TypeY_,
          typename TypeSize_>
        inline static void apply(
          TypeSize_ n,
          TypeY_ y[],
          const TypePtr_ row_ptr[],
          const TypePtr_ row_end[],
          const TypeIdx_ col_idx[],
          const TypeA_ a[],
          const TypeX_ x[],
          const TypeY_ alpha)
        {
          CONTEXT("LinAlg::Intern::CSRVecMult<true>::apply()");
          for(TypeSize_ i(0) ; i < n; ++i)
          {
            TypeY_ t = alpha * TypeY_(x[i]);
            for(TypePtr_ j(row_ptr[i]) ; j < row_end[i] ; ++j)
            {
              y[col_idx[j]] += TypeY_(a[j]) * t;
            }
          }
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Performs a matrix-vector-multiplication with a CSR matrix.
     *
     * \param[in] n
     * The number of rows of the sparse matrix.
     *
     * \param[in,out] y
     * The vector that recieves the result of the operation.
     *
     * \param[in] row_ptr
     * The row-pointer array of the sparse matrix.
     *
     * \param[in] row_end
     * The row-end-pointer array of the sparse matrix. May be set to \c nullptr if no row-end-pointer array is
     * available (classic CSR format).
     *
     * \param[in] col_idx
     * The column-index array of the sparse matrix.
     *
     * \param[in] a
     * The data array of the sparse matrix.
     *
     * \param[in] x
     * The vector that is to be multiplied by the sparse matrix.
     *
     * \param[in] alpha
     * The scaling factor for the product.
     *
     * \author Peter Zajac
     */
    template<
      bool trans_a_,
      typename TypePtr_,
      typename TypeIdx_,
      typename TypeA_,
      typename TypeX_,
      typename TypeY_,
      typename TypeSize_>
    inline void csr_vec_mult(
      TypeSize_ n,
      TypeY_ y[],
      const TypePtr_ row_ptr[],
      const TypePtr_ row_end[],
      const TypeIdx_ col_idx[],
      const TypeA_ a[],
      const TypeX_ x[],
      const TypeY_ alpha = TypeY_(1))
    {
      CONTEXT("LinAlg::csr_vec_mult()");
      ASSERT_(y != nullptr);
      ASSERT_(row_ptr != nullptr);
      ASSERT_(col_idx != nullptr);
      ASSERT_(a != nullptr);
      ASSERT_(x != nullptr);

      if(row_end == nullptr)
        row_end = &row_ptr[1];

      Intern::CSRVecMult<trans_a_>::apply(n, y, row_ptr, row_end, col_idx, a, x, alpha);
    }

  } // namespace LinAlg
} // namespace FEAST

#endif // KERNEL_LINEAR_ALGEBRA_HPP

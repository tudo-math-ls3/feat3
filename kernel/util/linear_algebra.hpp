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

// includes, system
#include <complex>
#include <cmath>

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

      for(TypeSize_ i(0) ; i < n ; ++i)
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
        r += std::abs(x[i]);
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
        r += std::abs(x[i]);
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

      return std::sqrt(r);
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

      return std::sqrt(r);
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
        t = std::abs(x[i]);
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

      return std::sqrt(r);
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

    template<
      int m_,
      int n_,
      typename TypeA_,
      typename TypeX_,
      typename TypeY_>
    inline void mat_vec_mult(
      TypeY_ (&y)[m_],
      const TypeA_ (&a)[m_][n_],
      const TypeX_ (&x)[n_],
      const TypeY_ alpha = TypeY_(1))
    {
      //Intern::MatVecMult<false>::apply(int(m_), int(n_), int(n_), &y[0], &a[0][0], &x[0], alpha);
      for(int i(0); i < m_; ++i)
      {
        y[i] = TypeY_(0);
        for(int j(0); j < n_; ++j)
        {
          y[i] += alpha * TypeY_(a[i][j] * x[j]);
        }
      }
    }

    template<
      int m_,
      int n_,
      typename TypeA_,
      typename TypeX_,
      typename TypeY_>
    inline void vec_mat_mult(
      TypeY_ (&y)[n_],
      const TypeA_ (&a)[m_][n_],
      const TypeX_ (&x)[m_],
      const TypeY_ alpha = TypeY_(1))
    {
      //Intern::MatVecMult<true>::apply(int(n_), int(m_), int(m_), &y[0], &a[0][0], &x[0], alpha);
      for(int i(0); i < n_; ++i)
      {
        y[i] = TypeY_(0);
        for(int j(0); j < m_; ++j)
        {
          y[i] += alpha * TypeY_(a[j][i] * x[j]);
        }
      }
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
        r = std::max(r, vec_norm_asum(n, &a[i * stride]));
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
        r = std::max(r, vec_norm_asum(n, &a[i * stride]));
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
      return std::sqrt(r);
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
      return std::sqrt(r);
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
        r = std::max(r, vec_norm_max(n, &a[i * stride]));
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
        r = std::max(r, vec_norm_max(n, &a[i * stride]));
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
      using std::abs;
      CONTEXT("LinAlg::mat_factorise()");
      ASSERT_(a != nullptr);
      ASSERT_(p != nullptr);
      ASSERT_(stride_a >= n);

      // loop over all rows of A
      for(TypeSize_ i(0) ; i < m ; ++i)
      {
        // find a pivot
        p[i] = TypeP_(i);
        TypeA_ tp = abs(a[i * (stride_a + 1)]);
        for(TypeSize_ j(i + 1) ; j < n ; ++j)
        {
          TypeA_ t = abs(a[j * stride_a + i]);
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
          a[i * (stride_a + 1)] = tp = TypeA_(1) / tp;
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
      template<
        int m_,
        int n_>
      struct MatDet;

      template<>
      struct MatDet<1,1>
      {
        template<typename M_>
        static M_ compute(const M_ m[])
        {
          return m[0];
        }
      };

      template<>
      struct MatDet<2,2>
      {
        template<typename M_>
        static M_ compute(const M_ m[])
        {
          return m[0]*m[3] - m[1]*m[2];
        }
      };

      template<>
      struct MatDet<3,3>
      {
        template<typename M_>
        static M_ compute(const M_ m[])
        {
          return
            + m[0]*(m[4]*m[8] - m[7]*m[5])
            + m[3]*(m[7]*m[2] - m[1]*m[8])
            + m[6]*(m[1]*m[5] - m[4]*m[2]);
        }
      };

      template<>
      struct MatDet<4,4>
      {
        template<typename M_>
        static M_ compute(const M_ m[])
        {
          // 2x2 determinants of rows 3-4
          M_ w[6] =
          {
            m[ 8]*m[13] - m[ 9]*m[12],
            m[ 8]*m[14] - m[10]*m[12],
            m[ 8]*m[15] - m[11]*m[12],
            m[ 9]*m[14] - m[10]*m[13],
            m[ 9]*m[15] - m[11]*m[13],
            m[10]*m[15] - m[11]*m[14]
          };

          return
            + m[0] * (m[5]*w[5] - m[6]*w[4] + m[7]*w[3])
            - m[1] * (m[4]*w[5] - m[6]*w[2] + m[7]*w[1])
            + m[2] * (m[4]*w[4] - m[5]*w[2] + m[7]*w[0])
            - m[3] * (m[4]*w[3] - m[5]*w[1] + m[6]*w[0]);
        }
      };

      template<>
      struct MatDet<5,5>
      {
        template<typename M_>
        static M_ compute(const M_ m[])
        {
          // 2x2 determinants of rows 4-5
          M_ v[10] =
          {
            m[15]*m[21] - m[16]*m[20],
            m[15]*m[22] - m[17]*m[20],
            m[15]*m[23] - m[18]*m[20],
            m[15]*m[24] - m[19]*m[20],
            m[16]*m[22] - m[17]*m[21],
            m[16]*m[23] - m[18]*m[21],
            m[16]*m[24] - m[19]*m[21],
            m[17]*m[23] - m[18]*m[22],
            m[17]*m[24] - m[19]*m[22],
            m[18]*m[24] - m[19]*m[23]
          };
          // 3x3 determinants of rows 3-4-5
          M_ w[10] =
          {
            m[10]*v[4] - m[11]*v[1] + m[12]*v[0],
            m[10]*v[5] - m[11]*v[2] + m[13]*v[0],
            m[10]*v[6] - m[11]*v[3] + m[14]*v[0],
            m[10]*v[7] - m[12]*v[2] + m[13]*v[1],
            m[10]*v[8] - m[12]*v[3] + m[14]*v[1],
            m[10]*v[9] - m[13]*v[3] + m[14]*v[2],
            m[11]*v[7] - m[12]*v[5] + m[13]*v[4],
            m[11]*v[8] - m[12]*v[6] + m[14]*v[4],
            m[11]*v[9] - m[13]*v[6] + m[14]*v[5],
            m[12]*v[9] - m[13]*v[8] + m[14]*v[7]
          };

          return
            + m[0]*(m[6]*w[9] - m[7]*w[8] + m[8]*w[7] - m[9]*w[6])
            - m[1]*(m[5]*w[9] - m[7]*w[5] + m[8]*w[4] - m[9]*w[3])
            + m[2]*(m[5]*w[8] - m[6]*w[5] + m[8]*w[2] - m[9]*w[1])
            - m[3]*(m[5]*w[7] - m[6]*w[4] + m[7]*w[2] - m[9]*w[0])
            + m[4]*(m[5]*w[6] - m[6]*w[3] + m[7]*w[1] - m[8]*w[0]);
        }
      };

      template<>
      struct MatDet<6,6>
      {
        template<typename M_>
        static M_ compute(const M_ m[])
        {
          // 2x2 determinants of rows 5-6
          M_ v[15] =
          {
            m[24]*m[31] - m[25]*m[30],
            m[24]*m[32] - m[26]*m[30],
            m[24]*m[33] - m[27]*m[30],
            m[24]*m[34] - m[28]*m[30],
            m[24]*m[35] - m[29]*m[30],
            m[25]*m[32] - m[26]*m[31],
            m[25]*m[33] - m[27]*m[31],
            m[25]*m[34] - m[28]*m[31],
            m[25]*m[35] - m[29]*m[31],
            m[26]*m[33] - m[27]*m[32],
            m[26]*m[34] - m[28]*m[32],
            m[26]*m[35] - m[29]*m[32],
            m[27]*m[34] - m[28]*m[33],
            m[27]*m[35] - m[29]*m[33],
            m[28]*m[35] - m[29]*m[34]
          };
          // 3x3 determinants of rows 4-5-6
          M_ w[20] =
          {
            m[18]*v[ 5] - m[19]*v[ 1] + m[20]*v[ 0],
            m[18]*v[ 6] - m[19]*v[ 2] + m[21]*v[ 0],
            m[18]*v[ 7] - m[19]*v[ 3] + m[22]*v[ 0],
            m[18]*v[ 8] - m[19]*v[ 4] + m[23]*v[ 0],
            m[18]*v[ 9] - m[20]*v[ 2] + m[21]*v[ 1],
            m[18]*v[10] - m[20]*v[ 3] + m[22]*v[ 1],
            m[18]*v[11] - m[20]*v[ 4] + m[23]*v[ 1],
            m[18]*v[12] - m[21]*v[ 3] + m[22]*v[ 2],
            m[18]*v[13] - m[21]*v[ 4] + m[23]*v[ 2],
            m[18]*v[14] - m[22]*v[ 4] + m[23]*v[ 3],
            m[19]*v[ 9] - m[20]*v[ 6] + m[21]*v[ 5],
            m[19]*v[10] - m[20]*v[ 7] + m[22]*v[ 5],
            m[19]*v[11] - m[20]*v[ 8] + m[23]*v[ 5],
            m[19]*v[12] - m[21]*v[ 7] + m[22]*v[ 6],
            m[19]*v[13] - m[21]*v[ 8] + m[23]*v[ 6],
            m[19]*v[14] - m[22]*v[ 8] + m[23]*v[ 7],
            m[20]*v[12] - m[21]*v[10] + m[22]*v[ 9],
            m[20]*v[13] - m[21]*v[11] + m[23]*v[ 9],
            m[20]*v[14] - m[22]*v[11] + m[23]*v[10],
            m[21]*v[14] - m[22]*v[13] + m[23]*v[12]
          };
          // 4x4 determinants of rows 3-4-5-6
          v[ 0] = m[12]*w[10] - m[13]*w[ 4] + m[14]*w[ 1] - m[15]*w[ 0];
          v[ 1] = m[12]*w[11] - m[13]*w[ 5] + m[14]*w[ 2] - m[16]*w[ 0];
          v[ 2] = m[12]*w[12] - m[13]*w[ 6] + m[14]*w[ 3] - m[17]*w[ 0];
          v[ 3] = m[12]*w[13] - m[13]*w[ 7] + m[15]*w[ 2] - m[16]*w[ 1];
          v[ 4] = m[12]*w[14] - m[13]*w[ 8] + m[15]*w[ 3] - m[17]*w[ 1];
          v[ 5] = m[12]*w[15] - m[13]*w[ 9] + m[16]*w[ 3] - m[17]*w[ 2];
          v[ 6] = m[12]*w[16] - m[14]*w[ 7] + m[15]*w[ 5] - m[16]*w[ 4];
          v[ 7] = m[12]*w[17] - m[14]*w[ 8] + m[15]*w[ 6] - m[17]*w[ 4];
          v[ 8] = m[12]*w[18] - m[14]*w[ 9] + m[16]*w[ 6] - m[17]*w[ 5];
          v[ 9] = m[12]*w[19] - m[15]*w[ 9] + m[16]*w[ 8] - m[17]*w[ 7];
          v[10] = m[13]*w[16] - m[14]*w[13] + m[15]*w[11] - m[16]*w[10];
          v[11] = m[13]*w[17] - m[14]*w[14] + m[15]*w[12] - m[17]*w[10];
          v[12] = m[13]*w[18] - m[14]*w[15] + m[16]*w[12] - m[17]*w[11];
          v[13] = m[13]*w[19] - m[15]*w[15] + m[16]*w[14] - m[17]*w[13];
          v[14] = m[14]*w[19] - m[15]*w[18] + m[16]*w[17] - m[17]*w[16];

          return
            + m[0]*(m[ 7]*v[14] - m[ 8]*v[13] + m[ 9]*v[12] - m[10]*v[11] + m[11]*v[10])
            - m[1]*(m[ 6]*v[14] - m[ 8]*v[ 9] + m[ 9]*v[ 8] - m[10]*v[ 7] + m[11]*v[ 6])
            + m[2]*(m[ 6]*v[13] - m[ 7]*v[ 9] + m[ 9]*v[ 5] - m[10]*v[ 4] + m[11]*v[ 3])
            - m[3]*(m[ 6]*v[12] - m[ 7]*v[ 8] + m[ 8]*v[ 5] - m[10]*v[ 2] + m[11]*v[ 1])
            + m[4]*(m[ 6]*v[11] - m[ 7]*v[ 7] + m[ 8]*v[ 4] - m[ 9]*v[ 2] + m[11]*v[ 0])
            - m[5]*(m[ 6]*v[10] - m[ 7]*v[ 6] + m[ 8]*v[ 3] - m[ 9]*v[ 1] + m[10]*v[ 0]);
        }
      };

      template<>
      struct MatDet<2, 1>
      {
        template<typename M_>
        static M_ compute(const M_ m[])
        {
          return std::sqrt(sqr(m[0]) + sqr(m[1]));
        }
      };

      template<>
      struct MatDet<3, 1>
      {
        template<typename M_>
        static M_ compute(const M_& m)
        {
          return std::sqrt(sqr(m[0]) + sqr(m[1]) + sqr(m[2]));
        }
      };

      template<>
      struct MatDet<3, 2>
      {
        template<typename M_>
        static M_ compute(const M_& m)
        {
          // This is the euclid norm of the 3D cross product of the two matrix columns.
          return std::sqrt(
            sqr(m[2]*m[5] - m[4]*m[3]) +
            sqr(m[4]*m[1] - m[0]*m[5]) +
            sqr(m[0]*m[3] - m[2]*m[1]));
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Calculates the (generalised) matrix determinant.
     *
     * For square matrices, this function computes the determinant of a matrix.\n
     * For non-square matrices, this function computes the generalised determinant of a matrix, i.e.
     * \f[ \sqrt{\det(A^\top \cdot A)} \f]
     *
     * Currently, this function supports only the following matrix dimensions:
     * \li Determinant: <c>n</c>-by-<c>n</c> for <c>1 <= n <= 6</c>
     * \li Generalised Determinant: 2x1, 3x1, 3x2
     *
     * \tparam m_, n_
     * The dimensions of the matrix whose determinant is to be computed.
     *
     * \param[in] a
     * The matrix whose (generalised) determinant is to be computed.
     *
     * \returns
     * The (generalised) determinant of \p a.
     *
     * \author Peter Zajac
     */
    template<
      int m_,
      int n_,
      typename TypeA_>
    inline TypeA_ mat_det(const TypeA_ a[])
    {
      CONTEXT("LinAlg::mat_det()");
      return Intern::MatDet<m_, n_>::compute(a);
    }

    /** \copydoc mat_det(const TypeA_ a[]) */
    template<
      int m_,
      int n_,
      typename TypeA_>
    inline TypeA_ mat_det(const TypeA_ (&a)[m_][n_])
    {
      CONTEXT("LinAlg::mat_det()");
      return Intern::MatDet<m_,n_>::compute(&a[0][0]);
    }

    /// \cond internal
    namespace Intern
    {
      // matrix inversion helper class
      template<
        bool force_plu_,
        int m_,
        int n_>
      struct MatInvert
      {
        template<typename Type_>
        inline static void apply(
          Type_ b[],
          const Type_ a[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<force_plu_,m_,n_>::apply()");
          throw InternalError("Cannot invert non-square matrix!");
        }
      };

      // note: setting force_plu_u = true will always call the factorisation code below rather than calling the
      //       direct inversion variants following below.
      template<
        bool force_plu_,
        int n_>
      struct MatInvert<force_plu_, n_, n_>
      {
        template<typename Type_>
        inline static void apply(
          Type_ b[],
          const Type_ a[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<force_plu_,n_,n_>::apply()");

          // This is a fallback implementation that uses the PLU factorisation.
          // This may be replaced by something more efficient later...

          // pivot array
          size_t p[n_];

          // temporary matrix for LU decomposition
          Type_ lu[n_ * n_];

          // copy A to LU
          mat_copy(n_, n_, n_, n_, lu, a);

          // factorise LU
          mat_factorise(n_, n_, n_, lu, p);

          // set B to identity
          mat_identity(n_, n_, b);

          // solve LU*X = B
          mat_solve_mat<false>(n_, n_, n_, b, n_, lu, p);
        }
      };

      template<>
      struct MatInvert<false, 1, 1>
      {
        // direct inversion for 1x1
        template<typename Type_>
        inline static void apply(
          Type_ B[],
          const Type_ A[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<false,1,1>::apply()");
          B[0] = Type_(1) / A[0];
        }
      };

      template<>
      struct MatInvert<false, 2, 2>
      {
        // direct inversion for 2x2
        template<typename Type_>
        inline static void apply(
          Type_ B[],
          const Type_ A[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<false,2,2>::apply()");
          Type_ d = Type_(1) / (A[0] * A[3] - A[1] * A[2]);
          B[0] =  d * A[3];
          B[1] = -d * A[1];
          B[2] = -d * A[2];
          B[3] =  d * A[0];
        }
      };

      template<>
      struct MatInvert<false, 3, 3>
      {
        // direct inversion for 3x3
        template<typename Type_>
        inline static void apply(
          Type_ B[],
          const Type_ A[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<false,3,3>::apply()");
          B[0] = A[4] * A[8] - A[5] * A[7];
          B[3] = A[5] * A[6] - A[3] * A[8];
          B[6] = A[3] * A[7] - A[4] * A[6];
          Type_ d = Type_(1) / (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
          B[0] *= d;
          B[3] *= d;
          B[6] *= d;
          B[1] = d * (A[2] * A[7] - A[1] * A[8]);
          B[4] = d * (A[0] * A[8] - A[2] * A[6]);
          B[7] = d * (A[1] * A[6] - A[0] * A[7]);
          B[2] = d * (A[1] * A[5] - A[2] * A[4]);
          B[5] = d * (A[2] * A[3] - A[0] * A[5]);
          B[8] = d * (A[0] * A[4] - A[1] * A[3]);
        }
      };

      template<>
      struct MatInvert<false, 4, 4>
      {
        // direct inversion for 4x4
        template<typename Type_>
        static void apply(
          Type_ B[],
          const Type_ A[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<false,4,4>::apply()");
          Type_ W[6];
          W[0] = A[8]*A[13]-A[9]*A[12];
          W[1] = A[8]*A[14]-A[10]*A[12];
          W[2] = A[8]*A[15]-A[11]*A[12];
          W[3] = A[9]*A[14]-A[10]*A[13];
          W[4] = A[9]*A[15]-A[11]*A[13];
          W[5] = A[10]*A[15]-A[11]*A[14];
          B[ 0] = A[5]*W[5]-A[6]*W[4]+A[7]*W[3];
          B[ 4] =-A[4]*W[5]+A[6]*W[2]-A[7]*W[1];
          B[ 8] = A[4]*W[4]-A[5]*W[2]+A[7]*W[0];
          B[12] =-A[4]*W[3]+A[5]*W[1]-A[6]*W[0];
          Type_ d = Type_(1) / (A[0]*B[0]+A[1]*B[4]+A[2]*B[8]+A[3]*B[12]);
          B[ 0] *= d;
          B[ 4] *= d;
          B[ 8] *= d;
          B[12] *= d;
          B[ 1] = d*(-A[1]*W[5]+A[2]*W[4]-A[3]*W[3]);
          B[ 5] = d*( A[0]*W[5]-A[2]*W[2]+A[3]*W[1]);
          B[ 9] = d*(-A[0]*W[4]+A[1]*W[2]-A[3]*W[0]);
          B[13] = d*( A[0]*W[3]-A[1]*W[1]+A[2]*W[0]);
          W[0] = A[0]*A[5]-A[1]*A[4];
          W[1] = A[0]*A[6]-A[2]*A[4];
          W[2] = A[0]*A[7]-A[3]*A[4];
          W[3] = A[1]*A[6]-A[2]*A[5];
          W[4] = A[1]*A[7]-A[3]*A[5];
          W[5] = A[2]*A[7]-A[3]*A[6];
          B[ 2] = d*( A[13]*W[5]-A[14]*W[4]+A[15]*W[3]);
          B[ 6] = d*(-A[12]*W[5]+A[14]*W[2]-A[15]*W[1]);
          B[10] = d*( A[12]*W[4]-A[13]*W[2]+A[15]*W[0]);
          B[14] = d*(-A[12]*W[3]+A[13]*W[1]-A[14]*W[0]);
          B[ 3] = d*(-A[9]*W[5]+A[10]*W[4]-A[11]*W[3]);
          B[ 7] = d*( A[8]*W[5]-A[10]*W[2]+A[11]*W[1]);
          B[11] = d*(-A[8]*W[4]+A[9]*W[2]-A[11]*W[0]);
          B[15] = d*( A[8]*W[3]-A[9]*W[1]+A[10]*W[0]);
        }
      };

      template<>
      struct MatInvert<false, 5, 5>
      {
        // direct inversion for 5x5
        template<typename Type_>
        static void apply(
          Type_ B[],
          const Type_ A[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<false,5,5>::apply()");
          Type_ W[20];
          W[ 0] = A[15]*A[21]-A[16]*A[20];
          W[ 1] = A[15]*A[22]-A[17]*A[20];
          W[ 2] = A[15]*A[23]-A[18]*A[20];
          W[ 3] = A[15]*A[24]-A[19]*A[20];
          W[ 4] = A[16]*A[22]-A[17]*A[21];
          W[ 5] = A[16]*A[23]-A[18]*A[21];
          W[ 6] = A[16]*A[24]-A[19]*A[21];
          W[ 7] = A[17]*A[23]-A[18]*A[22];
          W[ 8] = A[17]*A[24]-A[19]*A[22];
          W[ 9] = A[18]*A[24]-A[19]*A[23];
          W[10] = A[10]*W[4]-A[11]*W[1]+A[12]*W[0];
          W[11] = A[10]*W[5]-A[11]*W[2]+A[13]*W[0];
          W[12] = A[10]*W[6]-A[11]*W[3]+A[14]*W[0];
          W[13] = A[10]*W[7]-A[12]*W[2]+A[13]*W[1];
          W[14] = A[10]*W[8]-A[12]*W[3]+A[14]*W[1];
          W[15] = A[10]*W[9]-A[13]*W[3]+A[14]*W[2];
          W[16] = A[11]*W[7]-A[12]*W[5]+A[13]*W[4];
          W[17] = A[11]*W[8]-A[12]*W[6]+A[14]*W[4];
          W[18] = A[11]*W[9]-A[13]*W[6]+A[14]*W[5];
          W[19] = A[12]*W[9]-A[13]*W[8]+A[14]*W[7];
          B[ 0] = A[6]*W[19]-A[7]*W[18]+A[8]*W[17]-A[9]*W[16];
          B[ 5] =-A[5]*W[19]+A[7]*W[15]-A[8]*W[14]+A[9]*W[13];
          B[10] = A[5]*W[18]-A[6]*W[15]+A[8]*W[12]-A[9]*W[11];
          B[15] =-A[5]*W[17]+A[6]*W[14]-A[7]*W[12]+A[9]*W[10];
          B[20] = A[5]*W[16]-A[6]*W[13]+A[7]*W[11]-A[8]*W[10];
          Type_ d = Type_(1) / (A[0]*B[0]+A[1]*B[5]+A[2]*B[10]+A[3]*B[15]+A[4]*B[20]);
          B[ 0] *= d;
          B[ 5] *= d;
          B[10] *= d;
          B[15] *= d;
          B[20] *= d;
          B[ 1] = d*(-A[1]*W[19]+A[2]*W[18]-A[3]*W[17]+A[4]*W[16]);
          B[ 6] = d*( A[0]*W[19]-A[2]*W[15]+A[3]*W[14]-A[4]*W[13]);
          B[11] = d*(-A[0]*W[18]+A[1]*W[15]-A[3]*W[12]+A[4]*W[11]);
          B[16] = d*( A[0]*W[17]-A[1]*W[14]+A[2]*W[12]-A[4]*W[10]);
          B[21] = d*(-A[0]*W[16]+A[1]*W[13]-A[2]*W[11]+A[3]*W[10]);
          W[10] = A[5]*W[4]-A[6]*W[1]+A[7]*W[0];
          W[11] = A[5]*W[5]-A[6]*W[2]+A[8]*W[0];
          W[12] = A[5]*W[6]-A[6]*W[3]+A[9]*W[0];
          W[13] = A[5]*W[7]-A[7]*W[2]+A[8]*W[1];
          W[14] = A[5]*W[8]-A[7]*W[3]+A[9]*W[1];
          W[15] = A[5]*W[9]-A[8]*W[3]+A[9]*W[2];
          W[16] = A[6]*W[7]-A[7]*W[5]+A[8]*W[4];
          W[17] = A[6]*W[8]-A[7]*W[6]+A[9]*W[4];
          W[18] = A[6]*W[9]-A[8]*W[6]+A[9]*W[5];
          W[19] = A[7]*W[9]-A[8]*W[8]+A[9]*W[7];
          B[ 2] = d*( A[1]*W[19]-A[2]*W[18]+A[3]*W[17]-A[4]*W[16]);
          B[ 7] = d*(-A[0]*W[19]+A[2]*W[15]-A[3]*W[14]+A[4]*W[13]);
          B[12] = d*( A[0]*W[18]-A[1]*W[15]+A[3]*W[12]-A[4]*W[11]);
          B[17] = d*(-A[0]*W[17]+A[1]*W[14]-A[2]*W[12]+A[4]*W[10]);
          B[22] = d*( A[0]*W[16]-A[1]*W[13]+A[2]*W[11]-A[3]*W[10]);
          W[ 0] = A[0]*A[6]-A[1]*A[5];
          W[ 1] = A[0]*A[7]-A[2]*A[5];
          W[ 2] = A[0]*A[8]-A[3]*A[5];
          W[ 3] = A[0]*A[9]-A[4]*A[5];
          W[ 4] = A[1]*A[7]-A[2]*A[6];
          W[ 5] = A[1]*A[8]-A[3]*A[6];
          W[ 6] = A[1]*A[9]-A[4]*A[6];
          W[ 7] = A[2]*A[8]-A[3]*A[7];
          W[ 8] = A[2]*A[9]-A[4]*A[7];
          W[ 9] = A[3]*A[9]-A[4]*A[8];
          W[10] = A[10]*W[4]-A[11]*W[1]+A[12]*W[0];
          W[11] = A[10]*W[5]-A[11]*W[2]+A[13]*W[0];
          W[12] = A[10]*W[6]-A[11]*W[3]+A[14]*W[0];
          W[13] = A[10]*W[7]-A[12]*W[2]+A[13]*W[1];
          W[14] = A[10]*W[8]-A[12]*W[3]+A[14]*W[1];
          W[15] = A[10]*W[9]-A[13]*W[3]+A[14]*W[2];
          W[16] = A[11]*W[7]-A[12]*W[5]+A[13]*W[4];
          W[17] = A[11]*W[8]-A[12]*W[6]+A[14]*W[4];
          W[18] = A[11]*W[9]-A[13]*W[6]+A[14]*W[5];
          W[19] = A[12]*W[9]-A[13]*W[8]+A[14]*W[7];
          B[ 3] = d*( A[21]*W[19]-A[22]*W[18]+A[23]*W[17]-A[24]*W[16]);
          B[ 8] = d*(-A[20]*W[19]+A[22]*W[15]-A[23]*W[14]+A[24]*W[13]);
          B[13] = d*( A[20]*W[18]-A[21]*W[15]+A[23]*W[12]-A[24]*W[11]);
          B[18] = d*(-A[20]*W[17]+A[21]*W[14]-A[22]*W[12]+A[24]*W[10]);
          B[23] = d*( A[20]*W[16]-A[21]*W[13]+A[22]*W[11]-A[23]*W[10]);
          B[ 4] = d*(-A[16]*W[19]+A[17]*W[18]-A[18]*W[17]+A[19]*W[16]);
          B[ 9] = d*( A[15]*W[19]-A[17]*W[15]+A[18]*W[14]-A[19]*W[13]);
          B[14] = d*(-A[15]*W[18]+A[16]*W[15]-A[18]*W[12]+A[19]*W[11]);
          B[19] = d*( A[15]*W[17]-A[16]*W[14]+A[17]*W[12]-A[19]*W[10]);
          B[24] = d*(-A[15]*W[16]+A[16]*W[13]-A[17]*W[11]+A[18]*W[10]);
        }
      };

      template<>
      struct MatInvert<false, 6, 6>
      {
        // direct inversion for 6x6
        template<typename Type_>
        static void apply(
          Type_ B[],
          const Type_ A[])
        {
          CONTEXT("LinAlg::Intern::MatInvert<false,6,6>::apply()");
          Type_ W[35];
          W[ 0] = A[24]*A[31]-A[25]*A[30];
          W[ 1] = A[24]*A[32]-A[26]*A[30];
          W[ 2] = A[24]*A[33]-A[27]*A[30];
          W[ 3] = A[24]*A[34]-A[28]*A[30];
          W[ 4] = A[24]*A[35]-A[29]*A[30];
          W[ 5] = A[25]*A[32]-A[26]*A[31];
          W[ 6] = A[25]*A[33]-A[27]*A[31];
          W[ 7] = A[25]*A[34]-A[28]*A[31];
          W[ 8] = A[25]*A[35]-A[29]*A[31];
          W[ 9] = A[26]*A[33]-A[27]*A[32];
          W[10] = A[26]*A[34]-A[28]*A[32];
          W[11] = A[26]*A[35]-A[29]*A[32];
          W[12] = A[27]*A[34]-A[28]*A[33];
          W[13] = A[27]*A[35]-A[29]*A[33];
          W[14] = A[28]*A[35]-A[29]*A[34];
          W[15] = A[18]*W[5]-A[19]*W[1]+A[20]*W[0];
          W[16] = A[18]*W[6]-A[19]*W[2]+A[21]*W[0];
          W[17] = A[18]*W[7]-A[19]*W[3]+A[22]*W[0];
          W[18] = A[18]*W[8]-A[19]*W[4]+A[23]*W[0];
          W[19] = A[18]*W[9]-A[20]*W[2]+A[21]*W[1];
          W[20] = A[18]*W[10]-A[20]*W[3]+A[22]*W[1];
          W[21] = A[18]*W[11]-A[20]*W[4]+A[23]*W[1];
          W[22] = A[18]*W[12]-A[21]*W[3]+A[22]*W[2];
          W[23] = A[18]*W[13]-A[21]*W[4]+A[23]*W[2];
          W[24] = A[18]*W[14]-A[22]*W[4]+A[23]*W[3];
          W[25] = A[19]*W[9]-A[20]*W[6]+A[21]*W[5];
          W[26] = A[19]*W[10]-A[20]*W[7]+A[22]*W[5];
          W[27] = A[19]*W[11]-A[20]*W[8]+A[23]*W[5];
          W[28] = A[19]*W[12]-A[21]*W[7]+A[22]*W[6];
          W[29] = A[19]*W[13]-A[21]*W[8]+A[23]*W[6];
          W[30] = A[19]*W[14]-A[22]*W[8]+A[23]*W[7];
          W[31] = A[20]*W[12]-A[21]*W[10]+A[22]*W[9];
          W[32] = A[20]*W[13]-A[21]*W[11]+A[23]*W[9];
          W[33] = A[20]*W[14]-A[22]*W[11]+A[23]*W[10];
          W[34] = A[21]*W[14]-A[22]*W[13]+A[23]*W[12];
          W[ 0] = A[12]*W[25]-A[13]*W[19]+A[14]*W[16]-A[15]*W[15];
          W[ 1] = A[12]*W[26]-A[13]*W[20]+A[14]*W[17]-A[16]*W[15];
          W[ 2] = A[12]*W[27]-A[13]*W[21]+A[14]*W[18]-A[17]*W[15];
          W[ 3] = A[12]*W[28]-A[13]*W[22]+A[15]*W[17]-A[16]*W[16];
          W[ 4] = A[12]*W[29]-A[13]*W[23]+A[15]*W[18]-A[17]*W[16];
          W[ 5] = A[12]*W[30]-A[13]*W[24]+A[16]*W[18]-A[17]*W[17];
          W[ 6] = A[12]*W[31]-A[14]*W[22]+A[15]*W[20]-A[16]*W[19];
          W[ 7] = A[12]*W[32]-A[14]*W[23]+A[15]*W[21]-A[17]*W[19];
          W[ 8] = A[12]*W[33]-A[14]*W[24]+A[16]*W[21]-A[17]*W[20];
          W[ 9] = A[12]*W[34]-A[15]*W[24]+A[16]*W[23]-A[17]*W[22];
          W[10] = A[13]*W[31]-A[14]*W[28]+A[15]*W[26]-A[16]*W[25];
          W[11] = A[13]*W[32]-A[14]*W[29]+A[15]*W[27]-A[17]*W[25];
          W[12] = A[13]*W[33]-A[14]*W[30]+A[16]*W[27]-A[17]*W[26];
          W[13] = A[13]*W[34]-A[15]*W[30]+A[16]*W[29]-A[17]*W[28];
          W[14] = A[14]*W[34]-A[15]*W[33]+A[16]*W[32]-A[17]*W[31];
          B[ 0] =  A[7]*W[14]-A[8]*W[13]+A[9]*W[12]-A[10]*W[11]+A[11]*W[10];
          B[ 6] = -A[6]*W[14]+A[8]*W[9]-A[9]*W[8]+A[10]*W[7]-A[11]*W[6];
          B[12] =  A[6]*W[13]-A[7]*W[9]+A[9]*W[5]-A[10]*W[4]+A[11]*W[3];
          B[18] = -A[6]*W[12]+A[7]*W[8]-A[8]*W[5]+A[10]*W[2]-A[11]*W[1];
          B[24] =  A[6]*W[11]-A[7]*W[7]+A[8]*W[4]-A[9]*W[2]+A[11]*W[0];
          B[30] = -A[6]*W[10]+A[7]*W[6]-A[8]*W[3]+A[9]*W[1]-A[10]*W[0];
          Type_ d = Type_(1) / (A[0]*B[0] + A[1]*B[6] + A[2]*B[12]
              + A[3]*B[18] + A[4]*B[24] + A[5]*B[30]);
          B[ 0] *= d;
          B[ 6] *= d;
          B[12] *= d;
          B[18] *= d;
          B[24] *= d;
          B[30] *= d;
          B[ 1] = d*(-A[1]*W[14]+A[2]*W[13]-A[3]*W[12]+A[4]*W[11]-A[5]*W[10]);
          B[ 7] = d*( A[0]*W[14]-A[2]*W[9]+A[3]*W[8]-A[4]*W[7]+A[5]*W[6]);
          B[13] = d*(-A[0]*W[13]+A[1]*W[9]-A[3]*W[5]+A[4]*W[4]-A[5]*W[3]);
          B[19] = d*( A[0]*W[12]-A[1]*W[8]+A[2]*W[5]-A[4]*W[2]+A[5]*W[1]);
          B[25] = d*(-A[0]*W[11]+A[1]*W[7]-A[2]*W[4]+A[3]*W[2]-A[5]*W[0]);
          B[31] = d*( A[0]*W[10]-A[1]*W[6]+A[2]*W[3]-A[3]*W[1]+A[4]*W[0]);
          W[ 0] = A[24]*A[31]-A[25]*A[30];
          W[ 1] = A[24]*A[32]-A[26]*A[30];
          W[ 2] = A[24]*A[33]-A[27]*A[30];
          W[ 3] = A[24]*A[34]-A[28]*A[30];
          W[ 4] = A[24]*A[35]-A[29]*A[30];
          W[ 5] = A[25]*A[32]-A[26]*A[31];
          W[ 6] = A[25]*A[33]-A[27]*A[31];
          W[ 7] = A[25]*A[34]-A[28]*A[31];
          W[ 8] = A[25]*A[35]-A[29]*A[31];
          W[ 9] = A[26]*A[33]-A[27]*A[32];
          W[10] = A[26]*A[34]-A[28]*A[32];
          W[11] = A[26]*A[35]-A[29]*A[32];
          W[12] = A[27]*A[34]-A[28]*A[33];
          W[13] = A[27]*A[35]-A[29]*A[33];
          W[14] = A[28]*A[35]-A[29]*A[34];
          W[15] = A[6]*W[5]-A[7]*W[1]+A[8]*W[0];
          W[16] = A[6]*W[6]-A[7]*W[2]+A[9]*W[0];
          W[17] = A[6]*W[7]-A[7]*W[3]+A[10]*W[0];
          W[18] = A[6]*W[8]-A[7]*W[4]+A[11]*W[0];
          W[19] = A[6]*W[9]-A[8]*W[2]+A[9]*W[1];
          W[20] = A[6]*W[10]-A[8]*W[3]+A[10]*W[1];
          W[21] = A[6]*W[11]-A[8]*W[4]+A[11]*W[1];
          W[22] = A[6]*W[12]-A[9]*W[3]+A[10]*W[2];
          W[23] = A[6]*W[13]-A[9]*W[4]+A[11]*W[2];
          W[24] = A[6]*W[14]-A[10]*W[4]+A[11]*W[3];
          W[25] = A[7]*W[9]-A[8]*W[6]+A[9]*W[5];
          W[26] = A[7]*W[10]-A[8]*W[7]+A[10]*W[5];
          W[27] = A[7]*W[11]-A[8]*W[8]+A[11]*W[5];
          W[28] = A[7]*W[12]-A[9]*W[7]+A[10]*W[6];
          W[29] = A[7]*W[13]-A[9]*W[8]+A[11]*W[6];
          W[30] = A[7]*W[14]-A[10]*W[8]+A[11]*W[7];
          W[31] = A[8]*W[12]-A[9]*W[10]+A[10]*W[9];
          W[32] = A[8]*W[13]-A[9]*W[11]+A[11]*W[9];
          W[33] = A[8]*W[14]-A[10]*W[11]+A[11]*W[10];
          W[34] = A[9]*W[14]-A[10]*W[13]+A[11]*W[12];
          W[ 0] = A[0]*W[25]-A[1]*W[19]+A[2]*W[16]-A[3]*W[15];
          W[ 1] = A[0]*W[26]-A[1]*W[20]+A[2]*W[17]-A[4]*W[15];
          W[ 2] = A[0]*W[27]-A[1]*W[21]+A[2]*W[18]-A[5]*W[15];
          W[ 3] = A[0]*W[28]-A[1]*W[22]+A[3]*W[17]-A[4]*W[16];
          W[ 4] = A[0]*W[29]-A[1]*W[23]+A[3]*W[18]-A[5]*W[16];
          W[ 5] = A[0]*W[30]-A[1]*W[24]+A[4]*W[18]-A[5]*W[17];
          W[ 6] = A[0]*W[31]-A[2]*W[22]+A[3]*W[20]-A[4]*W[19];
          W[ 7] = A[0]*W[32]-A[2]*W[23]+A[3]*W[21]-A[5]*W[19];
          W[ 8] = A[0]*W[33]-A[2]*W[24]+A[4]*W[21]-A[5]*W[20];
          W[ 9] = A[0]*W[34]-A[3]*W[24]+A[4]*W[23]-A[5]*W[22];
          W[10] = A[1]*W[31]-A[2]*W[28]+A[3]*W[26]-A[4]*W[25];
          W[11] = A[1]*W[32]-A[2]*W[29]+A[3]*W[27]-A[5]*W[25];
          W[12] = A[1]*W[33]-A[2]*W[30]+A[4]*W[27]-A[5]*W[26];
          W[13] = A[1]*W[34]-A[3]*W[30]+A[4]*W[29]-A[5]*W[28];
          W[14] = A[2]*W[34]-A[3]*W[33]+A[4]*W[32]-A[5]*W[31];
          B[ 2] = d*( A[19]*W[14]-A[20]*W[13]+A[21]*W[12]-A[22]*W[11]+A[23]*W[10]);
          B[ 8] = d*(-A[18]*W[14]+A[20]*W[9]-A[21]*W[8]+A[22]*W[7]-A[23]*W[6]);
          B[14] = d*( A[18]*W[13]-A[19]*W[9]+A[21]*W[5]-A[22]*W[4]+A[23]*W[3]);
          B[20] = d*(-A[18]*W[12]+A[19]*W[8]-A[20]*W[5]+A[22]*W[2]-A[23]*W[1]);
          B[26] = d*( A[18]*W[11]-A[19]*W[7]+A[20]*W[4]-A[21]*W[2]+A[23]*W[0]);
          B[32] = d*(-A[18]*W[10]+A[19]*W[6]-A[20]*W[3]+A[21]*W[1]-A[22]*W[0]);
          B[ 3] = d*(-A[13]*W[14]+A[14]*W[13]-A[15]*W[12]+A[16]*W[11]-A[17]*W[10]);
          B[ 9] = d*( A[12]*W[14]-A[14]*W[9]+A[15]*W[8]-A[16]*W[7]+A[17]*W[6]);
          B[15] = d*(-A[12]*W[13]+A[13]*W[9]-A[15]*W[5]+A[16]*W[4]-A[17]*W[3]);
          B[21] = d*( A[12]*W[12]-A[13]*W[8]+A[14]*W[5]-A[16]*W[2]+A[17]*W[1]);
          B[27] = d*(-A[12]*W[11]+A[13]*W[7]-A[14]*W[4]+A[15]*W[2]-A[17]*W[0]);
          B[33] = d*( A[12]*W[10]-A[13]*W[6]+A[14]*W[3]-A[15]*W[1]+A[16]*W[0]);
          W[ 0] = A[12]*A[19]-A[13]*A[18];
          W[ 1] = A[12]*A[20]-A[14]*A[18];
          W[ 2] = A[12]*A[21]-A[15]*A[18];
          W[ 3] = A[12]*A[22]-A[16]*A[18];
          W[ 4] = A[12]*A[23]-A[17]*A[18];
          W[ 5] = A[13]*A[20]-A[14]*A[19];
          W[ 6] = A[13]*A[21]-A[15]*A[19];
          W[ 7] = A[13]*A[22]-A[16]*A[19];
          W[ 8] = A[13]*A[23]-A[17]*A[19];
          W[ 9] = A[14]*A[21]-A[15]*A[20];
          W[10] = A[14]*A[22]-A[16]*A[20];
          W[11] = A[14]*A[23]-A[17]*A[20];
          W[12] = A[15]*A[22]-A[16]*A[21];
          W[13] = A[15]*A[23]-A[17]*A[21];
          W[14] = A[16]*A[23]-A[17]*A[22];
          W[15] = A[6]*W[5]-A[7]*W[1]+A[8]*W[0];
          W[16] = A[6]*W[6]-A[7]*W[2]+A[9]*W[0];
          W[17] = A[6]*W[7]-A[7]*W[3]+A[10]*W[0];
          W[18] = A[6]*W[8]-A[7]*W[4]+A[11]*W[0];
          W[19] = A[6]*W[9]-A[8]*W[2]+A[9]*W[1];
          W[20] = A[6]*W[10]-A[8]*W[3]+A[10]*W[1];
          W[21] = A[6]*W[11]-A[8]*W[4]+A[11]*W[1];
          W[22] = A[6]*W[12]-A[9]*W[3]+A[10]*W[2];
          W[23] = A[6]*W[13]-A[9]*W[4]+A[11]*W[2];
          W[24] = A[6]*W[14]-A[10]*W[4]+A[11]*W[3];
          W[25] = A[7]*W[9]-A[8]*W[6]+A[9]*W[5];
          W[26] = A[7]*W[10]-A[8]*W[7]+A[10]*W[5];
          W[27] = A[7]*W[11]-A[8]*W[8]+A[11]*W[5];
          W[28] = A[7]*W[12]-A[9]*W[7]+A[10]*W[6];
          W[29] = A[7]*W[13]-A[9]*W[8]+A[11]*W[6];
          W[30] = A[7]*W[14]-A[10]*W[8]+A[11]*W[7];
          W[31] = A[8]*W[12]-A[9]*W[10]+A[10]*W[9];
          W[32] = A[8]*W[13]-A[9]*W[11]+A[11]*W[9];
          W[33] = A[8]*W[14]-A[10]*W[11]+A[11]*W[10];
          W[34] = A[9]*W[14]-A[10]*W[13]+A[11]*W[12];
          W[ 0] = A[0]*W[25]-A[1]*W[19]+A[2]*W[16]-A[3]*W[15];
          W[ 1] = A[0]*W[26]-A[1]*W[20]+A[2]*W[17]-A[4]*W[15];
          W[ 2] = A[0]*W[27]-A[1]*W[21]+A[2]*W[18]-A[5]*W[15];
          W[ 3] = A[0]*W[28]-A[1]*W[22]+A[3]*W[17]-A[4]*W[16];
          W[ 4] = A[0]*W[29]-A[1]*W[23]+A[3]*W[18]-A[5]*W[16];
          W[ 5] = A[0]*W[30]-A[1]*W[24]+A[4]*W[18]-A[5]*W[17];
          W[ 6] = A[0]*W[31]-A[2]*W[22]+A[3]*W[20]-A[4]*W[19];
          W[ 7] = A[0]*W[32]-A[2]*W[23]+A[3]*W[21]-A[5]*W[19];
          W[ 8] = A[0]*W[33]-A[2]*W[24]+A[4]*W[21]-A[5]*W[20];
          W[ 9] = A[0]*W[34]-A[3]*W[24]+A[4]*W[23]-A[5]*W[22];
          W[10] = A[1]*W[31]-A[2]*W[28]+A[3]*W[26]-A[4]*W[25];
          W[11] = A[1]*W[32]-A[2]*W[29]+A[3]*W[27]-A[5]*W[25];
          W[12] = A[1]*W[33]-A[2]*W[30]+A[4]*W[27]-A[5]*W[26];
          W[13] = A[1]*W[34]-A[3]*W[30]+A[4]*W[29]-A[5]*W[28];
          W[14] = A[2]*W[34]-A[3]*W[33]+A[4]*W[32]-A[5]*W[31];
          B[ 4] = d*( A[31]*W[14]-A[32]*W[13]+A[33]*W[12]-A[34]*W[11]+A[35]*W[10]);
          B[10] = d*(-A[30]*W[14]+A[32]*W[9]-A[33]*W[8]+A[34]*W[7]-A[35]*W[6]);
          B[16] = d*( A[30]*W[13]-A[31]*W[9]+A[33]*W[5]-A[34]*W[4]+A[35]*W[3]);
          B[22] = d*(-A[30]*W[12]+A[31]*W[8]-A[32]*W[5]+A[34]*W[2]-A[35]*W[1]);
          B[28] = d*( A[30]*W[11]-A[31]*W[7]+A[32]*W[4]-A[33]*W[2]+A[35]*W[0]);
          B[34] = d*(-A[30]*W[10]+A[31]*W[6]-A[32]*W[3]+A[33]*W[1]-A[34]*W[0]);
          B[ 5] = d*(-A[25]*W[14]+A[26]*W[13]-A[27]*W[12]+A[28]*W[11]-A[29]*W[10]);
          B[11] = d*( A[24]*W[14]-A[26]*W[9]+A[27]*W[8]-A[28]*W[7]+A[29]*W[6]);
          B[17] = d*(-A[24]*W[13]+A[25]*W[9]-A[27]*W[5]+A[28]*W[4]-A[29]*W[3]);
          B[23] = d*( A[24]*W[12]-A[25]*W[8]+A[26]*W[5]-A[28]*W[2]+A[29]*W[1]);
          B[29] = d*(-A[24]*W[11]+A[25]*W[7]-A[26]*W[4]+A[27]*W[2]-A[29]*W[0]);
          B[35] = d*( A[24]*W[10]-A[25]*W[6]+A[26]*W[3]-A[27]*W[1]+A[28]*W[0]);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Directly inverts a matrix.
     *
     * This function inverts a (small) regular matrix whose size is constant at compile time. Up to a size of 6x6
     * matrices this function uses a direct inversion approach based on <em>Cramer's rule</em> rather than performing
     * a pivoted LU-factorisation, which results in a both fast and stable calculation.
     * \see http://en.wikipedia.org/wiki/Cramer's_rule
     * \see http://mathworld.wolfram.com/CramersRule.html
     *
     * \tparam n_
     * The number of rows/columns of the matrix to be inverted. Has to be specified explicitly.
     *
     * \warning
     * Compilation of this function will abort with a static assertion failure if \p n_ is greater than 8. This is
     * not a bug, it's a feature to prevent stack overflows. If you feel like inverting bigger matrices, use the
     * #mat_factorise and #mat_solve_vec / #mat_solve_mat functions instead.
     *
     * \param[out] b
     * The n-by-n matrix \e B which recieves the inverse of \e A.
     *
     * \param[in] a
     * The n-by-n matrix \e A which is to be inverted.
     *
     * \warning
     * This function does not check the input matrix \p a for regularity, so use this function only when you are sure
     * that the input matrix is regular. In addition to that, it is recommended to call this function within a
     * \c try block to catch possibly thrown floating point exceptions.
     *
     * \author Peter Zajac
     */
    template<
      int n_,
      typename Type_>
    inline void mat_invert(
      Type_ b[],
      const Type_ a[])
    {
      CONTEXT("LinAlg::mat_invert()");
      // ensure that the matrix is "small"
      static_assert(n_ >= 1, "Invalid matrix size");
      static_assert(n_ <= 8, "Matrix size too big; use mat_factorise instead");
      Intern::MatInvert<false, n_, n_>::apply(b, a);
    }

    /** \copydoc mat_invert(Type_ b[], const Type_ a[]) */
    template<
      int m_,
      int n_,
      typename Type_>
    inline void mat_invert(
      Type_ (&b)[m_][n_],
      const Type_ (&a)[m_][n_])
    {
      CONTEXT("LinAlg::mat_invert()");
      // ensure that the matrix is "small"
      static_assert(m_ >= 1, "Invalid matrix size");
      static_assert(n_ >= 1, "Invalid matrix size");
      static_assert(n_ <= 8, "Matrix size too big; use mat_factorise instead");
      // Note: We do not check whether m_ == n_ at compile-time; this is handled by
      //       an exception within the following apply()-function at runtime.
      Intern::MatInvert<false, m_, n_>::apply(&b[0][0], &a[0][0]);
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

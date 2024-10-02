// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_TINY_ALGEBRA_HPP
#define KERNEL_UTIL_TINY_ALGEBRA_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#ifndef __CUDA_ARCH__
#include <kernel/util/math.hpp>
#endif

// includes, system
#ifndef __CUDACC__
#include <initializer_list>
#else
#include <kernel/util/cuda_math.cuh>
#endif

namespace FEAT
{
  /**
   * \brief Tiny namespace
   *
   * This namespace encapsulates Vector, Matrix and third-order Tensor classes, whose dimensions are given at
   * compile-time. These classes are commonly used by various other kernel components, most notably by the assembly.
   */
  namespace Tiny
  {
    /**
     * \brief Tiny Vector class template
     *
     * This class template implements a vector whose value type and size is given at compile-time. The value type can
     * be a primitive type, or some other object like (again) a Vector.
     *
     * \tparam T_
     * The value-type that the vector shall contain.
     *
     * \tparam n_
     * The length of the vector. Must be > 0.
     *
     * \tparam s_
     * The stride of the vector. Must be >= \p n_.
     *
     * \author Peter Zajac
     */
    template<
      typename T_,
      int n_,
      int s_ = n_>
    class Vector DOXY({});

    /**
     * \brief Tiny Matrix class template
     *
     * This class template implements a matrix whose value type and size is given at compile-time. The value type can
     * be a primitive type, or some other object like a Vector.
     *
     * \tparam T_
     * The datatype that the vector shall contain.
     *
     * \tparam m_, n_
     * The number rows and columns of the matrix. Must be > 0.
     *
     * \tparam sm_
     * The row stride of the matrix. Must be >= \p m_.
     *
     * \tparam sn_
     * The column stride of the matrix. Must be >= \p n_.
     *
     * \author Peter Zajac
     */
    template<
      typename T_,
      int m_,
      int n_,
      int sm_ = m_,
      int sn_ = n_>
    class Matrix DOXY({});

    /**
     * \brief Tiny Tensor3 class template
     *
     * This class template implements a 3rd-order tensor whose datatype and sizes are given at compile-time.
     * This class template implements a 3rd order tensor whose value type and size is given at compile-time. The value
     * type can be a primitive type, or some other object like a Vector.
     *
     * Technically, the Tensor3 class template realizes a l-tuple of m-by-n matrices.
     *
     * \tparam T_
     * The datatype that the vector shall contain.
     *
     * \tparam l_, m_, n_
     * The number tubes, rows and columns of the tensor. Must be > 0.
     *
     * \tparam sl_,
     * The tube stride of the tensor. Must be >= \p l_.
     *
     * \tparam sm_
     * The row stride of the tensor. Must be >= \p m_.
     *
     * \tparam sn_
     * The column stride of the tensor. Must be >= \p n_.
     *
     * \author Peter Zajac
     */
    template<
      typename T_,
      int l_,
      int m_,
      int n_,
      int sl_ = l_,
      int sm_ = m_,
      int sn_ = n_>
    class Tensor3 DOXY({});

    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Helper class that extracts the basic data type from another type
       *
       * \tparam T_
       * The type
       *
       * This is the end of the recursion where the type is something other than a FEAT::Tiny::... object.
       **/
      template<typename T_>
      struct DataTypeExtractor
      {
        /// The basic data type
        typedef T_ MyDataType;
        /// Level of recursion, see below
        static constexpr int level = 0;
      };

      /**
       * \brief Class template that recursively typedefs itself
       *
       * \tparam T_
       * The datatype that the vector shall contain.
       *
       * \tparam n_
       * The length of the vector. Must be > 0.
       *
       * \tparam s_
       * The stride of the vector. Must be >= \p n_.
       *
       * This goes on as long as the underlying type is still some kind of Vector.
       **/
      template<typename T_, int n_, int s_>
      struct DataTypeExtractor<Vector<T_, n_, s_>>
      {
        /// Recursive typedef
        typedef typename DataTypeExtractor<T_>::MyDataType MyDataType;
        /// The number of times the ValueType is some kind of Vector. A Vector is a tensor of order 1+level.
        static constexpr int level = DataTypeExtractor<T_>::level+1;
      };

      // Same for Matrix
      template<typename T_, int m_, int n_, int sm_, int sn_>
      struct DataTypeExtractor<Matrix<T_, m_, n_, sm_, sn_>>
      {
        /// Recursive typedef
        typedef typename DataTypeExtractor<T_>::MyDataType MyDataType;
        /// The number of times the ValueType is some kind of Vector. A Vector is a tensor of order 1+level.
        static constexpr int level = DataTypeExtractor<T_>::level+1;
      };

      // Same for Tensor3
      template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
      struct DataTypeExtractor<Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>>
      {
        /// Recursive typedef
        typedef typename DataTypeExtractor<T_>::MyDataType MyDataType;
        /// The number of times the ValueType is some kind of Vector. A Vector is a tensor of order 1+level.
        static constexpr int level = DataTypeExtractor<T_>::level+1;
      };

      // forward declarations of helper classes
      template<int m_, int n_, int sna_>
      struct DetHelper;

      template<int m_, int n_, int sna_>
      struct VolHelper;

      template<int m_, int n_, int sna_, int snb_>
      struct InverseHelper;

      template<int m_, int n_, int sna_, int snb_>
      struct CofactorHelper;
    } // namespace Intern
    /// \endcond

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    // Tiny Vector implementation
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    template<typename T_, int n_, int s_>
    class Vector
    {
      static_assert(n_ > 0, "invalid vector length");
      static_assert(s_ >= n_, "invalid vector stride");

    public:
      /// the length of the vector
      static constexpr int n = n_;
      /// the stride of the vector
      static constexpr int s = s_;

      /// the value type of the vector
      typedef T_ ValueType;
      /// The basic data type buried in the lowest level of the vector
      typedef typename Intern::DataTypeExtractor<ValueType>::MyDataType DataType;

      /// actual vector data
      T_ v[s_];

      /// default constructor
      CUDA_HOST_DEVICE Vector()
      {
      }

      /// \brief value-assignment constructor
      CUDA_HOST_DEVICE explicit Vector(DataType value)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = value;
        }
      }

      /// copy constructor
      template<int sx_>
      CUDA_HOST_DEVICE Vector(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = x.v[i];
        }
      }

      /**
       * \brief Initializer list constructor
       *
       * This operator allows to assign values to a new vector in a simple manner:
       * \code{.cpp}
         Tiny::Vector<double, 3> v1{0.1, 2.1, 7.3};
         Tiny::Vector<double, 3> v2 = {0.1, 2.1, 7.3};
       * \endcode
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE explicit Vector(const std::initializer_list<Tx_>& x)
      {
        XASSERTM(std::size_t(n_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < n_; ++i, ++it)
          v[i] = T_(*it);
      }

      /// value-assignment operator
      CUDA_HOST_DEVICE Vector& operator=(DataType value)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = value;
        }
        return *this;
      }

      /// copy-assignment operator
      template<int sx_>
      CUDA_HOST_DEVICE Vector& operator=(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = x.v[i];
        }
        return *this;
      }

      /**
       * \brief Initializer list assignment operator
       *
       * This operator allows to assign values to the vector in a simple manner:
       * \code{.cpp}
         Tiny::Vector<double, 3> v;
         v = {0.1, 2.1, 7.3};
       * \endcode
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       *
       * \returns *this
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE Vector& operator=(const std::initializer_list<Tx_>& x)
      {
        XASSERTM(std::size_t(n_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < n_; ++i, ++it)
          v[i] = T_(*it);
        return *this;
      }

      /**
       * \brief Access operator.
       *
       * \param[in] i
       * The index of the vector component to be returned.
       *
       * \returns A (const) reference to the <c>i</c>-th entry of the vector.
       */
      CUDA_HOST_DEVICE T_& operator()(int i)
      {
        ASSERTM((i >= 0) && (i < n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator()() */
      CUDA_HOST_DEVICE const T_& operator()(int i) const
      {
        ASSERTM((i >= 0) && (i < n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator()() */
      CUDA_HOST_DEVICE T_& operator[](int i)
      {
        ASSERTM((i >= 0) && (i < n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator[]() */
      CUDA_HOST_DEVICE const T_& operator[](int i) const
      {
        ASSERTM((i >= 0) && (i < n_), "index i out-of-bounds");
        return v[i];
      }

      /**
       * \brief Size-cast function.
       *
       * This function casts this vector's reference to another size.
       *
       * \tparam nn_
       * The length of the casted vector. Must be 0 < \p nn_ <= \p s_.
       *
       * \returns
       * A casted (const) reference of \p *this.
       */
      template<int nn_>
      CUDA_HOST_DEVICE Vector<T_, nn_, s_>& size_cast()
      {
        static_assert((nn_ > 0) && (nn_ <= s_), "invalid cast length");
        return reinterpret_cast<Vector<T_, nn_, s_>&>(*this);
      }

      /** \copydoc size_cast() */
      template<int nn_>
      CUDA_HOST_DEVICE const Vector<T_, nn_, s_>& size_cast() const
      {
        static_assert((nn_ > 0) && (nn_ <= s_), "invalid cast length");
        return reinterpret_cast<const Vector<T_, nn_, s_>&>(*this);
      }

      /// scalar-multiply operator
      CUDA_HOST_DEVICE Vector& operator*=(DataType alpha)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] *= alpha;
        }
        return *this;
      }

      /// element-wise-multiply operator
      template <int sx_>
      CUDA_HOST_DEVICE Vector& operator*=(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] *= x.v[i];
        }
        return *this;
      }

      /// vector-add operator
      template<int sx_>
      CUDA_HOST_DEVICE Vector& operator+=(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] += x.v[i];
        }
        return *this;
      }

      /// vector-subtract operator
      template<int sx_>
      CUDA_HOST_DEVICE Vector& operator-=(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] -= x.v[i];
        }
        return *this;
      }

      /**
       * \brief Format the vector.
       *
       * \param[in] alpha
       * The value that the vector is to be set to.
       */
      CUDA_HOST_DEVICE void format(DataType alpha = DataType(0))
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = alpha;
        }
      }

      /**
       * \brief Normalizes this vector.
       *
       * \returns \c *this
       */
      CUDA_HOST_DEVICE Vector& normalize()
      {
      #ifndef __CUDACC__
        const DataType norm2(this->norm_euclid());
        ASSERTM(norm2 > Math::eps<DataType>(), "Trying to normalize a null vector!");
        return ((*this) *= (DataType(1)/norm2));
      #else
        const DataType norm2_sqr(this->norm_euclid_sqr());
        ASSERTM(norm2_sqr > CudaMath::cuda_get_eps<DataType>(), "Trying to normalize a null vector!");
        return ((*this) *= CudaMath::cuda_rsqrt(norm2_sqr));
      #endif
      }

      /**
       * \brief Negates the vector, i.e. effectively multiplies all components by -1
       *
       * \returns \c *this
       */
      CUDA_HOST_DEVICE Vector& negate()
      {
        for(int i(0); i < n_; ++i)
          v[i] = -v[i];
        return *this;
      }

      /**
       * \brief Adds another scaled vector onto this vector.
       *
       * \param[in] x
       * The vector to be added onto this vector.
       *
       * \param[in] alpha
       * The scaling parameter for the axpy.
       *
       * \returns \c *this
       */
      template<int snx_>
      CUDA_HOST_DEVICE Vector& axpy(DataType alpha, const Vector<T_, n_, snx_>& x)
      {
        for(int i(0); i < n_; ++i)
          v[i] += alpha * x.v[i];
        return *this;
      }

      /**
       * \brief Sets this vector to the convex combination of two other vectors.
       *
       * Let \e y denote \c this vector, then this function computes:
       * \f[ y \leftarrow (1-\alpha)\cdot a + \alpha\cdot b \f]
       *
       * \param[in] alpha
       * The interpolation parameter for the convex combination. Should be 0 <= alpha <= 1
       *
       * \param[in] a
       * The first vector for the convex combination.
       *
       * \param[in] b
       * The second vector for the convex combination.
       *
       * \returns \c *this
       */
      template<int sna_, int snb_>
      CUDA_HOST_DEVICE Vector& set_convex(DataType alpha, const Vector<T_, n_, sna_>& a, const Vector<T_, n_, snb_>& b)
      {
        for(int i(0); i < n_; ++i)
          v[i] = (T_(1) - alpha) * a.v[i] + alpha * b.v[i];
        return *this;
      }

      /**
       * \brief Sets this vector to the result of a matrix-vector product.
       *
       * Let \e y denote \c this vector, and let \e A denote the input matrix and \e x in the input vector, then
       * this function computes:
       * \f[ y \leftarrow A\cdot x \f]
       *
       * \param[in] a
       * The matrix for the product.
       *
       * \param[in] x
       * The (right) multiplicand vector for the product.
       *
       * \returns \p *this
       */
      template<int m_, int sma_, int sna_, int sx_>
      CUDA_HOST_DEVICE Vector& set_mat_vec_mult(const Matrix<T_, n_, m_, sma_, sna_>& a, const Vector<T_, m_, sx_>& x)
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&x, "result vector and multiplicand vector 'x' must be different objects");

        for(int i(0); i < n_; ++i)
        {
          v[i] = T_(0);
          for(int j(0); j < m_; ++j)
          {
            v[i] += a.v[i][j] * x.v[j];
          }
        }
        return *this;
      }

      /**
       * \brief Sets this vector to the result of a vector-matrix product.
       *
       * Let \e y denote \c this vector, and let \e x in the input vector and \e A denote the input matrix, then
       * this function computes:
       * \f[ y^\top \leftarrow x^\top \cdot A \Longleftrightarrow y \leftarrow A^\top\cdot x \f]
       *
       * \param[in] x
       * The (left) multiplicand vector for the product.
       *
       * \param[in] a
       * The matrix for the product.
       *
       * \returns \p *this
       */
      template<int m_, int sma_, int sna_, int sx_>
      CUDA_HOST_DEVICE Vector& set_vec_mat_mult(const Vector<T_, m_, sx_>& x, const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&x, "result vector and multiplicand vector 'x' must be different objects");

        for(int j(0); j < n_; ++j)
        {
          v[j] = T_(0);
          for(int i(0); i < m_; ++i)
          {
            v[j] += a.v[i][j] * x.v[i];
          }
        }
        return *this;
      }

      /**
       * \brief Adds the result of a matrix-vector product onto this vector.
       *
       * Let \e y denote \c this vector, and let \e A denote the input matrix and \e x in the input vector, then
       * this function computes:
       * \f[ y \leftarrow y + \alpha\cdot A\cdot x \f]
       *
       * \param[in] a
       * The matrix for the product.
       *
       * \param[in] x
       * The (right) multiplicand vector for the product.
       *
       * \param[in] alpha
       * The scaling parameter for the product.
       *
       * \returns \p *this
       */
      template<int m_, int sma_, int sna_, int sx_>
      CUDA_HOST_DEVICE Vector& add_mat_vec_mult(const Matrix<T_, n_, m_, sma_, sna_>& a, const Vector<T_, m_, sx_>& x, DataType alpha = DataType(1))
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&x, "result vector and multiplicand vector 'x' must be different objects");

        for(int i(0); i < n_; ++i)
        {
          for(int j(0); j < m_; ++j)
          {
            v[i] += alpha * a.v[i][j] * x.v[j];
          }
        }
        return *this;
      }

      /**
       * \brief Adds the result of a vector-matrix product onto this vector.
       *
       * Let \e y denote \c this vector, and let \e x in the input vector and \e A denote the input matrix, then
       * this function computes:
       * \f[ y^\top \leftarrow y^\top + \alpha\cdot x^\top \cdot A \Longleftrightarrow y \leftarrow y + \alpha\cdot A^\top\cdot x \f]
       *
       * \param[in] x
       * The (left) multiplicand vector for the product.
       *
       * \param[in] a
       * The matrix for the product.
       *
       * \param[in] alpha
       * The scaling parameter for the product.
       *
       * \returns \p *this
       */
      template<int m_, int sma_, int sna_, int sx_>
      CUDA_HOST_DEVICE Vector& add_vec_mat_mult(const Vector<T_, m_, sx_>& x, const Matrix<T_, m_, n_, sma_, sna_>& a, DataType alpha = DataType(1))
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&x, "result vector and multiplicand vector 'x' must be different objects");

        for(int j(0); j < n_; ++j)
        {
          for(int i(0); i < m_; ++i)
          {
            v[j] += alpha * a.v[i][j] * x.v[i];
          }
        }
        return *this;
      }

      /**
       * \brief Computes the squared euclid norm of the vector.
       *
       * \returns
       * The squared euclid norm of the vector.
       */
      CUDA_HOST_DEVICE DataType norm_euclid_sqr() const
      {
        DataType r(DataType(0));
        for(int i(0); i < n_; ++i)
          #ifndef __CUDACC__
          r += Math::sqr(v[i]);
          #else
          r += CudaMath::cuda_sqr(v[i]);
          #endif
        return r;
      }

      /**
       * \brief Computes the euclid norm of the vector.
       *
       * \returns
       * The euclid norm of the vector.
       */
      CUDA_HOST_DEVICE DataType norm_euclid() const
      {
      #ifndef __CUDACC__
        return Math::sqrt(norm_euclid_sqr());
      #else
        return CudaMath::cuda_sqrt(norm_euclid_sqr());
      #endif
      }

      /**
       * \brief Computes the l1-norm of the vector.
       *
       * \returns
       * The l1-norm of the vector.
       */
      CUDA_HOST_DEVICE DataType norm_l1() const
      {
        DataType r(DataType(0));
        for(int i(0); i < n_; ++i)
          #ifndef __CUDACC__
          r += Math::abs(v[i]);
          #else
          r += CudaMath::cuda_abs(v[i]);
          #endif
        return r;
      }

      /**
       * \brief Computes the max-norm of the vector.
       *
       * \returns
       * The max-norm of the vector.
       */
      CUDA_HOST_DEVICE DataType norm_max() const
      {
        DataType r(DataType(0));
        for(int i(0); i < n_; ++i)
        #ifndef __CUDACC__
          r = Math::max(r, Math::abs(v[i]));
        #else
          r = CudaMath::cuda_max(r, CudaMath::cuda_abs(v[i]));
        #endif
        return r;
      }

      /**
       * \brief Returns a null-vector.
       */
      CUDA_HOST_DEVICE static Vector null()
      {
        return Vector(DataType(0));
      }

      /**
       * \brief Tiny::Vector streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The vector to be streamed.
       */
      CUDA_HOST friend std::ostream & operator<< (std::ostream & lhs, const Vector & b)
      {
        lhs << "[";
        for (int i(0) ; i < b.n ; ++i)
        {
          lhs << "  " << stringify(b(i));
        }
        lhs << "]";

        return lhs;
      }
    }; // class Vector

    /*
    template<typename T_, int sx_, int sa_>
    inline void cross(Vector<T_, 2, sx_>& x, const Vector<T_, 2, sa_>& a)
    {
      x.v[0] =  a.v[1];
      x.v[1] = -a.v[0];
    }

    template<typename T_, int sx_, int sa_, int sb_>
    inline void cross(Vector<T_, 3, sx_>& x, const Vector<T_, 3, sa_>& a, const Vector<T_, 3, sb_>& b)
    {
      x.v[0] = a.v[1]*b.v[2] - a.v[2]*b.v[1];
      x.v[1] = a.v[2]*b.v[0] - a.v[0]*b.v[2];
      x.v[2] = a.v[0]*b.v[1] - a.v[1]*b.v[0];
    }*/

    /// scalar-left-multiplication operator
    template<typename T_, int n_, int s_>
    CUDA_HOST_DEVICE inline Vector<T_, n_> operator*(typename Vector<T_, n_>::DataType alpha, const Vector<T_, n_, s_>& x)
    {
      return Vector<T_, n_>(x) *= alpha;
    }

    /// scalar-right-multiplication operator
    template<typename T_, int n_, int s_>
    CUDA_HOST_DEVICE inline Vector<T_, n_> operator*(const Vector<T_, n_, s_>& x, typename Vector<T_, n_>::DataType alpha)
    {
      return Vector<T_, n_>(x) *= alpha;
    }

    /// vector element-wise-product operator
    template<typename T_, int n_, int sa_, int sb_>
    CUDA_HOST_DEVICE inline Vector<T_, n_> component_product(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      return Vector<T_, n_>(a) *= b;
    }

    /// vector addition operator
    template<typename T_, int n_, int sa_, int sb_>
    CUDA_HOST_DEVICE inline Vector<T_, n_> operator+(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      return Vector<T_, n_>(a) += b;
    }

    /// vector subtraction operator
    template<typename T_, int n_, int sa_, int sb_>
    CUDA_HOST_DEVICE inline Vector<T_, n_> operator-(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      return Vector<T_, n_>(a) -= b;
    }

      /**
       * \brief Calculates the counter-clockwise opening angle between two 2D vectors
       *
       * \param[in] x Vector x
       * \param[in] y Vector y
       */
      template<typename T_>
      CUDA_HOST inline T_ calculate_opening_angle(const Vector<T_,2>& x, const Vector<T_, 2>& y)
      {
        #ifdef __CUDACC__
        return T_(0);
        #else
        return Math::calc_opening_angle(x[0], x[1], y[0], y[1]);
        #endif
      }


    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    // Tiny Matrix implementation
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    template<typename T_, int m_, int n_, int sm_, int sn_>
    class Matrix
    {
      static_assert(m_ > 0, "invalid row count");
      static_assert(n_ > 0, "invalid column count");
      static_assert(sm_ >= m_, "invalid row stride");
      static_assert(sn_ >= n_, "invalid column stride");

    public:
      /// the row count of the matrix
      static constexpr int m = m_;
      /// the column count of the matrix
      static constexpr int n = n_;
      /// the row stride of the matrix
      static constexpr int sm = sm_;
      /// the column stride of the matrix
      static constexpr int sn = sn_;

      /// the data type of the matrix
      typedef T_ ValueType;
      /// The basic data type buried in the lowest level of the vector
      typedef typename Intern::DataTypeExtractor<ValueType>::MyDataType DataType;

      /// the type of a single matrix row
      typedef Vector<T_, n_, sn_> RowType;
      /// actual matrix data; that's an array of vectors
      RowType v[sm_];

      /// default constructor
      CUDA_HOST_DEVICE Matrix()
      {
      }

      /// value-assignment constructor
      CUDA_HOST_DEVICE explicit Matrix(DataType value)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] = value;
        }
      }

      /// copy constructor
      template<typename T2_, int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix(const Matrix<T2_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] = T_(a.v[i][j]);
          }
        }
      }

      /**
       * \brief Initializer list of Tiny::Vector constructor
       *
       * \code{.cpp}
       * Tiny::Vector<double, 2> v{...}, w{...};
       * Tiny::Matrix<double, 2, 2> m{v, w};
       * \endcode
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE explicit Matrix(const std::initializer_list<Tx_>& x)
      {
        XASSERTM(std::size_t(m_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < m_; ++i, ++it)
          v[i] = *it;
      }

      /**
       * \brief Initializer list constructor
       *
       * \note This overload seems to be necessary, because the example code below
       * does not compile otherwise for some reason that I do not understand...
       *
       * \code{.cpp}
       * Tiny::Matrix<double, 2, 2> v{{0.0, 1.0}, {2.0, 3.0}};
       * \endcode
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE explicit Matrix(const std::initializer_list<std::initializer_list<Tx_>>& x)
      {
        XASSERTM(std::size_t(m_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < m_; ++i, ++it)
          v[i] = *it;
      }

      /// value-assignment operator
      CUDA_HOST_DEVICE Matrix& operator=(DataType value)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] = value;
        }
        return *this;
      }

      /// assignment operator
      template<int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix& operator=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] = a.v[i];
        }
        return *this;
      }

      /**
       * \brief Initializer list assignment operator
       *
       * \code{.cpp}
       * Tiny::Vector<double, 2> v{1.0, 2.0}, w{3.0, 4.0};
       * Tiny::Matrix<double, 2, 2> m;
         m = {v, w};
       * \endcode
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE Matrix& operator=(const std::initializer_list<Tx_>& x)
      {
        XASSERTM(std::size_t(m_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < m_; ++i, ++it)
          v[i] = *it;
        return *this;
      }

      /**
       * \brief Initializer list assignment operator
       *
       * \code{.cpp}
       * Tiny::Matrix<double, 2, 2> m;
         m = {{0.0, 1.0}, {2.0, 3.0}};
       * \endcode
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE Matrix& operator=(const std::initializer_list<std::initializer_list<Tx_>>& x)
      {
        XASSERTM(std::size_t(m_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < m_; ++i, ++it)
          v[i] = *it;
        return *this;
      }

      /**
       * \brief Access operator.
       *
       * \param[in] i,j
       * The indices of the matrix entry that is to be returned.
       *
       * \returns
       * A (const) reference to the matrix entry at position (i,j).
       */
      CUDA_HOST_DEVICE T_& operator()(int i, int j)
      {
        ASSERTM( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERTM( (j >= 0) && (j < n_), "index j out-of-bounds");
        return v[i][j];
      }

      /** \copydoc operator()() */
      CUDA_HOST_DEVICE const T_& operator()(int i, int j) const
      {
        ASSERTM( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERTM( (j >= 0) && (j < n_), "index j out-of-bounds");
        return v[i][j];
      }

      /**
       * \brief Row-Access operator.
       *
       * \param[in] i
       * The index of the row that is to be returned.
       *
       * \returns
       * A (const) reference to the <c>i</c>-th row of the matrix.S
       */
      CUDA_HOST_DEVICE RowType& operator[](int i)
      {
        ASSERTM( (i >= 0) && (i <m_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator[]() */
      CUDA_HOST_DEVICE const RowType& operator[](int i) const
      {
        ASSERTM( (i >= 0) && (i <m_), "index i out-of-bounds");
        return v[i];
      }

      /**
       * \brief Size-cast function.
       *
       * This function casts this matrix's reference to another size.
       *
       * \tparam mm_, nn_
       * The dimensions of the casted matrix. Must be 0 < \p mm_ <= \p sm_ and 0 < \p nn_ <= \p sn_.
       *
       * \returns
       * A casted (const) reference of \p *this.
       */
      template<int mm_, int nn_>
      CUDA_HOST_DEVICE Matrix<T_, mm_, nn_, sm_, sn_>& size_cast()
      {
        static_assert((mm_ > 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ > 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<Matrix<T_, mm_, nn_, sm_, sn_>&>(*this);
      }

      /** \copydoc size_cast() */
      template<int mm_, int nn_>
      CUDA_HOST_DEVICE const Matrix<T_, mm_, nn_, sm_, sn_>& size_cast() const
      {
        static_assert((mm_ > 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ > 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<const Matrix<T_, mm_, nn_, sm_, sn_>&>(*this);
      }

      /// scalar-right-multiply-by operator
      CUDA_HOST_DEVICE Matrix& operator*=(DataType alpha)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] *= alpha;
        }
        return *this;
      }

      /// matrix component-wise addition operator
      template<int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix& operator+=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] += a.v[i];
        }
        return *this;
      }

      /// matrix component-wise subtraction operator
      template<int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix& operator-=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] -= a.v[i];
        }
        return *this;
      }

      /**
       * \brief Formats the matrix.
       *
       * \param[in] alpha
       * The value that the matrix is to be set to.
       */
      CUDA_HOST_DEVICE void format(DataType alpha = DataType(0))
      {
        for(int i(0); i < m_; ++i)
        {
          v[i].format(alpha);
        }
      }

      /**
       * \brief Returns the Hessian norm square.
       *
       * This function computes and returns
       * \f[ \sum_{i=0}^{m-1}\sum_{j=0}^{n-1} K_{ij}\cdot (A_{ij})^2 \f]
       * where K_ij is 1 for i=j and 1/2 otherwise.
       *
       * \note This function is used for the computation of the H^2 semi-norm of a function.
       *
       * \returns The Hessian norm square of the matrix.
       */
      CUDA_HOST_DEVICE DataType norm_hessian_sqr() const
      {
        DataType r(0);
        for(int i(0); i < m_; ++i)
        {
          #ifndef __CUDACC__
          r += Math::sqr(v[i][i]);
          #else
          r += CudaMath::cuda_sqr(v[i][i]);
          #endif
          for(int j(0); j < n_; ++j)
          {
            #ifndef __CUDACC__
            r += Math::sqr(v[i][j]);
            #else
            r += CudaMath::cuda_sqr(v[i][j]);
            #endif
          }
        }
        return r / DataType(2);
      }

      /**
       * \brief Returns the Frobenius norm of the matrix.
       *
       * This function computes and returns
       * \f[ \Big(\sum_{i=0}^{m-1}\sum_{j=0}^{n-1} (A_{ij})^2\Big)^{\frac{1}{2}} \f]
       *
       * \returns The Frobenius norm of the matrix.
       */
      CUDA_HOST_DEVICE DataType norm_frobenius() const
      {
        #ifndef __CUDACC__
        return Math::sqrt(norm_frobenius_sqr());
        #else
        return CudaMath::cuda_sqrt(norm_frobenius_sqr());
        #endif
      }

      /**
       * \brief Returns the Frobenius norm squared of the matrix.
       *
       * This function computes and returns
       * \f[ \sum_{i=0}^{m-1}\sum_{j=0}^{n-1} (A_{ij})^2\f]
       *
       * \returns The Frobenius norm of the matrix.
       */
      CUDA_HOST_DEVICE DataType norm_frobenius_sqr() const
      {
        DataType r(0);
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            #ifndef __CUDACC__
            r += Math::sqr(v[i][j]);
            #else
            r += CudaMath::cuda_sqr(v[i][j]);
            #endif
          }
        }
        return r;
      }

      /**
       * \brief Returns the Frobenius norm of the difference of this matrix and the identity matrix.
       *
       * \returns The Frobenius norm of the difference of this matrix and the identity matrix.
       */
      CUDA_HOST_DEVICE DataType norm_sub_id_frobenius() const
      {
        DataType r(0);
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            #ifndef __CUDACC__
            r += Math::sqr(v[i][j] - DataType(i == j ? 1 : 0));
            #else
            r += CudaMath::cuda_sqr(v[i][j] - DataType(i == j ? 1 : 0));
            #endif
          }
        }
        #ifndef __CUDACC__
        return Math::sqrt(r);
        #else
        return CudaMath::cuda_sqrt(r);
        #endif
      }

      /**
       * \brief Returns the trace of the matrix.
       *
       * This function computes and returns
       * \f[ \sum_{i=0}^{\min(m,n)} A_{ii}\f]
       * i.e. the sum of all main diagonal elements.
       *
       * \returns The trace of the matrix.
       */
      CUDA_HOST_DEVICE DataType trace() const
      {
        #ifndef __CUDACC__
        int k = Math::min(m_, n_);
        #else
        int k = CudaMath::cuda_min(m_, n_);
        #endif
        DataType r(0);
        for(int i(0); i < k; ++i)
        {
          r += v[i][i];
        }
        return r;
      }

      /**
       * \brief Returns the determinant of the matrix.
       *
       * \warning This function only works for \p m_ = \p n_ and will intentionally fail to compile in any other case.
       *
       * \returns The determinant of the matrix.
       */
      CUDA_HOST_DEVICE DataType det() const
      {
        return Intern::DetHelper<m_, n_, sn_>::compute(&this->v[0].v[0]);
      }

      /**
       * \brief Returns the volume of the matrix.
       *
       * The volume of a matrix \e A is defined as:
       * \f[ \textnormal{vol}(A) := \sqrt{\textnormal{det}(A^\top\cdot A)} \f]
       *
       * \note For \p m_ = \p n_, the volume of the matrix equals the absolute of its determinant.
       * \note For \p m_ < \p n_, the volume of the matrix is always zero.
       *
       * \returns The volume of the matrix.
       */
      CUDA_HOST_DEVICE DataType vol() const
      {
        return Intern::VolHelper<m_, n_, sn_>::compute(&this->v[0].v[0]);
      }

      /**
       * \brief Sets this matrix to the inverse of another matrix.
       *
       * \warning This function only works for \p m_ = \p n_ and will intentionally fail to compile in any other case.
       *
       * \param[in] a
       * The matrix whose inverse is to be stored in this matrix.
       *
       * \returns \p *this
       */
      template<int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix& set_inverse(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and input matrix 'a' must be different objects");
        Intern::InverseHelper<m_, n_, sn_, sna_>::compute(&this->v[0].v[0], &a.v[0].v[0]);
        return *this;
      }

      #ifdef __CUDACC__
      /**
       * \brief Sets this matrix to the inverse of another matrix.
       *
       * \warning This function only works for \p m_ = \p n_ and will intentionally fail to compile in any other case.
       * \tparam ThreadGroup_ The active threadgroup.
       *
       * \param[in] tg
       * A reference to the threadgroups.
       * \param[in] a
       * The matrix whose inverse is to be stored in this matrix.
       * \param[in] det
       * The determinant of the matrix.
       *
       * \returns \p *this
       */
      template<typename ThreadGroup_, int sma_, int sna_>
      CUDA_HOST_DEVICE __forceinline__ Matrix& grouped_set_inverse(const ThreadGroup_& tg, const Matrix<T_, m_, n_, sma_, sna_>& a, const T_& det)
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and input matrix 'a' must be different objects");
        Intern::InverseHelper<m_, n_, sn_, sna_>::grouped_compute(tg, &this->v[0].v[0], &a.v[0].v[0], det);
        return *this;
      }
      #endif

      /**
       * \brief Sets this matrix to the cofactor matrix of another matrix.
       *
       * \f$ \mathrm{Cof}(A) \f$ is the <em>cofactor matrix</em> of \f$ A \f$. If \f$A \in \mathbb{R}^{n \times n}\f$,
       * define \f$A^{(i,j)} \in \mathbb{R}^{n-1 \times n-1}\f$ as the matrix obtained by deleting the \f$ i \f$th column
       * and \f$ j \f$th row of \f$ A \f$. Then
       * \f[
       *   \mathrm{Cof}(A)_{i,j} = (-1)^{i+j} \det(A^{(i,j)})
       * \f] and if \f$ A \in \mathrm{GL}_n\f$, \f$ \mathrm{Cof}(A) = \det(A) A^{-T} \f$.
       *
       * \warning This function only works for \p m_ = \p n_ and will intentionally fail to compile in any other case.
       *
       * \param[in] a
       * The matrix whose cofactor matrix is to be stored in this matrix.
       *
       * \returns \p *this
       */
      template<int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix& set_cofactor(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and input matrix 'a' must be different objects");
        Intern::CofactorHelper<m_, n_, sn_, sna_>::compute(&this->v[0].v[0], &a.v[0].v[0]);
        return *this;
      }

      /**
       * \brief Sets this matrix to the transpose of another matrix.
       *
       * \param[in] a
       * The matrix whose transpose is to be stored in this matrix.
       *
       * \returns \p *this
       */
      template<int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix& set_transpose(const Matrix<T_, n_, m_, sma_, sna_>& a)
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and input matrix 'a' must be different objects");

        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] = a.v[j][i];
          }
        }
        return *this;
      }

      /**
       * \brief Sets this matrix to the Gram matrix of another matrix.
       *
       * Let \e C denote \c this n-by-n matrix, and let \e A denote the l-by-n input matrix,
       * then this function computes:
       * \f[ C\leftarrow A^\top\cdot A \f]
       *
       * \param[in] a
       * The matrix whose Gram matrix is to be stored in this matrix.
       *
       * \returns \p *this
       */
      template<int l_, int sla_, int sna_>
      CUDA_HOST_DEVICE Matrix& set_gram(const Matrix<T_, l_, n_, sla_, sna_>& a)
      {
        static_assert(m_ == n_, "Gram matrices must be square");

        format();

        for(int k(0); k < l_; ++k)
        {
          for(int i(0); i < n_; ++i)
          {
            for(int j(0); j < n_; ++j)
            {
              v[i][j] += a.v[k][i] * a.v[k][j];
            }
          }
        }

        return *this;
      }

      /**
       * \brief Computes the scalar product of two vectors with this matrix.
       *
       * This function returns
       * \f[ x^\top\cdot A\cdot y = \sum_{i=0}^{m-1}\sum_{j=0}^{n-1} x_i\cdot A_{ij}\cdot y_j\f]
       *
       * \param[in] x
       * The left muliplicant vector of size \p m_.
       *
       * \param[in] y
       * The right multiplicand vector of size \p n_.
       *
       * \returns
       * The scalar product of \p x and \p y with this matrix.
       */
      template<int snx_, int sny_>
      CUDA_HOST_DEVICE DataType scalar_product(const Vector<T_, m_, snx_>& x, const Vector<T_, n_, sny_>& y) const
      {
        DataType r(DataType(0));
        for(int i(0); i < m_; ++i)
        {
          r += x[i] * dot(v[i], y);
        }
        return r;
      }

      /**
       * \brief Adds the outer product of two vectors onto the matrix.
       *
       * This function performs:
       * \f[ a_{ij} \leftarrow a_{ij} + \alpha x_i y_j \f]
       *
       * \param[in] x
       * The left multiplicand vector of size \p m_.
       *
       * \param[in] y
       * The right multiplicand vector of size \p n_.
       *
       * \param[in] alpha
       * The scaling factor for the outer product.
       *
       * \returns \c *this
       */
      template<int snx_, int sny_>
      CUDA_HOST_DEVICE Matrix& add_outer_product(
        const Vector<T_, m_, snx_>& x,
        const Vector<T_, n_, sny_>& y,
        const DataType alpha = DataType(1))
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] += alpha * x[i] * y[j];
          }
        }
        return *this;
      }

      /**
       * \brief Sets this matrix to the outer product of two vectors.
       *
       * This function performs:
       * \f[ a_{ij} \leftarrow x_i y_j \f]
       *
       * \param[in] x
       * The left multiplicand vector of size \p m_.
       *
       * \param[in] y
       * The right multiplicand vector of size \p n_.
       *
       * \returns \c *this
       */
      template<int snx_, int sny_>
      CUDA_HOST_DEVICE Matrix& set_outer_product(
        const Vector<T_, m_, snx_>& x,
        const Vector<T_, n_, sny_>& y)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] = x[i] * y[j];
          }
        }
        return *this;
      }

      /**
       * \brief Adds another scaled matrix onto this matrix.
       *
       * \param[in] alpha
       * The scaling parameter for the axpy.
       *
       * \param[in] a
       * The matrix to be added onto this matrix.
       *
       * \returns \c *this
       */
      template<int sma_, int sna_>
      CUDA_HOST_DEVICE Matrix& axpy(DataType alpha, const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] += alpha * a.v[i][j];
          }
        }
        return *this;
      }

      /**
       * \brief Adds a value onto the matrix's main diagonal.
       *
       * \param[in] alpha
       * The value that is to be added onto the main diagonal.
       */
      CUDA_HOST_DEVICE Matrix& add_scalar_main_diag(DataType alpha)
      {
        for(int i(0); (i < m_) && (i < n_); ++i)
          v[i][i] += alpha;
        return *this;
      }

      /**
       * \brief Adds the algebraic matrix-product of two other matrices onto this matrix.
       *
       * Let \e C denote \c this matrix, and let \e A denote the left m-by-l matrix and \e B the right l-by-n matrix,
       * then this function computes:
       * \f[ C\leftarrow C + \alpha A\cdot B \f]
       *
       * \param[in] a
       * The left m-by-l multiplicand matrix.
       *
       * \param[in] b
       * The right l-by-n multiplicand matrix.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int la_, int lb_, int sma_, int sna_, int smb_, int snb_>
      CUDA_HOST_DEVICE Matrix& add_mat_mat_mult(
        const Matrix<T_, m_, la_, sma_, sna_>& a,
        const Matrix<T_, lb_, n_, smb_, snb_>& b,
        DataType alpha = DataType(1))
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and multiplicand matrix 'a' must be different objects");
        ASSERTM((const void*)this != (const void*)&b, "result matrix and multiplicand matrix 'b' must be different objects");
        ASSERTM(la_ == lb_, "second dimension of a must be equal to first dimension of b");

        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DataType r(0);
            for(int k(0); k < la_; ++k)
            {
              r += a.v[i][k] * b.v[k][j];
            }
            v[i][j] += alpha * r;
          }
        }
        return *this;
      }

      /**
       * \brief Sets this matrix to the algebraic matrix-product of two other matrices.
       *
       * Let \e C denote \c this matrix, and let \e A denote the left m-by-l matrix and \e B the right l-by-n matrix,
       * then this function computes:
       * \f[ C\leftarrow A\cdot B \f]
       *
       * \param[in] a
       * The left m-by-l multiplicand matrix.
       *
       * \param[in] b
       * The right l-by-n multiplicand matrix.
       *
       * \returns \p *this
       */
      template<int la_, int lb_, int sma_, int sna_, int smb_, int snb_>
      CUDA_HOST_DEVICE Matrix& set_mat_mat_mult(const Matrix<T_, m_, la_, sma_, sna_>& a, const Matrix<T_, lb_, n_, smb_, snb_>& b)
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and multiplicand matrix 'a' must be different objects");
        ASSERTM((const void*)this != (const void*)&b, "result matrix and multiplicand matrix 'b' must be different objects");
        ASSERTM(la_ == lb_, "second dimension of a must be equal to first dimension of b");

        format();
        return add_mat_mat_mult(a, b);
      }

      /**
       * \brief Adds the algebraic matrix double-product of three other matrices onto this matrix.
       *
       * Let \e C denote \c this m-by-n matrix, \e B the left k-by-m input matrix, \e D the right l-by-n input matrix
       * and \e A the inner k-by-l input matrix, then this operation computes:
       *   \f[ C \leftarrow C + \alpha B^\top\cdot A\cdot D\f]
       *
       * \param[in] a
       * The inner k-by-l multiplicand matrix.
       *
       * \param[in] b
       * The left k-by-m multiplicand matrix.
       *
       * \param[in] d
       * The right l-by-n multiplicand matrix.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int k_, int l_, int sma_, int sna_, int smb_, int snb_, int smd_, int snd_>
      CUDA_HOST_DEVICE Matrix& add_double_mat_mult(
        const Matrix<T_, k_, l_, sma_, sna_>& a,
        const Matrix<T_, k_, m_, smb_, snb_>& b,
        const Matrix<T_, l_, n_, smd_, snd_>& d,
        DataType alpha = DataType(1))
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and multiplicand matrix 'a' must be different objects");
        ASSERTM((const void*)this != (const void*)&b, "result matrix and multiplicand matrix 'b' must be different objects");
        ASSERTM((const void*)this != (const void*)&d, "result matrix and multiplicand matrix 'd' must be different objects");

        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DataType r(0);
            for(int p(0); p < k_; ++p)
            {
              DataType t(0);
              for(int q(0); q < l_; ++q)
              {
                t += a(p,q) * d(q,j);
              }
              r += b(p,i)*t;
            }
            v[i][j] += alpha * r;
          }
        }
        return *this;
      }

      /**
       * \brief Sets this matrix to the algebraic matrix double-product of three other matrices.
       *
       * Let \e C denote \c this m-by-n matrix, \e B the left k-by-m input matrix, \e D the right l-by-n input matrix
       * and \e A the inner k-by-l input matrix, then this operation computes:
       *   \f[ C \leftarrow B^\top\cdot A\cdot D\f]
       *
       * \param[in] a
       * The inner k-by-l multiplicand matrix.
       *
       * \param[in] b
       * The left k-by-m multiplicand matrix.
       *
       * \param[in] d
       * The right l-by-n multiplicand matrix.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int k_, int l_, int sma_, int sna_, int smb_, int snb_, int smd_, int snd_>
      CUDA_HOST_DEVICE Matrix& set_double_mat_mult(
        const Matrix<T_, k_, l_, sma_, sna_>& a,
        const Matrix<T_, k_, m_, smb_, snb_>& b,
        const Matrix<T_, l_, n_, smd_, snd_>& d,
        T_ alpha = T_(1))
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&a, "result matrix and multiplicand matrix 'a' must be different objects");
        ASSERTM((const void*)this != (const void*)&b, "result matrix and multiplicand matrix 'b' must be different objects");
        ASSERTM((const void*)this != (const void*)&d, "result matrix and multiplicand matrix 'd' must be different objects");

        format();
        return add_double_mat_mult(a, b, d, alpha);
      }

      /**
       * \brief Adds the result of a vector-tensor left-product onto this matrix.
       *
       * Let \e A denote \c this m-by-n matrix, \e v the l-size input vector and \e T the l-by-m-by-b input tensor,
       * then this operation computes:
       * \f[ \forall i\in\{0,...,m-1\},j\in\{0,...,n-1\}:~ A_{ij} \leftarrow A_{ij} + \alpha \sum_{k=0}^{l-1} v_k\cdot T_{kij}\f]
       *
       * \note This function is used for the computation of second-order derivatives by the chain rule.
       *
       * \param[in] x
       * The l-size vector that serves as a left multiplicand.
       *
       * \param[in] t
       * The l-by-m-by-n tensor that serves as a right multiplicand.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int l_, int snv_, int slt_, int smt_, int snt_>
      CUDA_HOST_DEVICE Matrix& add_vec_tensor_mult(
        const Vector<T_, l_, snv_>& x,
        const Tensor3<T_, l_, m_, n_, slt_, smt_, snt_>& t,
        DataType alpha = DataType(1))
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DataType r(0);
            for(int k(0); k < l_; ++k)
            {
              r += x(k) * t(k,i,j);
            }
            v[i][j] += alpha * r;
          }
        }
        return *this;
      }

      /**
       * \brief Sets this matrix to the result of a vector-tensor left-product.
       *
       * Let \e A denote \c this m-by-n matrix, \e v the l-size input vector and \e T the l-by-m-by-b input tensor,
       * then this operation computes:
       * \f[ \forall i\in\{0,...,m-1\},j\in\{0,...,n-1\}:~ A_{ij} \leftarrow \alpha \sum_{k=0}^{l-1} v_k\cdot T_{kij}\f]
       *
       * \note This function is used for the computation of second-order derivatives by the chain rule.
       *
       * \param[in] x
       * The l-size vector that serves as a left multiplicand.
       *
       * \param[in] t
       * The l-by-m-by-n tensor that serves as a right multiplicand.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int l_, int snv_, int slt_, int smt_, int snt_>
      CUDA_HOST_DEVICE Matrix& set_vec_tensor_mult(
        const Vector<T_, l_, snv_>& x,
        const Tensor3<T_, l_, m_, n_, slt_, smt_, snt_>& t,
        DataType alpha = DataType(1))
      {
        format();
        return add_vec_tensor_mult(x, t, alpha);
      }

      /**
       * \brief Sets this matrix to the identity matrix.
       *
       * \returns \p *this
       */
      CUDA_HOST_DEVICE Matrix& set_identity()
      {
        for(int i(0); i < m_; ++i)
          for(int j(0); j < n_; ++j)
            v[i][j] = (i == j ? T_(1) : T_(0));
        return *this;
      }

      /**
       * \brief Sets this matrix to a 2D rotation matrix.
       *
       * \param[in] angle
       * Specifies the rotation angle in radians
       *
       * \returns \p *this
       */
      CUDA_HOST_DEVICE Matrix& set_rotation_2d(T_ angle)
      {
        static_assert((m_ == 2) && (n_ == 2), "this function works only for 2x2 matrices");
        #ifndef __CUDACC__
        v[0][0] =  (v[1][1] = Math::cos(angle));
        v[0][1] = -(v[1][0] = Math::sin(angle));
        #else
        v[0][0] =  (v[1][1] = CudaMath::cuda_cos(angle));
        v[0][1] = -(v[1][0] = CudaMath::cuda_sin(angle));
        #endif
        return *this;
      }

      /**
       * \brief Sets this matrix to a 3D yaw-pitch-roll rotation matrix.
       *
       * This function sets a 3D rotation matrix by specifying the yaw, pitch and roll angles.
       *
       * The convention for the rotation is <b>Z-Y'-X''</b>, also known as <em>nautical angles</em>.
       *
       * \see
       * - https://en.wikipedia.org/wiki/Euler_angles
       * - https://en.wikipedia.org/wiki/Aircraft_principal_axes
       * - http://planning.cs.uiuc.edu/node102.html
       *
       * \returns \p *this
       */
      CUDA_HOST_DEVICE Matrix& set_rotation_3d(T_ yaw, T_ pitch, T_ roll)
      {
        static_assert((m_ == 3) && (n_ == 3), "this function works only for 3x3 matrices");
        #ifndef __CUDACC__
        const T_ cy = Math::cos(yaw);
        const T_ sy = Math::sin(yaw);
        const T_ cp = Math::cos(pitch);
        const T_ sp = Math::sin(pitch);
        const T_ cr = Math::cos(roll);
        const T_ sr = Math::sin(roll);
        #else
        const T_ cy = CudaMath::cuda_cos(yaw);
        const T_ sy = CudaMath::cuda_sin(yaw);
        const T_ cp = CudaMath::cuda_cos(pitch);
        const T_ sp = CudaMath::cuda_sin(pitch);
        const T_ cr = CudaMath::cuda_cos(roll);
        const T_ sr = CudaMath::cuda_sin(roll);
        #endif
        v[0][0] = cy*cp;
        v[0][1] = cy*sp*sr - sy*cr;
        v[0][2] = cy*sp*cr + sy*sr;
        v[1][0] = sy*cp;
        v[1][1] = sy*sp*sr + cy*cr;
        v[1][2] = sy*sp*cr - cy*sr;
        v[2][0] = -sp;
        v[2][1] = cp*sr;
        v[2][2] = cp*cr;
        return *this;
      }

      /**
       * \brief Returns a null-matrix.
       */
      CUDA_HOST_DEVICE static Matrix null()
      {
        return Matrix(DataType(0));
      }

      /**
       * \brief Tiny::Matrix streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] A The matrix to be streamed.
       */
      CUDA_HOST friend std::ostream & operator<< (std::ostream & lhs, const Matrix& A)
      {
        for (int i(0) ; i < m-1 ; ++i)
        {
          lhs << A[i] << std::endl;
        }
        lhs << A[m-1];

        return lhs;
      }
    }; // class Matrix

#ifdef DOXYGEN
    /**
     * \brief Computes the positively oriented orthogonal vector to the columns of a m_ x (m_-1) Matrix
     *
     * \tparam T_
     * Data type.
     *
     * \tparam m_
     * Number of rows of the matrix.
     *
     * \tparam sm_
     * Row stride of the matrix.
     *
     * \tparam sn_
     * Column stride of the matrix.
     *
     * \param[in] tau
     * The m_ x (m_-1) matrix.
     *
     * \returns
     * A vector that is orthogonal to the m_-1 columns of the input matrix, but not normalized.
     *
     * \note So far, this is only implemented for m_ = 2,3.
     */
    template<typename T_, int m_, int sm_, int sn_>
    Vector<T_, m_> orthogonal(const Matrix<T_, m_, m_-1, sm_, sn_>& tau);
#endif

    /// \cond internal
    template<typename T_, int sm_, int sn_>
    CUDA_HOST_DEVICE Vector<T_, 2> orthogonal(const Matrix<T_, 2, 1, sm_, sn_>& tau)
    {
      Vector<T_, 2, sm_> nu(T_(0));

      // 2d "cross" product. The sign has to be on the second component so the input is rotated in negative direction
      nu[0] =  tau[1][0];
      nu[1] = -tau[0][0];

      return nu;
    }

    template<typename T_, int sm_, int sn_>
    CUDA_HOST_DEVICE Vector<T_, 3> orthogonal(const Matrix<T_, 3, 2, sm_, sn_>& tau)
    {
      Vector<T_, 3, sm_> nu(T_(0));

      // 3d cross product
      nu[0] = tau[1][0]*tau[2][1] - tau[2][0]*tau[1][1];
      nu[1] = tau[2][0]*tau[0][1] - tau[0][0]*tau[2][1];
      nu[2] = tau[0][0]*tau[1][1] - tau[1][0]*tau[0][1];

      return nu;
    }
    /// \endcond

    /// matrix-vector-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_, int sx_>
    CUDA_HOST_DEVICE inline Vector<T_, m_> operator*(const Matrix<T_, m_, n_, sm_, sn_>& a, const Vector<T_, n_, sx_>& x)
    {
      return Vector<T_, m_>().set_mat_vec_mult(a, x);
    }

    /// vector-matrix-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_, int sx_>
    CUDA_HOST_DEVICE inline Vector<T_, n_> operator*(const Vector<T_, m_, sx_>& x, const Matrix<T_, m_, n_, sm_, sn_>& a)
    {
      return Vector<T_, n_>().set_vec_mat_mult(x, a);
    }

    /// scalar-left-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline Matrix<T_, m_, n_> operator*(typename Matrix<T_, m_, n_>::DataType alpha, const Matrix<T_, m_, n_, sm_, sn_>& a)
    {
      return Matrix<T_, m_, n_>(a) *= alpha;
    }

    /// scalar-right-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline Matrix<T_, m_, n_, sm_, sn_> operator*(const Matrix<T_, m_, n_, sm_, sn_>& a, typename Matrix<T_, m_, n_>::DataType alpha)
    {
      return Matrix<T_, m_, n_>(a) *= alpha;
    }

    /// algebraic matrix-matrix-multiply operator
    template<typename T_, int m_, int n_, int l_, int sma_, int sna_, int smb_, int snb_>
    CUDA_HOST_DEVICE inline Matrix<T_, m_, n_> operator*(const Matrix<T_, m_, l_, sma_, sna_>& a, const Matrix<T_, l_, n_, smb_, snb_>& b)
    {
      return Matrix<T_, m_, n_>().set_mat_mat_mult(a, b);
    }

    /// matrix addition operator
    template<typename T_, int m_, int n_,int sma_, int sna_, int smb_, int snb_>
    CUDA_HOST_DEVICE inline Matrix<T_, m_, n_> operator+(const Matrix<T_, m_, n_, sma_, sna_>& a, const Matrix<T_, m_, n_, smb_, snb_>& b)
    {
      return Matrix<T_, m_, n_>(a) += b;
    }

    /// matrix subtraction operator
    template<typename T_, int m_, int n_,int sma_, int sna_, int smb_, int snb_>
    CUDA_HOST_DEVICE inline Matrix<T_, m_, n_> operator-(const Matrix<T_, m_, n_, sma_, sna_>& a, const Matrix<T_, m_, n_, smb_, snb_>& b)
    {
      return Matrix<T_, m_, n_>(a) -= b;
    }

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    // Tiny Tensor3 implementation
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
    class Tensor3
    {
      static_assert(l_ > 0, "invalid tube count");
      static_assert(m_ > 0, "invalid row count");
      static_assert(n_ > 0, "invalid column count");
      static_assert(sl_ >= l_, "invalid tube stride");
      static_assert(sm_ >= m_, "invalid row stride");
      static_assert(sn_ >= n_, "invalid column stride");

    public:
      /// the tube count of the tensor
      static constexpr int l = l_;
      /// the row count of the tensor
      static constexpr int m = m_;
      /// the column count of the tensor
      static constexpr int n = n_;
      /// the tube stride of the tensor
      static constexpr int sl = sl_;
      /// the row stride of the tensor
      static constexpr int sm = sm_;
      /// the column stride of the tensor
      static constexpr int sn = sn_;

      /// the data type of the tensor
      typedef T_ ValueType;
      /// The basic data type buried in the lowest level of the vector
      typedef typename Intern::DataTypeExtractor<ValueType>::MyDataType DataType;

      /// Type of tensor data; that's an array of matrices
      typedef Matrix<T_, m_, n_, sm_, sn_> PlaneType;
      /// Actual tensor data
      PlaneType v[sl_];

      /// default constructor
      CUDA_HOST_DEVICE Tensor3()
      {
      }

      /// value-assignment constructor
      CUDA_HOST_DEVICE explicit Tensor3(DataType value)
      {
        for(int i(0); i < l_; ++i)
          v[i] = value;
      }

      /**
       * \brief Initializer list constructor
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE explicit Tensor3(const std::initializer_list<Tx_>& x)
      {
        XASSERTM(std::size_t(l_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < l_; ++i, ++it)
          v[i] = *it;
      }

      /**
       * \brief Initializer list constructor
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE explicit Tensor3(const std::initializer_list<std::initializer_list<std::initializer_list<Tx_>>>& x)
      {
        XASSERTM(std::size_t(l_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < l_; ++i, ++it)
          v[i] = *it;
      }

      /// copy-constructor
      template<int sla_, int sma_, int sna_>
      CUDA_HOST_DEVICE Tensor3(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] = a.v[i];
      }

      /// value-assignment operator
      CUDA_HOST_DEVICE Tensor3& operator=(DataType value)
      {
        for(int i(0); i < l_; ++i)
          v[i] = value;
        return *this;
      }

      /// copy-assignment operator
      template<int sla_, int sma_, int sna_>
      CUDA_HOST_DEVICE Tensor3& operator=(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] = a.v[i];
        return *this;
      }

      /**
       * \brief Initializer list assignment operator
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       *
       * \returns *this
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE Tensor3& operator=(const std::initializer_list<Tx_>& x)
      {
        XASSERTM(std::size_t(l_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < l_; ++i, ++it)
          v[i] = *it;
        return *this;
      }

      /**
       * \brief Initializer list assignment operator
       *
       * \param[in] x
       * The initializer list whose elements are to be assigned.
       *
       * \returns *this
       */
      template<typename Tx_>
      CUDA_HOST_DEVICE Tensor3& operator=(const std::initializer_list<std::initializer_list<std::initializer_list<Tx_>>>& x)
      {
        XASSERTM(std::size_t(l_) == x.size(), "invalid initializer list size");
        auto it(x.begin());
        for(int i(0); i < l_; ++i, ++it)
          v[i] = *it;
        return *this;
      }

      /**
       * \brief Access operator
       *
       * \param[in] h,i,j
       * The indices of the tensor entry that is to be returned.
       *
       * \returns
       * A (const) reference to the tensor entry at position (h,i,j).
       */
      CUDA_HOST_DEVICE T_& operator()(int h, int i, int j)
      {
        ASSERTM( (h >= 0) && (h < l_), "index h out-of-bounds");
        ASSERTM( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERTM( (j >= 0) && (j < n_), "index j out-of-bounds");
        return v[h](i,j);
      }

      /** \copydoc operator()() */
      CUDA_HOST_DEVICE const T_& operator()(int h, int i, int j) const
      {
        ASSERTM( (h >= 0) && (h < l_), "index h out-of-bounds");
        ASSERTM( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERTM( (j >= 0) && (j < n_), "index j out-of-bounds");
        return v[h](i,j);
      }

      /**
       * \brief Plane-Access operator.
       *
       * \param[in] h
       * The index of the plane that is to be returned.
       *
       * \returns
       * A (const) reference to the matrix representing the <c>h</c>-th plane of the tensor.
       */
      CUDA_HOST_DEVICE PlaneType& operator[](int h)
      {
        ASSERTM( (h >= 0) && (h < l_), "index h out-of-bounds");
        return v[h];
      }

      /** \copydoc operator[]() */
      CUDA_HOST_DEVICE const PlaneType& operator[](int h) const
      {
        ASSERTM( (h >= 0) && (h < l_), "index h out-of-bounds");
        return v[h];
      }

      /**
       * \brief Size-cast function.
       *
       * This function casts this tensor's reference to another size.
       *
       * \tparam ll_, mm_, nn_
       * The dimensions of the casted matrix. Must be 0 < \p ll_ <= \p sl_, 0 < \p mm_ <= \p sm_
       * and 0 < \p nn_ <= \p sn_.
       *
       * \returns
       * A casted (const) reference of \p *this.
       */
      template<int ll_, int mm_, int nn_>
      CUDA_HOST_DEVICE Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>& size_cast()
      {
        static_assert((ll_ >= 0) && (ll_ <= sl_), "invalid cast tube count");
        static_assert((mm_ >= 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ >= 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>&>(*this);
      }

      /** \copydoc size_cast() */
      template<int ll_, int mm_, int nn_>
      CUDA_HOST_DEVICE const Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>& size_cast() const
      {
        static_assert((ll_ >= 0) && (ll_ <= sl_), "invalid cast tube count");
        static_assert((mm_ >= 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ >= 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>&>(*this);
      }

      /// scalar right-multiply-by operator
      CUDA_HOST_DEVICE Tensor3& operator*=(DataType alpha)
      {
        for(int i(0); i < l_; ++i)
          v[i] *= alpha;
        return *this;
      }

      /// tensor component-wise addition operator
      template<int sla_, int sma_, int sna_>
      CUDA_HOST_DEVICE Tensor3& operator+=(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] += a.v[i];
        return *this;
      }

      /// tensor component-wise subtraction operator
      template<int sla_, int sma_, int sna_>
      CUDA_HOST_DEVICE Tensor3& operator-=(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] -= a.v[i];
        return *this;
      }

      /// formats the tensor
      CUDA_HOST_DEVICE void format(DataType alpha = DataType(0))
      {
        (*this) = alpha;
      }

      /**
       * \brief Adds the result of a matrix-tensor product onto this tensor.
       *
       * Let \e K denote this tensor, \e A the input matrix and \e T the input tensor, then this operation computes:
       * \f[ \forall h\in\{0,...,l-1\}, i\in\{0,...,m-1\},j\in\{0,...,n-1\}:~ K_{hij} \leftarrow K_{hij} +
       *     \alpha \sum_{p=0}^{k-1} A_{hp}\cdot T_{pij}\f]
       *
       * \note This function is used for the computation of second-order derivatives by the chain rule.
       *
       * \param[in] a
       * The l-by-k matrix that serves as the left multiplicand.
       *
       * \param[in] t
       * The k-by-m-by-n tensor that serves as the right multiplicand.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int k_, int sma_, int sna_, int slt_, int smt_, int snt_>
      CUDA_HOST_DEVICE Tensor3& add_mat_tensor_mult(
        const Matrix<T_, l_, k_, sma_, sna_>& a,
        const Tensor3<T_, k_, m_, n_, slt_, smt_, snt_>& t,
        DataType alpha = DataType(1))
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&t, "result tensor and multiplicand tensor 't' must be different objects");

        for(int h(0); h < l_; ++h)
        {
          for(int i(0); i < m_; ++i)
          {
            for(int j(0); j < n_; ++j)
            {
              DataType r(0);
              for(int p(0); p < k_; ++p)
              {
                r += a(h,p) * t(p,i,j);
              }
              operator()(h,i,j) += alpha * r;
            }
          }
        }
        return *this;
      }

      /**
       * \brief Adds the result of a matrix-tensor-matrix double-product onto this tensor.
       *
       * Let \e K denote this l-by-m-by-n tensor, \e T the l-by-m'-by-n' input tensor,
       * \e B the m'-by-m and \e D the n'-by-n input matrices, then this operation computes:
       * \f[ \forall h\in\{0,...,l-1\}, i\in\{0,...,m-1\}, j\in\{0,...,n-1\}:~ K_{hij} \leftarrow K_{hij} +
       *     \alpha \sum_{p=0}^{m'-1}\sum_{q=0}^{n'-1} T_{hpq}\cdot B_{pi}\cdot D_{qj} \f]
       *
       * \note This function is used for the computation of second-order derivatives by the chain rule.
       *
       * \param[in] t
       * The l-by-m'-by-n' tensor that serves as the inner multiplicand.
       *
       * \param[in] b
       * The m'-by-m matrix that serves as the left multiplicand.
       *
       * \param[in] d
       * The n'-by-n matrix that serves as the right multiplicand.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<
        int lt_, int mt_, int nt_, // input tensor dimensions
        int slt_, int smt_, int snt_, // input tensor strides
        int smb_, int snb_, int smd_, int snd_> // input matrix strides
      CUDA_HOST_DEVICE Tensor3& add_double_mat_mult(
        const Tensor3<T_, lt_, mt_, nt_, slt_, smt_, snt_>& t,
        const Matrix<T_, nt_, n_, smb_, snb_>& b,
        const Matrix<T_, mt_, m_, smd_, snd_>& d,
        DataType alpha = DataType(1))
      {
        // we have to compare void* addresses here, because we might get a type mismatch error otherwise
        ASSERTM((const void*)this != (const void*)&t, "result tensor and multiplicand tensor 't' must be different objects");

        for(int h(0); h < l_; ++h)
        {
          for(int i(0); i < m_; ++i)
          {
            for(int j(0); j < n_; ++j)
            {
              DataType r(0);
              for(int p(0); p < mt_; ++p)
              {
                for(int q(0); q < nt_; ++q)
                {
                  r += t(h,p,q) * b(p,i) * d(q,j);
                }
              }
              operator()(h,i,j) += alpha * r;
            }
          }
        }
        return *this;
      }

      /**
       * \brief Adds the result of a vector-matrix outer product onto this tensor.
       *
       * Let \e K denote this tensor, \e x the input vector and \e A the input matrix, then this operation computes:
       * \f[ \forall h\in\{0,...,l-1\}, i\in\{0,...,m-1\},j\in\{0,...,n-1\}:~ K_{hij} \leftarrow K_{hij} +
       *     \alpha x_h A_{ij}\f]
       *
       * \note This function is used for the computation of vector-valued hessian tensors.
       *
       * \param[in] x
       * The l-length vector that serves as the left multiplicand.
       *
       * \param[in] a
       * The m-by-n matrix that serves as the right multiplicand.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int slx_, int sma_, int sna_>
      CUDA_HOST_DEVICE Tensor3& add_vec_mat_outer_product(
        const Vector<T_, l_, slx_>& x,
        const Matrix<T_, m_, n_, sma_, sna_>& a,
        DataType alpha = DataType(1))
      {
        for(int h(0); h < l_; ++h)
        {
          for(int i(0); i < m_; ++i)
          {
            for(int j(0); j < n_; ++j)
            {
              operator()(h,i,j) += alpha * x(h) * a(i, j);
            }
          }
        }
        return *this;
      }

      /**
       * \brief Returns a null-tensor.
       */
      CUDA_HOST_DEVICE static Tensor3 null()
      {
        return Tensor3(DataType(0));
      }
    }; // class Tensor3<...>

    /// scalar left-multiply operator
    template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline Tensor3<T_, l_, m_, n_, sl_, sm_, sn_> operator*(
      typename Tensor3<T_, l_, m_, n_>::DataType alpha, const Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>& a)
    {
      return Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>(a) *= alpha;
    }

    /// scalar right-multiply operator
    template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline Tensor3<T_, l_, m_, n_, sl_, sm_, sn_> operator*(
      const Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>& a, typename Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>::DataType alpha)
    {
      return Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>(a) *= alpha;
    }

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    // Various helper functions
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    /**
     * \brief Computes the dot-product of two scalars.
     *
     * \param[in] a,b
     * The scalar whose product is to be computed.
     *
     * \returns The product of \p a and \p b.
     *
     * \note
     * This overload is disabled by SFINAE for Tiny::Vector, Tiny::Matrix and Tiny::Tensor3 types.
     */
    template<typename T_, typename std::enable_if<Intern::DataTypeExtractor<T_>::level == 0, bool>::type = true>
    CUDA_HOST_DEVICE inline T_ dot(const T_& a, const T_& b)
    {
      return a*b;
    }

    /**
     * \brief Computes the dot-product of two vectors.
     *
     * This function returns the dot-product of \e a and \e b:
     * \f[\sum_{k=0}^{n-1} a_k\cdot b_k\f]
     *
     * \param[in] a,b
     * The vectors whose dot-product is to be computed.
     *
     * \returns The dot-product of \p a and \p b.
     */
    template<typename T_, int n_, int sa_, int sb_>
    CUDA_HOST_DEVICE inline typename Vector<T_, n_>::DataType dot(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      typename Vector<T_, n_>::DataType r(0);

      for(int i(0); i < n_; ++i)
      {
        r += Tiny::dot(a.v[i], b.v[i]);
      }
      return r;
    }

    /**
     * \brief Computes the dot-product of two matrices.
     *
     * This function returns the dot-product of \e A and \e B:
     * \f[\sum_{i=0}^{m-1} \sum_{j=0}^{n-1} a_{ij} \cdot b_{ij}\f]
     *
     * \param[in] a,b
     * The matrices whose dot-product is to be computed.
     *
     * \returns The dot-product of \p a and \p b.
     */
    template<typename T_, int m_, int n_, int sma_, int sna_, int smb_, int snb_>
    CUDA_HOST_DEVICE inline typename Matrix<T_, m_, n_>::DataType dot(const Matrix<T_, m_, n_, sma_, sna_>& a, const Matrix<T_, m_, n_, smb_, snb_>& b)
    {
      typename Matrix<T_, m_, n_>::DataType r(0);
      for(int i(0); i < m_; ++i)
      {
        r += Tiny::dot(a.v[i], b.v[i]);
      }
      return r;
    }

    /**
     * \brief Computes the dot-product of two tensor3 objects.
     *
     * \param[in] a,b
     * The tensor3 whose dot-product is to be computed.
     *
     * \returns The dot-product of \p a and \p b.
     */
    template<typename T_, int l_, int m_, int n_, int sla_ ,int sma_, int sna_, int slb_, int smb_, int snb_>
    CUDA_HOST_DEVICE inline typename Tensor3<T_, l_, m_, n_>::DataType dot(
      const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a,
      const Tensor3<T_, l_, m_, n_, slb_, smb_, snb_>& b)
    {
      typename Tensor3<T_, l_, m_, n_>::DataType r(0);
      for(int i(0); i < l_; ++i)
      {
        r += Tiny::dot(a.v[i], b.v[i]);
      }
      return r;
    }

    /**
     * \brief Adds a scaled identity onto a scalar.
     *
     * \param[in,out] x
     * The scalar onto which to add to
     *
     * \param[in] alpha
     * The scalar to be added onto \p x
     */
    template<typename T_>
    CUDA_HOST_DEVICE inline void add_id(T_& x, const T_& alpha)
    {
      x += alpha;
    }

    /**
     * \brief Adds a scaled identity onto a vector.
     *
     * \param[in,out] x
     * The vector onto which to add to
     *
     * \param[in] alpha
     * The scalar to be added onto \p x
     */
    template<typename T_, int n_, int sn_>
    CUDA_HOST_DEVICE inline void add_id(Vector<T_, n_, sn_>& x, const typename Vector<T_, n_, sn_>::DataType& alpha)
    {
      for(int i(0); i < n_; ++i)
        add_id(x(i), alpha);
    }

    /**
     * \brief Adds a scaled identity onto a matrix.
     *
     * \param[in,out] x
     * The matrix onto which to add to
     *
     * \param[in] alpha
     * The scalar to be added onto the main diagonal of \p x
     */
    template<typename T_, int n_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline void add_id(Matrix<T_, n_, n_, sm_, sn_>& x, const typename Matrix<T_, n_, n_, sm_, sn_>::DataType& alpha)
    {
      for(int i(0); i < n_; ++i)
        add_id(x(i,i), alpha);
    }

    /**
     * \brief Adds a scaled identity onto a tensor3.
     *
     * \param[in,out] x
     * The tensor3 onto which to add to
     *
     * \param[in] alpha
     * The scalar to be added onto the main diagonal of \p x
     */
    template<typename T_, int n_, int sl_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline void add_id(Tensor3<T_, n_, n_, n_, sl_, sm_, sn_>& x, const typename Tensor3<T_, n_, n_, n_, sl_, sm_, sn_>::DataType& alpha)
    {
      for(int i(0); i < n_; ++i)
        add_id(x(i,i,i), alpha);
    }

    /**
     * \brief Performs an AXPY of two scalars
     *
     * This function performs: y += alpha*x
     *
     * \param[in,out] y
     * The object that receives the AXPY result
     *
     * \param[in] x
     * The object that is to be added onto \p y
     *
     * \param[in] alpha
     * The scaling factor for \p x
     */
    template<typename T_>
    CUDA_HOST_DEVICE inline void axpy(T_& y, const T_& x, const T_& alpha)
    {
      y += alpha*x;
    }

    /**
     * \brief Performs an AXPY of two vectors
     *
     * This function performs: y += alpha*x
     *
     * \param[in,out] y
     * The object that receives the AXPY result
     *
     * \param[in] x
     * The object that is to be added onto \p y
     *
     * \param[in] alpha
     * The scaling factor for \p x
     */
    template<typename T_, int n_, int sn_>
    CUDA_HOST_DEVICE inline void axpy(
      Vector<T_, n_, sn_>& y,
      const Vector<T_, n_, sn_>& x,
      const typename Vector<T_, n_, sn_>::DataType& alpha)
    {
      for(int i(0); i < n_; ++i)
        axpy(y.v[i], x.v[i], alpha);
    }

    /**
     * \brief Performs an AXPY of two matrices
     *
     * This function performs: y += alpha*x
     *
     * \param[in,out] y
     * The object that receives the AXPY result
     *
     * \param[in] x
     * The object that is to be added onto \p y
     *
     * \param[in] alpha
     * The scaling factor for \p x
     */
    template<typename T_, int m_, int n_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline void axpy(
      Matrix<T_, m_, n_, sm_, sn_>& y,
      const Matrix<T_, m_, n_, sm_, sn_>& x,
      const typename Matrix<T_, m_, n_, sm_, sn_>::DataType& alpha)
    {
      for(int i(0); i < m_; ++i)
        axpy(y.v[i], x.v[i], alpha);
    }

    /**
     * \brief Performs an AXPY of two tensor3
     *
     * This function performs: y += alpha*x
     *
     * \param[in,out] y
     * The object that receives the AXPY result
     *
     * \param[in] x
     * The object that is to be added onto \p y
     *
     * \param[in] alpha
     * The scaling factor for \p x
     */
    template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
    CUDA_HOST_DEVICE inline void axpy(
      Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>& y,
      const Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>& x,
      const typename Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>::DataType& alpha)
    {
      for(int i(0); i < l_; ++i)
        axpy(y.v[i], x.v[i], alpha);
    }

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    // Internal helpers implementation
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /// \cond internal
    namespace Intern
    {
      // generic square matrix inversion:
      template<int n_, int sna_>
      struct DetHelper<n_,n_,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // create temporary copy of a and a pivot array
          T_ b[n_*n_];
          int p[n_];

          // copy matrix a to b
          for (int i(0); i < n_; ++i)
          {
            for (int j(0); j < n_; ++j)
            {
              b[i*n_+j] = a[i*sna_+j];
            }
          }

          // perform matrix inversion which returns the determinant
          #ifndef __CUDACC__
          const T_ det = Math::invert_matrix(n_, n_, b, p);
          #else
          const T_ det = CudaMath::cuda_invert_matrix(n_, n_, b, p);
          #endif


          // if the returned value is not normal, we can assume that the matrix is singular
          #ifndef __CUDACC__
          return Math::isnormal(det) ? det : T_(0);
          #else
          return CudaMath::cuda_isnormal(det) ? det : T_(0);
          #endif
        }
      };

      template<int sna_>
      struct DetHelper<1,1,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          return a[0];
        }
      };

      template<int sna_>
      struct DetHelper<2,2,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          return a[0*sna_+0]*a[1*sna_+1] - a[0*sna_+1]*a[1*sna_+0];
        }
      };

      template<int sna_>
      struct DetHelper<3,3,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          return  a[0*sna_+0]*(a[1*sna_+1]*a[2*sna_+2] - a[1*sna_+2]*a[2*sna_+1])
                + a[0*sna_+1]*(a[1*sna_+2]*a[2*sna_+0] - a[1*sna_+0]*a[2*sna_+2])
                + a[0*sna_+2]*(a[1*sna_+0]*a[2*sna_+1] - a[1*sna_+1]*a[2*sna_+0]);
        }
      };

      template<int sna_>
      struct DetHelper<4,4,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // 2x2 determinants of rows 3-4
          T_ w[6] =
          {
            a[2*sna_+0]*a[3*sna_+1] - a[2*sna_+1]*a[3*sna_+0],
            a[2*sna_+0]*a[3*sna_+2] - a[2*sna_+2]*a[3*sna_+0],
            a[2*sna_+0]*a[3*sna_+3] - a[2*sna_+3]*a[3*sna_+0],
            a[2*sna_+1]*a[3*sna_+2] - a[2*sna_+2]*a[3*sna_+1],
            a[2*sna_+1]*a[3*sna_+3] - a[2*sna_+3]*a[3*sna_+1],
            a[2*sna_+2]*a[3*sna_+3] - a[2*sna_+3]*a[3*sna_+2]
          };

          return
            + a[0*sna_+0] * (a[1*sna_+1]*w[5] - a[1*sna_+2]*w[4] + a[1*sna_+3]*w[3])
            - a[0*sna_+1] * (a[1*sna_+0]*w[5] - a[1*sna_+2]*w[2] + a[1*sna_+3]*w[1])
            + a[0*sna_+2] * (a[1*sna_+0]*w[4] - a[1*sna_+1]*w[2] + a[1*sna_+3]*w[0])
            - a[0*sna_+3] * (a[1*sna_+0]*w[3] - a[1*sna_+1]*w[1] + a[1*sna_+2]*w[0]);
        }
      };

      template<int sna_>
      struct DetHelper<5,5,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // 2x2 determinants of rows 4-5
          T_ v[10] =
          {
            a[3*sna_+0]*a[4*sna_+1] - a[3*sna_+1]*a[4*sna_+0],
            a[3*sna_+0]*a[4*sna_+2] - a[3*sna_+2]*a[4*sna_+0],
            a[3*sna_+0]*a[4*sna_+3] - a[3*sna_+3]*a[4*sna_+0],
            a[3*sna_+0]*a[4*sna_+4] - a[3*sna_+4]*a[4*sna_+0],
            a[3*sna_+1]*a[4*sna_+2] - a[3*sna_+2]*a[4*sna_+1],
            a[3*sna_+1]*a[4*sna_+3] - a[3*sna_+3]*a[4*sna_+1],
            a[3*sna_+1]*a[4*sna_+4] - a[3*sna_+4]*a[4*sna_+1],
            a[3*sna_+2]*a[4*sna_+3] - a[3*sna_+3]*a[4*sna_+2],
            a[3*sna_+2]*a[4*sna_+4] - a[3*sna_+4]*a[4*sna_+2],
            a[3*sna_+3]*a[4*sna_+4] - a[3*sna_+4]*a[4*sna_+3]
          };
          // 3x3 determinants of rows 3-4-5
          T_ w[10] =
          {
            a[2*sna_+0]*v[4] - a[2*sna_+1]*v[1] + a[2*sna_+2]*v[0],
            a[2*sna_+0]*v[5] - a[2*sna_+1]*v[2] + a[2*sna_+3]*v[0],
            a[2*sna_+0]*v[6] - a[2*sna_+1]*v[3] + a[2*sna_+4]*v[0],
            a[2*sna_+0]*v[7] - a[2*sna_+2]*v[2] + a[2*sna_+3]*v[1],
            a[2*sna_+0]*v[8] - a[2*sna_+2]*v[3] + a[2*sna_+4]*v[1],
            a[2*sna_+0]*v[9] - a[2*sna_+3]*v[3] + a[2*sna_+4]*v[2],
            a[2*sna_+1]*v[7] - a[2*sna_+2]*v[5] + a[2*sna_+3]*v[4],
            a[2*sna_+1]*v[8] - a[2*sna_+2]*v[6] + a[2*sna_+4]*v[4],
            a[2*sna_+1]*v[9] - a[2*sna_+3]*v[6] + a[2*sna_+4]*v[5],
            a[2*sna_+2]*v[9] - a[2*sna_+3]*v[8] + a[2*sna_+4]*v[7]
          };

          return
            + a[0*sna_+0]*(a[1*sna_+1]*w[9] - a[1*sna_+2]*w[8] + a[1*sna_+3]*w[7] - a[1*sna_+4]*w[6])
            - a[0*sna_+1]*(a[1*sna_+0]*w[9] - a[1*sna_+2]*w[5] + a[1*sna_+3]*w[4] - a[1*sna_+4]*w[3])
            + a[0*sna_+2]*(a[1*sna_+0]*w[8] - a[1*sna_+1]*w[5] + a[1*sna_+3]*w[2] - a[1*sna_+4]*w[1])
            - a[0*sna_+3]*(a[1*sna_+0]*w[7] - a[1*sna_+1]*w[4] + a[1*sna_+2]*w[2] - a[1*sna_+4]*w[0])
            + a[0*sna_+4]*(a[1*sna_+0]*w[6] - a[1*sna_+1]*w[3] + a[1*sna_+2]*w[1] - a[1*sna_+3]*w[0]);
        }
      };

      template<int sna_>
      struct DetHelper<6,6,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // 2x2 determinants of rows 5-6
          T_ v[15] =
          {
            a[4*sna_+0]*a[5*sna_+1] - a[4*sna_+1]*a[5*sna_+0],
            a[4*sna_+0]*a[5*sna_+2] - a[4*sna_+2]*a[5*sna_+0],
            a[4*sna_+0]*a[5*sna_+3] - a[4*sna_+3]*a[5*sna_+0],
            a[4*sna_+0]*a[5*sna_+4] - a[4*sna_+4]*a[5*sna_+0],
            a[4*sna_+0]*a[5*sna_+5] - a[4*sna_+5]*a[5*sna_+0],
            a[4*sna_+1]*a[5*sna_+2] - a[4*sna_+2]*a[5*sna_+1],
            a[4*sna_+1]*a[5*sna_+3] - a[4*sna_+3]*a[5*sna_+1],
            a[4*sna_+1]*a[5*sna_+4] - a[4*sna_+4]*a[5*sna_+1],
            a[4*sna_+1]*a[5*sna_+5] - a[4*sna_+5]*a[5*sna_+1],
            a[4*sna_+2]*a[5*sna_+3] - a[4*sna_+3]*a[5*sna_+2],
            a[4*sna_+2]*a[5*sna_+4] - a[4*sna_+4]*a[5*sna_+2],
            a[4*sna_+2]*a[5*sna_+5] - a[4*sna_+5]*a[5*sna_+2],
            a[4*sna_+3]*a[5*sna_+4] - a[4*sna_+4]*a[5*sna_+3],
            a[4*sna_+3]*a[5*sna_+5] - a[4*sna_+5]*a[5*sna_+3],
            a[4*sna_+4]*a[5*sna_+5] - a[4*sna_+5]*a[5*sna_+4]
          };
          // 3x3 determinants of rows 4-5-6
          T_ w[20] =
          {
            a[3*sna_+0]*v[ 5] - a[3*sna_+1]*v[ 1] + a[3*sna_+2]*v[ 0],
            a[3*sna_+0]*v[ 6] - a[3*sna_+1]*v[ 2] + a[3*sna_+3]*v[ 0],
            a[3*sna_+0]*v[ 7] - a[3*sna_+1]*v[ 3] + a[3*sna_+4]*v[ 0],
            a[3*sna_+0]*v[ 8] - a[3*sna_+1]*v[ 4] + a[3*sna_+5]*v[ 0],
            a[3*sna_+0]*v[ 9] - a[3*sna_+2]*v[ 2] + a[3*sna_+3]*v[ 1],
            a[3*sna_+0]*v[10] - a[3*sna_+2]*v[ 3] + a[3*sna_+4]*v[ 1],
            a[3*sna_+0]*v[11] - a[3*sna_+2]*v[ 4] + a[3*sna_+5]*v[ 1],
            a[3*sna_+0]*v[12] - a[3*sna_+3]*v[ 3] + a[3*sna_+4]*v[ 2],
            a[3*sna_+0]*v[13] - a[3*sna_+3]*v[ 4] + a[3*sna_+5]*v[ 2],
            a[3*sna_+0]*v[14] - a[3*sna_+4]*v[ 4] + a[3*sna_+5]*v[ 3],
            a[3*sna_+1]*v[ 9] - a[3*sna_+2]*v[ 6] + a[3*sna_+3]*v[ 5],
            a[3*sna_+1]*v[10] - a[3*sna_+2]*v[ 7] + a[3*sna_+4]*v[ 5],
            a[3*sna_+1]*v[11] - a[3*sna_+2]*v[ 8] + a[3*sna_+5]*v[ 5],
            a[3*sna_+1]*v[12] - a[3*sna_+3]*v[ 7] + a[3*sna_+4]*v[ 6],
            a[3*sna_+1]*v[13] - a[3*sna_+3]*v[ 8] + a[3*sna_+5]*v[ 6],
            a[3*sna_+1]*v[14] - a[3*sna_+4]*v[ 8] + a[3*sna_+5]*v[ 7],
            a[3*sna_+2]*v[12] - a[3*sna_+3]*v[10] + a[3*sna_+4]*v[ 9],
            a[3*sna_+2]*v[13] - a[3*sna_+3]*v[11] + a[3*sna_+5]*v[ 9],
            a[3*sna_+2]*v[14] - a[3*sna_+4]*v[11] + a[3*sna_+5]*v[10],
            a[3*sna_+3]*v[14] - a[3*sna_+4]*v[13] + a[3*sna_+5]*v[12]
          };
          // 4x4 determinants of rows 3-4-5-6
          v[ 0] = a[2*sna_+0]*w[10] - a[2*sna_+1]*w[ 4] + a[2*sna_+2]*w[ 1] - a[2*sna_+3]*w[ 0];
          v[ 1] = a[2*sna_+0]*w[11] - a[2*sna_+1]*w[ 5] + a[2*sna_+2]*w[ 2] - a[2*sna_+4]*w[ 0];
          v[ 2] = a[2*sna_+0]*w[12] - a[2*sna_+1]*w[ 6] + a[2*sna_+2]*w[ 3] - a[2*sna_+5]*w[ 0];
          v[ 3] = a[2*sna_+0]*w[13] - a[2*sna_+1]*w[ 7] + a[2*sna_+3]*w[ 2] - a[2*sna_+4]*w[ 1];
          v[ 4] = a[2*sna_+0]*w[14] - a[2*sna_+1]*w[ 8] + a[2*sna_+3]*w[ 3] - a[2*sna_+5]*w[ 1];
          v[ 5] = a[2*sna_+0]*w[15] - a[2*sna_+1]*w[ 9] + a[2*sna_+4]*w[ 3] - a[2*sna_+5]*w[ 2];
          v[ 6] = a[2*sna_+0]*w[16] - a[2*sna_+2]*w[ 7] + a[2*sna_+3]*w[ 5] - a[2*sna_+4]*w[ 4];
          v[ 7] = a[2*sna_+0]*w[17] - a[2*sna_+2]*w[ 8] + a[2*sna_+3]*w[ 6] - a[2*sna_+5]*w[ 4];
          v[ 8] = a[2*sna_+0]*w[18] - a[2*sna_+2]*w[ 9] + a[2*sna_+4]*w[ 6] - a[2*sna_+5]*w[ 5];
          v[ 9] = a[2*sna_+0]*w[19] - a[2*sna_+3]*w[ 9] + a[2*sna_+4]*w[ 8] - a[2*sna_+5]*w[ 7];
          v[10] = a[2*sna_+1]*w[16] - a[2*sna_+2]*w[13] + a[2*sna_+3]*w[11] - a[2*sna_+4]*w[10];
          v[11] = a[2*sna_+1]*w[17] - a[2*sna_+2]*w[14] + a[2*sna_+3]*w[12] - a[2*sna_+5]*w[10];
          v[12] = a[2*sna_+1]*w[18] - a[2*sna_+2]*w[15] + a[2*sna_+4]*w[12] - a[2*sna_+5]*w[11];
          v[13] = a[2*sna_+1]*w[19] - a[2*sna_+3]*w[15] + a[2*sna_+4]*w[14] - a[2*sna_+5]*w[13];
          v[14] = a[2*sna_+2]*w[19] - a[2*sna_+3]*w[18] + a[2*sna_+4]*w[17] - a[2*sna_+5]*w[16];

          return
            + a[0*sna_+0]*(a[1*sna_+1]*v[14] - a[1*sna_+2]*v[13] + a[1*sna_+3]*v[12] - a[1*sna_+4]*v[11] + a[1*sna_+5]*v[10])
            - a[0*sna_+1]*(a[1*sna_+0]*v[14] - a[1*sna_+2]*v[ 9] + a[1*sna_+3]*v[ 8] - a[1*sna_+4]*v[ 7] + a[1*sna_+5]*v[ 6])
            + a[0*sna_+2]*(a[1*sna_+0]*v[13] - a[1*sna_+1]*v[ 9] + a[1*sna_+3]*v[ 5] - a[1*sna_+4]*v[ 4] + a[1*sna_+5]*v[ 3])
            - a[0*sna_+3]*(a[1*sna_+0]*v[12] - a[1*sna_+1]*v[ 8] + a[1*sna_+2]*v[ 5] - a[1*sna_+4]*v[ 2] + a[1*sna_+5]*v[ 1])
            + a[0*sna_+4]*(a[1*sna_+0]*v[11] - a[1*sna_+1]*v[ 7] + a[1*sna_+2]*v[ 4] - a[1*sna_+3]*v[ 2] + a[1*sna_+5]*v[ 0])
            - a[0*sna_+5]*(a[1*sna_+0]*v[10] - a[1*sna_+1]*v[ 6] + a[1*sna_+2]*v[ 3] - a[1*sna_+3]*v[ 1] + a[1*sna_+4]*v[ 0]);
        }
      };

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<int m_, int n_, int sna_>
      struct VolHelper
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // generic fallback implementation: compute b := a^T * a and return sqrt(det(b))
          T_ b[n_*n_];
          for(int i(0); i < n_; ++i)
          {
            for(int j(0); j < n_; ++j)
            {
              b[i*n_+j] = T_(0);
              for(int k(0); k < m_; ++k)
              {
                b[i*n_+j] += a[k*n_+i]*a[k*n_+j];
              }
            }
          }
          #ifndef __CUDACC__
          return Math::sqrt(DetHelper<n_,n_,n_>::compute(b));
          #else
          return CudaMath::cuda_sqrt(DetHelper<n_,n_,n_>::compute(b));
          #endif
        }
      };

      template<int n_, int sna_>
      struct VolHelper<n_,n_,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // square matrix special case: vol(a) = abs(det(a))
          #ifndef __CUDACC__
          return Math::abs(DetHelper<n_,n_,n_>::compute(a));
          #else
          return CudaMath::cuda_abs(DetHelper<n_,n_,n_>::compute(a));
          #endif
        }
      };

      template<int sna_>
      struct VolHelper<2,1,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // This is the euclid norm of the only matrix column.
          #ifndef __CUDACC__
          return Math::sqrt(Math::sqr(a[0*sna_]) + Math::sqr(a[1*sna_]));
          #else
          return CudaMath::cuda_sqrt(CudaMath::cuda_sqr(a[0*sna_]) + CudaMath::cuda_sqr(a[1*sna_]));
          #endif
        }
      };

      template<int sna_>
      struct VolHelper<3,1,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // This is the euclid norm of the only matrix column.
          #ifndef __CUDACC__
          return Math::sqrt(Math::sqr(a[0*sna_]) + Math::sqr(a[1*sna_]) + Math::sqr(a[2*sna_]));
          #else
          return CudaMath::cuda_sqrt(CudaMath::cuda_sqr(a[0*sna_]) + CudaMath::cuda_sqr(a[1*sna_]) + CudaMath::cuda_sqr(a[2*sna_]));
          #endif
        }
      };

      template<int sna_>
      struct VolHelper<3,2,sna_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(const T_* a)
        {
          // This is the euclid norm of the 3D cross product of the two matrix columns.
          #ifndef __CUDACC__
          return Math::sqrt(
            Math::sqr(a[1*sna_+0]*a[2*sna_+1] - a[2*sna_+0]*a[1*sna_+1]) +
            Math::sqr(a[2*sna_+0]*a[0*sna_+1] - a[0*sna_+0]*a[2*sna_+1]) +
            Math::sqr(a[0*sna_+0]*a[1*sna_+1] - a[1*sna_+0]*a[0*sna_+1]));
          #else
          return CudaMath::cuda_sqrt(
            CudaMath::cuda_sqr(a[1*sna_+0]*a[2*sna_+1] - a[2*sna_+0]*a[1*sna_+1]) +
            CudaMath::cuda_sqr(a[2*sna_+0]*a[0*sna_+1] - a[0*sna_+0]*a[2*sna_+1]) +
            CudaMath::cuda_sqr(a[0*sna_+0]*a[1*sna_+1] - a[1*sna_+0]*a[0*sna_+1]));
          #endif
        }
      };

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<int sna_, int snb_>
      struct InverseHelper<1,1,sna_,snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(T_* b, const T_* a)
        {
          b[0] = T_(1) / a[0];
          return a[0];
        }
      };

      template<int sna_, int snb_>
      struct InverseHelper<2,2,sna_,snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(T_* b, const T_* a)
        {
          T_ det = a[0*sna_+0]*a[1*sna_+1] - a[0*sna_+1]*a[1*sna_+0];
          T_ d = T_(1) / det;
          b[0*snb_+0] =  d*a[1*sna_+1];
          b[0*snb_+1] = -d*a[0*sna_+1];
          b[1*snb_+0] = -d*a[1*sna_+0];
          b[1*snb_+1] =  d*a[0*sna_+0];
          return det;
        }

        #ifdef __CUDACC__
        template<typename T_, typename ThreadGroup_>
        CUDA_DEVICE static __forceinline__ void grouped_compute(const ThreadGroup_& tg, T_* b, const T_* a, const T_& det)
        {
          //i think an if else cascade could do better than heavy modulo operations...
          T_ d = T_(1) / det;

          for(int idx = tg.thread_rank(); idx < 4; idx += tg.num_threads())
          {
            const int i = idx /2;
            const int j = idx %2;
            b[i*snb_+0] = ((i+j)==1 ? (-1) : (1)) * d * a[(1-j)*sna_+(1-i)];

          }
          b[0*snb_+0] =  d*a[1*sna_+1];
          b[0*snb_+1] = -d*a[0*sna_+1];
          b[1*snb_+0] = -d*a[1*sna_+0];
          b[1*snb_+1] =  d*a[0*sna_+0];
        }
        #endif
      };

      template<int sna_, int snb_>
      struct InverseHelper<3,3,sna_,snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(T_* b, const T_* a)
        {
          b[0*snb_+0] = a[1*sna_+1]*a[2*sna_+2] - a[1*sna_+2]*a[2*sna_+1];
          b[1*snb_+0] = a[1*sna_+2]*a[2*sna_+0] - a[1*sna_+0]*a[2*sna_+2];
          b[2*snb_+0] = a[1*sna_+0]*a[2*sna_+1] - a[1*sna_+1]*a[2*sna_+0];
          T_ det = a[0*sna_+0]*b[0*snb_+0] + a[0*sna_+1]*b[1*snb_+0] + a[0*sna_+2]*b[2*snb_+0];
          T_ d = T_(1) / det;
          b[0*snb_+0] *= d;
          b[1*snb_+0] *= d;
          b[2*snb_+0] *= d;
          b[0*snb_+1] = d*(a[0*sna_+2]*a[2*sna_+1] - a[0*sna_+1]*a[2*sna_+2]);
          b[1*snb_+1] = d*(a[0*sna_+0]*a[2*sna_+2] - a[0*sna_+2]*a[2*sna_+0]);
          b[2*snb_+1] = d*(a[0*sna_+1]*a[2*sna_+0] - a[0*sna_+0]*a[2*sna_+1]);
          b[0*snb_+2] = d*(a[0*sna_+1]*a[1*sna_+2] - a[0*sna_+2]*a[1*sna_+1]);
          b[1*snb_+2] = d*(a[0*sna_+2]*a[1*sna_+0] - a[0*sna_+0]*a[1*sna_+2]);
          b[2*snb_+2] = d*(a[0*sna_+0]*a[1*sna_+1] - a[0*sna_+1]*a[1*sna_+0]);
          return det;
        }

        #ifdef __CUDACC__
        template<typename T_, typename ThreadGroup_>
        CUDA_DEVICE static __forceinline__ void grouped_compute(const ThreadGroup_& tg, T_* b, const T_* a, const T_& det)
        {
          // for(int i = tg.thread_rank(); i < 3; i += tg.num_threads())
          // {
          //   b[i*snb_+0] = a[1*sna_+1+i-3*(i/2)]*a[2*sna_+(2*(1-i+i/2))] - a[1*sna_+(2*(1-i+i/2))]*a[2*sna_+1+i-3*(i/2)];
          // }
          // tg.sync();
          // T_ det = a[0*sna_+0]*b[0*snb_+0] + a[0*sna_+1]*b[1*snb_+0] + a[0*sna_+2]*b[2*snb_+0];
          // T_ d = T_(1) / det;
          // for(int idx = tg.thread_rank(); idx < 3; idx += tg.num_threads())
          // {
          //   const int i = idx/3;
          //   const int j = idx%3;
          //   b[i*snb_+j] = d*(a[()*sna_+1]*a[2*sna_+2] - a[1*sna_+2]*a[2*sna_+1])

          // }
          //i think an if else cascade could do better than heavy modulo operations...
          T_ d = T_(1) / det;

          for(int idx = tg.thread_rank(); idx < 9; idx += tg.num_threads())
          {
            if(idx == 0)
              b[0*snb_+0] = d*(a[1*sna_+1]*a[2*sna_+2] - a[1*sna_+2]*a[2*sna_+1]);
            else if(idx == 3)
              b[1*snb_+0] = d*(a[1*sna_+2]*a[2*sna_+0] - a[1*sna_+0]*a[2*sna_+2]);
            else if(idx == 6)
              b[2*snb_+0] = d*(a[1*sna_+0]*a[2*sna_+1] - a[1*sna_+1]*a[2*sna_+0]);
            else if(idx == 1)
              b[0*snb_+1] = d*(a[0*sna_+2]*a[2*sna_+1] - a[0*sna_+1]*a[2*sna_+2]);
            else if(idx == 4)
              b[1*snb_+1] = d*(a[0*sna_+0]*a[2*sna_+2] - a[0*sna_+2]*a[2*sna_+0]);
            else if(idx == 7)
              b[2*snb_+1] = d*(a[0*sna_+1]*a[2*sna_+0] - a[0*sna_+0]*a[2*sna_+1]);
            else if(idx == 2)
              b[0*snb_+2] = d*(a[0*sna_+1]*a[1*sna_+2] - a[0*sna_+2]*a[1*sna_+1]);
            else if(idx == 5)
              b[1*snb_+2] = d*(a[0*sna_+2]*a[1*sna_+0] - a[0*sna_+0]*a[1*sna_+2]);
            else if(idx == 8)
              b[2*snb_+2] = d*(a[0*sna_+0]*a[1*sna_+1] - a[0*sna_+1]*a[1*sna_+0]);
          }
        }
        #endif
      };

      template<int sna_, int snb_>
      struct InverseHelper<4,4,sna_,snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(T_* b, const T_* a)
        {
          T_ w[6];
          w[0] = a[2*sna_+0]*a[3*sna_+1]-a[2*sna_+1]*a[3*sna_+0];
          w[1] = a[2*sna_+0]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+0];
          w[2] = a[2*sna_+0]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+0];
          w[3] = a[2*sna_+1]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+1];
          w[4] = a[2*sna_+1]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+1];
          w[5] = a[2*sna_+2]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+2];
          b[0*snb_+0] = a[1*sna_+1]*w[5]-a[1*sna_+2]*w[4]+a[1*sna_+3]*w[3];
          b[1*snb_+0] =-a[1*sna_+0]*w[5]+a[1*sna_+2]*w[2]-a[1*sna_+3]*w[1];
          b[2*snb_+0] = a[1*sna_+0]*w[4]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          b[3*snb_+0] =-a[1*sna_+0]*w[3]+a[1*sna_+1]*w[1]-a[1*sna_+2]*w[0];
          T_ det = a[0*sna_+0]*b[0*snb_+0]+a[0*sna_+1]*b[1*snb_+0]+a[0*sna_+2]*b[2*snb_+0]+a[0*sna_+3]*b[3*snb_+0];
          T_ d = T_(1) / det;
          b[0*snb_+0] *= d;
          b[1*snb_+0] *= d;
          b[2*snb_+0] *= d;
          b[3*snb_+0] *= d;
          b[0*snb_+1] = d*(-a[0*sna_+1]*w[5]+a[0*sna_+2]*w[4]-a[0*sna_+3]*w[3]);
          b[1*snb_+1] = d*( a[0*sna_+0]*w[5]-a[0*sna_+2]*w[2]+a[0*sna_+3]*w[1]);
          b[2*snb_+1] = d*(-a[0*sna_+0]*w[4]+a[0*sna_+1]*w[2]-a[0*sna_+3]*w[0]);
          b[3*snb_+1] = d*( a[0*sna_+0]*w[3]-a[0*sna_+1]*w[1]+a[0*sna_+2]*w[0]);
          w[0] = a[0*sna_+0]*a[1*sna_+1]-a[0*sna_+1]*a[1*sna_+0];
          w[1] = a[0*sna_+0]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+0];
          w[2] = a[0*sna_+0]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+0];
          w[3] = a[0*sna_+1]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+1];
          w[4] = a[0*sna_+1]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+1];
          w[5] = a[0*sna_+2]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+2];
          b[0*snb_+2] = d*( a[3*sna_+1]*w[5]-a[3*sna_+2]*w[4]+a[3*sna_+3]*w[3]);
          b[1*snb_+2] = d*(-a[3*sna_+0]*w[5]+a[3*sna_+2]*w[2]-a[3*sna_+3]*w[1]);
          b[2*snb_+2] = d*( a[3*sna_+0]*w[4]-a[3*sna_+1]*w[2]+a[3*sna_+3]*w[0]);
          b[3*snb_+2] = d*(-a[3*sna_+0]*w[3]+a[3*sna_+1]*w[1]-a[3*sna_+2]*w[0]);
          b[0*snb_+3] = d*(-a[2*sna_+1]*w[5]+a[2*sna_+2]*w[4]-a[2*sna_+3]*w[3]);
          b[1*snb_+3] = d*( a[2*sna_+0]*w[5]-a[2*sna_+2]*w[2]+a[2*sna_+3]*w[1]);
          b[2*snb_+3] = d*(-a[2*sna_+0]*w[4]+a[2*sna_+1]*w[2]-a[2*sna_+3]*w[0]);
          b[3*snb_+3] = d*( a[2*sna_+0]*w[3]-a[2*sna_+1]*w[1]+a[2*sna_+2]*w[0]);
          return det;
        }

        #ifdef __CUDACC__
        template<typename T_, typename ThreadGroup_>
        CUDA_DEVICE static __forceinline__ void grouped_compute(const ThreadGroup_& tg, T_* b, const T_* a, const T_&)
        {
          constexpr int n_ = 4;
          // copy matrix a to b
          for (int i = tg.thread_rank(); i < n_*n_; i += tg.num_threads())
          {
            const int row = i/n_;
            const int col = i%n_;
            b[row*snb_+col] = a[row*sna_+col];
          }

          // create shared pivot array
          __shared__ int p[n_];
          for (int i = tg.thread_rank(); i < n_; i += tg.num_threads())
          {
            p[i] = 0;
          }
          tg.sync();
          // call grouped invert from first warp subgroup
          if(tg.thread_rank() < 32)
          {
            auto a_g = cg::coalesced_threads();

            CudaMath::cuda_grouped_invert_matrix(a_g, n_, snb_, b, p);
          }
          tg.sync();

        }
        #endif
      };

      template<int sna_, int snb_>
      struct InverseHelper<5,5,sna_,snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(T_* b, const T_* a)
        {
          T_ w[20];
          w[ 0] = a[3*sna_+0]*a[4*sna_+1]-a[3*sna_+1]*a[4*sna_+0];
          w[ 1] = a[3*sna_+0]*a[4*sna_+2]-a[3*sna_+2]*a[4*sna_+0];
          w[ 2] = a[3*sna_+0]*a[4*sna_+3]-a[3*sna_+3]*a[4*sna_+0];
          w[ 3] = a[3*sna_+0]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+0];
          w[ 4] = a[3*sna_+1]*a[4*sna_+2]-a[3*sna_+2]*a[4*sna_+1];
          w[ 5] = a[3*sna_+1]*a[4*sna_+3]-a[3*sna_+3]*a[4*sna_+1];
          w[ 6] = a[3*sna_+1]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+1];
          w[ 7] = a[3*sna_+2]*a[4*sna_+3]-a[3*sna_+3]*a[4*sna_+2];
          w[ 8] = a[3*sna_+2]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+2];
          w[ 9] = a[3*sna_+3]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+3];
          w[10] = a[2*sna_+0]*w[4]-a[2*sna_+1]*w[1]+a[2*sna_+2]*w[0];
          w[11] = a[2*sna_+0]*w[5]-a[2*sna_+1]*w[2]+a[2*sna_+3]*w[0];
          w[12] = a[2*sna_+0]*w[6]-a[2*sna_+1]*w[3]+a[2*sna_+4]*w[0];
          w[13] = a[2*sna_+0]*w[7]-a[2*sna_+2]*w[2]+a[2*sna_+3]*w[1];
          w[14] = a[2*sna_+0]*w[8]-a[2*sna_+2]*w[3]+a[2*sna_+4]*w[1];
          w[15] = a[2*sna_+0]*w[9]-a[2*sna_+3]*w[3]+a[2*sna_+4]*w[2];
          w[16] = a[2*sna_+1]*w[7]-a[2*sna_+2]*w[5]+a[2*sna_+3]*w[4];
          w[17] = a[2*sna_+1]*w[8]-a[2*sna_+2]*w[6]+a[2*sna_+4]*w[4];
          w[18] = a[2*sna_+1]*w[9]-a[2*sna_+3]*w[6]+a[2*sna_+4]*w[5];
          w[19] = a[2*sna_+2]*w[9]-a[2*sna_+3]*w[8]+a[2*sna_+4]*w[7];
          b[0*snb_+0] = a[1*sna_+1]*w[19]-a[1*sna_+2]*w[18]+a[1*sna_+3]*w[17]-a[1*sna_+4]*w[16];
          b[1*snb_+0] =-a[1*sna_+0]*w[19]+a[1*sna_+2]*w[15]-a[1*sna_+3]*w[14]+a[1*sna_+4]*w[13];
          b[2*snb_+0] = a[1*sna_+0]*w[18]-a[1*sna_+1]*w[15]+a[1*sna_+3]*w[12]-a[1*sna_+4]*w[11];
          b[3*snb_+0] =-a[1*sna_+0]*w[17]+a[1*sna_+1]*w[14]-a[1*sna_+2]*w[12]+a[1*sna_+4]*w[10];
          b[4*snb_+0] = a[1*sna_+0]*w[16]-a[1*sna_+1]*w[13]+a[1*sna_+2]*w[11]-a[1*sna_+3]*w[10];
          T_ det = a[0*sna_+0]*b[0*snb_+0]+a[0*sna_+1]*b[1*snb_+0]+a[0*sna_+2]*b[2*snb_+0]+a[0*sna_+3]*b[3*snb_+0]+a[0*sna_+4]*b[4*snb_+0];
          T_ d = T_(1) / det;
          b[0*snb_+0] *= d;
          b[1*snb_+0] *= d;
          b[2*snb_+0] *= d;
          b[3*snb_+0] *= d;
          b[4*snb_+0] *= d;
          b[0*snb_+1] = d*(-a[0*sna_+1]*w[19]+a[0*sna_+2]*w[18]-a[0*sna_+3]*w[17]+a[0*sna_+4]*w[16]);
          b[1*snb_+1] = d*( a[0*sna_+0]*w[19]-a[0*sna_+2]*w[15]+a[0*sna_+3]*w[14]-a[0*sna_+4]*w[13]);
          b[2*snb_+1] = d*(-a[0*sna_+0]*w[18]+a[0*sna_+1]*w[15]-a[0*sna_+3]*w[12]+a[0*sna_+4]*w[11]);
          b[3*snb_+1] = d*( a[0*sna_+0]*w[17]-a[0*sna_+1]*w[14]+a[0*sna_+2]*w[12]-a[0*sna_+4]*w[10]);
          b[4*snb_+1] = d*(-a[0*sna_+0]*w[16]+a[0*sna_+1]*w[13]-a[0*sna_+2]*w[11]+a[0*sna_+3]*w[10]);
          w[10] = a[1*sna_+0]*w[4]-a[1*sna_+1]*w[1]+a[1*sna_+2]*w[0];
          w[11] = a[1*sna_+0]*w[5]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          w[12] = a[1*sna_+0]*w[6]-a[1*sna_+1]*w[3]+a[1*sna_+4]*w[0];
          w[13] = a[1*sna_+0]*w[7]-a[1*sna_+2]*w[2]+a[1*sna_+3]*w[1];
          w[14] = a[1*sna_+0]*w[8]-a[1*sna_+2]*w[3]+a[1*sna_+4]*w[1];
          w[15] = a[1*sna_+0]*w[9]-a[1*sna_+3]*w[3]+a[1*sna_+4]*w[2];
          w[16] = a[1*sna_+1]*w[7]-a[1*sna_+2]*w[5]+a[1*sna_+3]*w[4];
          w[17] = a[1*sna_+1]*w[8]-a[1*sna_+2]*w[6]+a[1*sna_+4]*w[4];
          w[18] = a[1*sna_+1]*w[9]-a[1*sna_+3]*w[6]+a[1*sna_+4]*w[5];
          w[19] = a[1*sna_+2]*w[9]-a[1*sna_+3]*w[8]+a[1*sna_+4]*w[7];
          b[0*snb_+2] = d*( a[0*sna_+1]*w[19]-a[0*sna_+2]*w[18]+a[0*sna_+3]*w[17]-a[0*sna_+4]*w[16]);
          b[1*snb_+2] = d*(-a[0*sna_+0]*w[19]+a[0*sna_+2]*w[15]-a[0*sna_+3]*w[14]+a[0*sna_+4]*w[13]);
          b[2*snb_+2] = d*( a[0*sna_+0]*w[18]-a[0*sna_+1]*w[15]+a[0*sna_+3]*w[12]-a[0*sna_+4]*w[11]);
          b[3*snb_+2] = d*(-a[0*sna_+0]*w[17]+a[0*sna_+1]*w[14]-a[0*sna_+2]*w[12]+a[0*sna_+4]*w[10]);
          b[4*snb_+2] = d*( a[0*sna_+0]*w[16]-a[0*sna_+1]*w[13]+a[0*sna_+2]*w[11]-a[0*sna_+3]*w[10]);
          w[ 0] = a[0*sna_+0]*a[1*sna_+1]-a[0*sna_+1]*a[1*sna_+0];
          w[ 1] = a[0*sna_+0]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+0];
          w[ 2] = a[0*sna_+0]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+0];
          w[ 3] = a[0*sna_+0]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+0];
          w[ 4] = a[0*sna_+1]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+1];
          w[ 5] = a[0*sna_+1]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+1];
          w[ 6] = a[0*sna_+1]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+1];
          w[ 7] = a[0*sna_+2]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+2];
          w[ 8] = a[0*sna_+2]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+2];
          w[ 9] = a[0*sna_+3]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+3];
          w[10] = a[2*sna_+0]*w[4]-a[2*sna_+1]*w[1]+a[2*sna_+2]*w[0];
          w[11] = a[2*sna_+0]*w[5]-a[2*sna_+1]*w[2]+a[2*sna_+3]*w[0];
          w[12] = a[2*sna_+0]*w[6]-a[2*sna_+1]*w[3]+a[2*sna_+4]*w[0];
          w[13] = a[2*sna_+0]*w[7]-a[2*sna_+2]*w[2]+a[2*sna_+3]*w[1];
          w[14] = a[2*sna_+0]*w[8]-a[2*sna_+2]*w[3]+a[2*sna_+4]*w[1];
          w[15] = a[2*sna_+0]*w[9]-a[2*sna_+3]*w[3]+a[2*sna_+4]*w[2];
          w[16] = a[2*sna_+1]*w[7]-a[2*sna_+2]*w[5]+a[2*sna_+3]*w[4];
          w[17] = a[2*sna_+1]*w[8]-a[2*sna_+2]*w[6]+a[2*sna_+4]*w[4];
          w[18] = a[2*sna_+1]*w[9]-a[2*sna_+3]*w[6]+a[2*sna_+4]*w[5];
          w[19] = a[2*sna_+2]*w[9]-a[2*sna_+3]*w[8]+a[2*sna_+4]*w[7];
          b[0*snb_+3] = d*( a[4*sna_+1]*w[19]-a[4*sna_+2]*w[18]+a[4*sna_+3]*w[17]-a[4*sna_+4]*w[16]);
          b[1*snb_+3] = d*(-a[4*sna_+0]*w[19]+a[4*sna_+2]*w[15]-a[4*sna_+3]*w[14]+a[4*sna_+4]*w[13]);
          b[2*snb_+3] = d*( a[4*sna_+0]*w[18]-a[4*sna_+1]*w[15]+a[4*sna_+3]*w[12]-a[4*sna_+4]*w[11]);
          b[3*snb_+3] = d*(-a[4*sna_+0]*w[17]+a[4*sna_+1]*w[14]-a[4*sna_+2]*w[12]+a[4*sna_+4]*w[10]);
          b[4*snb_+3] = d*( a[4*sna_+0]*w[16]-a[4*sna_+1]*w[13]+a[4*sna_+2]*w[11]-a[4*sna_+3]*w[10]);
          b[0*snb_+4] = d*(-a[3*sna_+1]*w[19]+a[3*sna_+2]*w[18]-a[3*sna_+3]*w[17]+a[3*sna_+4]*w[16]);
          b[1*snb_+4] = d*( a[3*sna_+0]*w[19]-a[3*sna_+2]*w[15]+a[3*sna_+3]*w[14]-a[3*sna_+4]*w[13]);
          b[2*snb_+4] = d*(-a[3*sna_+0]*w[18]+a[3*sna_+1]*w[15]-a[3*sna_+3]*w[12]+a[3*sna_+4]*w[11]);
          b[3*snb_+4] = d*( a[3*sna_+0]*w[17]-a[3*sna_+1]*w[14]+a[3*sna_+2]*w[12]-a[3*sna_+4]*w[10]);
          b[4*snb_+4] = d*(-a[3*sna_+0]*w[16]+a[3*sna_+1]*w[13]-a[3*sna_+2]*w[11]+a[3*sna_+3]*w[10]);
          return det;
        }

        #ifdef __CUDACC__
        template<typename T_, typename ThreadGroup_>
        CUDA_DEVICE static __forceinline__ void grouped_compute(const ThreadGroup_& tg, T_* b, const T_* a, const T_&)
        {
          constexpr int n_ = 5;
          // copy matrix a to b
          for (int i = tg.thread_rank(); i < n_*n_; i += tg.num_threads())
          {
            const int row = i/n_;
            const int col = i%n_;
            b[row*snb_+col] = a[row*sna_+col];
          }

          // create shared pivot array
          __shared__ int p[n_];
          for (int i = tg.thread_rank(); i < n_; i += tg.num_threads())
          {
            p[i] = 0;
          }
          tg.sync();
          // call grouped invert from first warp subgroup
          if(tg.thread_rank() < 32)
          {
            auto a_g = cg::coalesced_threads();

            CudaMath::cuda_grouped_invert_matrix(a_g, n_, snb_, b, p);
          }
          tg.sync();

        }
        #endif
      };

      template<int sna_, int snb_>
      struct InverseHelper<6,6,sna_,snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(T_* b, const T_* a)
        {
          T_ w[35];
          w[ 0] = a[4*sna_+0]*a[5*sna_+1]-a[4*sna_+1]*a[5*sna_+0];
          w[ 1] = a[4*sna_+0]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+0];
          w[ 2] = a[4*sna_+0]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+0];
          w[ 3] = a[4*sna_+0]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+0];
          w[ 4] = a[4*sna_+0]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+0];
          w[ 5] = a[4*sna_+1]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+1];
          w[ 6] = a[4*sna_+1]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+1];
          w[ 7] = a[4*sna_+1]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+1];
          w[ 8] = a[4*sna_+1]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+1];
          w[ 9] = a[4*sna_+2]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+2];
          w[10] = a[4*sna_+2]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+2];
          w[11] = a[4*sna_+2]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+2];
          w[12] = a[4*sna_+3]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+3];
          w[13] = a[4*sna_+3]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+3];
          w[14] = a[4*sna_+4]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+4];
          w[15] = a[3*sna_+0]*w[5]-a[3*sna_+1]*w[1]+a[3*sna_+2]*w[0];
          w[16] = a[3*sna_+0]*w[6]-a[3*sna_+1]*w[2]+a[3*sna_+3]*w[0];
          w[17] = a[3*sna_+0]*w[7]-a[3*sna_+1]*w[3]+a[3*sna_+4]*w[0];
          w[18] = a[3*sna_+0]*w[8]-a[3*sna_+1]*w[4]+a[3*sna_+5]*w[0];
          w[19] = a[3*sna_+0]*w[9]-a[3*sna_+2]*w[2]+a[3*sna_+3]*w[1];
          w[20] = a[3*sna_+0]*w[10]-a[3*sna_+2]*w[3]+a[3*sna_+4]*w[1];
          w[21] = a[3*sna_+0]*w[11]-a[3*sna_+2]*w[4]+a[3*sna_+5]*w[1];
          w[22] = a[3*sna_+0]*w[12]-a[3*sna_+3]*w[3]+a[3*sna_+4]*w[2];
          w[23] = a[3*sna_+0]*w[13]-a[3*sna_+3]*w[4]+a[3*sna_+5]*w[2];
          w[24] = a[3*sna_+0]*w[14]-a[3*sna_+4]*w[4]+a[3*sna_+5]*w[3];
          w[25] = a[3*sna_+1]*w[9]-a[3*sna_+2]*w[6]+a[3*sna_+3]*w[5];
          w[26] = a[3*sna_+1]*w[10]-a[3*sna_+2]*w[7]+a[3*sna_+4]*w[5];
          w[27] = a[3*sna_+1]*w[11]-a[3*sna_+2]*w[8]+a[3*sna_+5]*w[5];
          w[28] = a[3*sna_+1]*w[12]-a[3*sna_+3]*w[7]+a[3*sna_+4]*w[6];
          w[29] = a[3*sna_+1]*w[13]-a[3*sna_+3]*w[8]+a[3*sna_+5]*w[6];
          w[30] = a[3*sna_+1]*w[14]-a[3*sna_+4]*w[8]+a[3*sna_+5]*w[7];
          w[31] = a[3*sna_+2]*w[12]-a[3*sna_+3]*w[10]+a[3*sna_+4]*w[9];
          w[32] = a[3*sna_+2]*w[13]-a[3*sna_+3]*w[11]+a[3*sna_+5]*w[9];
          w[33] = a[3*sna_+2]*w[14]-a[3*sna_+4]*w[11]+a[3*sna_+5]*w[10];
          w[34] = a[3*sna_+3]*w[14]-a[3*sna_+4]*w[13]+a[3*sna_+5]*w[12];
          w[ 0] = a[2*sna_+0]*w[25]-a[2*sna_+1]*w[19]+a[2*sna_+2]*w[16]-a[2*sna_+3]*w[15];
          w[ 1] = a[2*sna_+0]*w[26]-a[2*sna_+1]*w[20]+a[2*sna_+2]*w[17]-a[2*sna_+4]*w[15];
          w[ 2] = a[2*sna_+0]*w[27]-a[2*sna_+1]*w[21]+a[2*sna_+2]*w[18]-a[2*sna_+5]*w[15];
          w[ 3] = a[2*sna_+0]*w[28]-a[2*sna_+1]*w[22]+a[2*sna_+3]*w[17]-a[2*sna_+4]*w[16];
          w[ 4] = a[2*sna_+0]*w[29]-a[2*sna_+1]*w[23]+a[2*sna_+3]*w[18]-a[2*sna_+5]*w[16];
          w[ 5] = a[2*sna_+0]*w[30]-a[2*sna_+1]*w[24]+a[2*sna_+4]*w[18]-a[2*sna_+5]*w[17];
          w[ 6] = a[2*sna_+0]*w[31]-a[2*sna_+2]*w[22]+a[2*sna_+3]*w[20]-a[2*sna_+4]*w[19];
          w[ 7] = a[2*sna_+0]*w[32]-a[2*sna_+2]*w[23]+a[2*sna_+3]*w[21]-a[2*sna_+5]*w[19];
          w[ 8] = a[2*sna_+0]*w[33]-a[2*sna_+2]*w[24]+a[2*sna_+4]*w[21]-a[2*sna_+5]*w[20];
          w[ 9] = a[2*sna_+0]*w[34]-a[2*sna_+3]*w[24]+a[2*sna_+4]*w[23]-a[2*sna_+5]*w[22];
          w[10] = a[2*sna_+1]*w[31]-a[2*sna_+2]*w[28]+a[2*sna_+3]*w[26]-a[2*sna_+4]*w[25];
          w[11] = a[2*sna_+1]*w[32]-a[2*sna_+2]*w[29]+a[2*sna_+3]*w[27]-a[2*sna_+5]*w[25];
          w[12] = a[2*sna_+1]*w[33]-a[2*sna_+2]*w[30]+a[2*sna_+4]*w[27]-a[2*sna_+5]*w[26];
          w[13] = a[2*sna_+1]*w[34]-a[2*sna_+3]*w[30]+a[2*sna_+4]*w[29]-a[2*sna_+5]*w[28];
          w[14] = a[2*sna_+2]*w[34]-a[2*sna_+3]*w[33]+a[2*sna_+4]*w[32]-a[2*sna_+5]*w[31];
          b[0*snb_+0] =  a[1*sna_+1]*w[14]-a[1*sna_+2]*w[13]+a[1*sna_+3]*w[12]-a[1*sna_+4]*w[11]+a[1*sna_+5]*w[10];
          b[1*snb_+0] = -a[1*sna_+0]*w[14]+a[1*sna_+2]*w[9]-a[1*sna_+3]*w[8]+a[1*sna_+4]*w[7]-a[1*sna_+5]*w[6];
          b[2*snb_+0] =  a[1*sna_+0]*w[13]-a[1*sna_+1]*w[9]+a[1*sna_+3]*w[5]-a[1*sna_+4]*w[4]+a[1*sna_+5]*w[3];
          b[3*snb_+0] = -a[1*sna_+0]*w[12]+a[1*sna_+1]*w[8]-a[1*sna_+2]*w[5]+a[1*sna_+4]*w[2]-a[1*sna_+5]*w[1];
          b[4*snb_+0] =  a[1*sna_+0]*w[11]-a[1*sna_+1]*w[7]+a[1*sna_+2]*w[4]-a[1*sna_+3]*w[2]+a[1*sna_+5]*w[0];
          b[5*snb_+0] = -a[1*sna_+0]*w[10]+a[1*sna_+1]*w[6]-a[1*sna_+2]*w[3]+a[1*sna_+3]*w[1]-a[1*sna_+4]*w[0];
          T_ det = a[0*sna_+0]*b[0*snb_+0] + a[0*sna_+1]*b[1*snb_+0] + a[0*sna_+2]*b[2*snb_+0]
                 + a[0*sna_+3]*b[3*snb_+0] + a[0*sna_+4]*b[4*snb_+0] + a[0*sna_+5]*b[5*snb_+0];
          T_ d = T_(1) / det;
          b[0*snb_+0] *= d;
          b[1*snb_+0] *= d;
          b[2*snb_+0] *= d;
          b[3*snb_+0] *= d;
          b[4*snb_+0] *= d;
          b[5*snb_+0] *= d;
          b[0*snb_+1] = d*(-a[0*sna_+1]*w[14]+a[0*sna_+2]*w[13]-a[0*sna_+3]*w[12]+a[0*sna_+4]*w[11]-a[0*sna_+5]*w[10]);
          b[1*snb_+1] = d*( a[0*sna_+0]*w[14]-a[0*sna_+2]*w[9]+a[0*sna_+3]*w[8]-a[0*sna_+4]*w[7]+a[0*sna_+5]*w[6]);
          b[2*snb_+1] = d*(-a[0*sna_+0]*w[13]+a[0*sna_+1]*w[9]-a[0*sna_+3]*w[5]+a[0*sna_+4]*w[4]-a[0*sna_+5]*w[3]);
          b[3*snb_+1] = d*( a[0*sna_+0]*w[12]-a[0*sna_+1]*w[8]+a[0*sna_+2]*w[5]-a[0*sna_+4]*w[2]+a[0*sna_+5]*w[1]);
          b[4*snb_+1] = d*(-a[0*sna_+0]*w[11]+a[0*sna_+1]*w[7]-a[0*sna_+2]*w[4]+a[0*sna_+3]*w[2]-a[0*sna_+5]*w[0]);
          b[5*snb_+1] = d*( a[0*sna_+0]*w[10]-a[0*sna_+1]*w[6]+a[0*sna_+2]*w[3]-a[0*sna_+3]*w[1]+a[0*sna_+4]*w[0]);
          w[ 0] = a[4*sna_+0]*a[5*sna_+1]-a[4*sna_+1]*a[5*sna_+0];
          w[ 1] = a[4*sna_+0]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+0];
          w[ 2] = a[4*sna_+0]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+0];
          w[ 3] = a[4*sna_+0]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+0];
          w[ 4] = a[4*sna_+0]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+0];
          w[ 5] = a[4*sna_+1]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+1];
          w[ 6] = a[4*sna_+1]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+1];
          w[ 7] = a[4*sna_+1]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+1];
          w[ 8] = a[4*sna_+1]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+1];
          w[ 9] = a[4*sna_+2]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+2];
          w[10] = a[4*sna_+2]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+2];
          w[11] = a[4*sna_+2]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+2];
          w[12] = a[4*sna_+3]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+3];
          w[13] = a[4*sna_+3]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+3];
          w[14] = a[4*sna_+4]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+4];
          w[15] = a[1*sna_+0]*w[5]-a[1*sna_+1]*w[1]+a[1*sna_+2]*w[0];
          w[16] = a[1*sna_+0]*w[6]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          w[17] = a[1*sna_+0]*w[7]-a[1*sna_+1]*w[3]+a[1*sna_+4]*w[0];
          w[18] = a[1*sna_+0]*w[8]-a[1*sna_+1]*w[4]+a[1*sna_+5]*w[0];
          w[19] = a[1*sna_+0]*w[9]-a[1*sna_+2]*w[2]+a[1*sna_+3]*w[1];
          w[20] = a[1*sna_+0]*w[10]-a[1*sna_+2]*w[3]+a[1*sna_+4]*w[1];
          w[21] = a[1*sna_+0]*w[11]-a[1*sna_+2]*w[4]+a[1*sna_+5]*w[1];
          w[22] = a[1*sna_+0]*w[12]-a[1*sna_+3]*w[3]+a[1*sna_+4]*w[2];
          w[23] = a[1*sna_+0]*w[13]-a[1*sna_+3]*w[4]+a[1*sna_+5]*w[2];
          w[24] = a[1*sna_+0]*w[14]-a[1*sna_+4]*w[4]+a[1*sna_+5]*w[3];
          w[25] = a[1*sna_+1]*w[9]-a[1*sna_+2]*w[6]+a[1*sna_+3]*w[5];
          w[26] = a[1*sna_+1]*w[10]-a[1*sna_+2]*w[7]+a[1*sna_+4]*w[5];
          w[27] = a[1*sna_+1]*w[11]-a[1*sna_+2]*w[8]+a[1*sna_+5]*w[5];
          w[28] = a[1*sna_+1]*w[12]-a[1*sna_+3]*w[7]+a[1*sna_+4]*w[6];
          w[29] = a[1*sna_+1]*w[13]-a[1*sna_+3]*w[8]+a[1*sna_+5]*w[6];
          w[30] = a[1*sna_+1]*w[14]-a[1*sna_+4]*w[8]+a[1*sna_+5]*w[7];
          w[31] = a[1*sna_+2]*w[12]-a[1*sna_+3]*w[10]+a[1*sna_+4]*w[9];
          w[32] = a[1*sna_+2]*w[13]-a[1*sna_+3]*w[11]+a[1*sna_+5]*w[9];
          w[33] = a[1*sna_+2]*w[14]-a[1*sna_+4]*w[11]+a[1*sna_+5]*w[10];
          w[34] = a[1*sna_+3]*w[14]-a[1*sna_+4]*w[13]+a[1*sna_+5]*w[12];
          w[ 0] = a[0*sna_+0]*w[25]-a[0*sna_+1]*w[19]+a[0*sna_+2]*w[16]-a[0*sna_+3]*w[15];
          w[ 1] = a[0*sna_+0]*w[26]-a[0*sna_+1]*w[20]+a[0*sna_+2]*w[17]-a[0*sna_+4]*w[15];
          w[ 2] = a[0*sna_+0]*w[27]-a[0*sna_+1]*w[21]+a[0*sna_+2]*w[18]-a[0*sna_+5]*w[15];
          w[ 3] = a[0*sna_+0]*w[28]-a[0*sna_+1]*w[22]+a[0*sna_+3]*w[17]-a[0*sna_+4]*w[16];
          w[ 4] = a[0*sna_+0]*w[29]-a[0*sna_+1]*w[23]+a[0*sna_+3]*w[18]-a[0*sna_+5]*w[16];
          w[ 5] = a[0*sna_+0]*w[30]-a[0*sna_+1]*w[24]+a[0*sna_+4]*w[18]-a[0*sna_+5]*w[17];
          w[ 6] = a[0*sna_+0]*w[31]-a[0*sna_+2]*w[22]+a[0*sna_+3]*w[20]-a[0*sna_+4]*w[19];
          w[ 7] = a[0*sna_+0]*w[32]-a[0*sna_+2]*w[23]+a[0*sna_+3]*w[21]-a[0*sna_+5]*w[19];
          w[ 8] = a[0*sna_+0]*w[33]-a[0*sna_+2]*w[24]+a[0*sna_+4]*w[21]-a[0*sna_+5]*w[20];
          w[ 9] = a[0*sna_+0]*w[34]-a[0*sna_+3]*w[24]+a[0*sna_+4]*w[23]-a[0*sna_+5]*w[22];
          w[10] = a[0*sna_+1]*w[31]-a[0*sna_+2]*w[28]+a[0*sna_+3]*w[26]-a[0*sna_+4]*w[25];
          w[11] = a[0*sna_+1]*w[32]-a[0*sna_+2]*w[29]+a[0*sna_+3]*w[27]-a[0*sna_+5]*w[25];
          w[12] = a[0*sna_+1]*w[33]-a[0*sna_+2]*w[30]+a[0*sna_+4]*w[27]-a[0*sna_+5]*w[26];
          w[13] = a[0*sna_+1]*w[34]-a[0*sna_+3]*w[30]+a[0*sna_+4]*w[29]-a[0*sna_+5]*w[28];
          w[14] = a[0*sna_+2]*w[34]-a[0*sna_+3]*w[33]+a[0*sna_+4]*w[32]-a[0*sna_+5]*w[31];
          b[0*snb_+2] = d*( a[3*sna_+1]*w[14]-a[3*sna_+2]*w[13]+a[3*sna_+3]*w[12]-a[3*sna_+4]*w[11]+a[3*sna_+5]*w[10]);
          b[1*snb_+2] = d*(-a[3*sna_+0]*w[14]+a[3*sna_+2]*w[9]-a[3*sna_+3]*w[8]+a[3*sna_+4]*w[7]-a[3*sna_+5]*w[6]);
          b[2*snb_+2] = d*( a[3*sna_+0]*w[13]-a[3*sna_+1]*w[9]+a[3*sna_+3]*w[5]-a[3*sna_+4]*w[4]+a[3*sna_+5]*w[3]);
          b[3*snb_+2] = d*(-a[3*sna_+0]*w[12]+a[3*sna_+1]*w[8]-a[3*sna_+2]*w[5]+a[3*sna_+4]*w[2]-a[3*sna_+5]*w[1]);
          b[4*snb_+2] = d*( a[3*sna_+0]*w[11]-a[3*sna_+1]*w[7]+a[3*sna_+2]*w[4]-a[3*sna_+3]*w[2]+a[3*sna_+5]*w[0]);
          b[5*snb_+2] = d*(-a[3*sna_+0]*w[10]+a[3*sna_+1]*w[6]-a[3*sna_+2]*w[3]+a[3*sna_+3]*w[1]-a[3*sna_+4]*w[0]);
          b[0*snb_+3] = d*(-a[2*sna_+1]*w[14]+a[2*sna_+2]*w[13]-a[2*sna_+3]*w[12]+a[2*sna_+4]*w[11]-a[2*sna_+5]*w[10]);
          b[1*snb_+3] = d*( a[2*sna_+0]*w[14]-a[2*sna_+2]*w[9]+a[2*sna_+3]*w[8]-a[2*sna_+4]*w[7]+a[2*sna_+5]*w[6]);
          b[2*snb_+3] = d*(-a[2*sna_+0]*w[13]+a[2*sna_+1]*w[9]-a[2*sna_+3]*w[5]+a[2*sna_+4]*w[4]-a[2*sna_+5]*w[3]);
          b[3*snb_+3] = d*( a[2*sna_+0]*w[12]-a[2*sna_+1]*w[8]+a[2*sna_+2]*w[5]-a[2*sna_+4]*w[2]+a[2*sna_+5]*w[1]);
          b[4*snb_+3] = d*(-a[2*sna_+0]*w[11]+a[2*sna_+1]*w[7]-a[2*sna_+2]*w[4]+a[2*sna_+3]*w[2]-a[2*sna_+5]*w[0]);
          b[5*snb_+3] = d*( a[2*sna_+0]*w[10]-a[2*sna_+1]*w[6]+a[2*sna_+2]*w[3]-a[2*sna_+3]*w[1]+a[2*sna_+4]*w[0]);
          w[ 0] = a[2*sna_+0]*a[3*sna_+1]-a[2*sna_+1]*a[3*sna_+0];
          w[ 1] = a[2*sna_+0]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+0];
          w[ 2] = a[2*sna_+0]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+0];
          w[ 3] = a[2*sna_+0]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+0];
          w[ 4] = a[2*sna_+0]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+0];
          w[ 5] = a[2*sna_+1]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+1];
          w[ 6] = a[2*sna_+1]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+1];
          w[ 7] = a[2*sna_+1]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+1];
          w[ 8] = a[2*sna_+1]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+1];
          w[ 9] = a[2*sna_+2]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+2];
          w[10] = a[2*sna_+2]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+2];
          w[11] = a[2*sna_+2]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+2];
          w[12] = a[2*sna_+3]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+3];
          w[13] = a[2*sna_+3]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+3];
          w[14] = a[2*sna_+4]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+4];
          w[15] = a[1*sna_+0]*w[5]-a[1*sna_+1]*w[1]+a[1*sna_+2]*w[0];
          w[16] = a[1*sna_+0]*w[6]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          w[17] = a[1*sna_+0]*w[7]-a[1*sna_+1]*w[3]+a[1*sna_+4]*w[0];
          w[18] = a[1*sna_+0]*w[8]-a[1*sna_+1]*w[4]+a[1*sna_+5]*w[0];
          w[19] = a[1*sna_+0]*w[9]-a[1*sna_+2]*w[2]+a[1*sna_+3]*w[1];
          w[20] = a[1*sna_+0]*w[10]-a[1*sna_+2]*w[3]+a[1*sna_+4]*w[1];
          w[21] = a[1*sna_+0]*w[11]-a[1*sna_+2]*w[4]+a[1*sna_+5]*w[1];
          w[22] = a[1*sna_+0]*w[12]-a[1*sna_+3]*w[3]+a[1*sna_+4]*w[2];
          w[23] = a[1*sna_+0]*w[13]-a[1*sna_+3]*w[4]+a[1*sna_+5]*w[2];
          w[24] = a[1*sna_+0]*w[14]-a[1*sna_+4]*w[4]+a[1*sna_+5]*w[3];
          w[25] = a[1*sna_+1]*w[9]-a[1*sna_+2]*w[6]+a[1*sna_+3]*w[5];
          w[26] = a[1*sna_+1]*w[10]-a[1*sna_+2]*w[7]+a[1*sna_+4]*w[5];
          w[27] = a[1*sna_+1]*w[11]-a[1*sna_+2]*w[8]+a[1*sna_+5]*w[5];
          w[28] = a[1*sna_+1]*w[12]-a[1*sna_+3]*w[7]+a[1*sna_+4]*w[6];
          w[29] = a[1*sna_+1]*w[13]-a[1*sna_+3]*w[8]+a[1*sna_+5]*w[6];
          w[30] = a[1*sna_+1]*w[14]-a[1*sna_+4]*w[8]+a[1*sna_+5]*w[7];
          w[31] = a[1*sna_+2]*w[12]-a[1*sna_+3]*w[10]+a[1*sna_+4]*w[9];
          w[32] = a[1*sna_+2]*w[13]-a[1*sna_+3]*w[11]+a[1*sna_+5]*w[9];
          w[33] = a[1*sna_+2]*w[14]-a[1*sna_+4]*w[11]+a[1*sna_+5]*w[10];
          w[34] = a[1*sna_+3]*w[14]-a[1*sna_+4]*w[13]+a[1*sna_+5]*w[12];
          w[ 0] = a[0*sna_+0]*w[25]-a[0*sna_+1]*w[19]+a[0*sna_+2]*w[16]-a[0*sna_+3]*w[15];
          w[ 1] = a[0*sna_+0]*w[26]-a[0*sna_+1]*w[20]+a[0*sna_+2]*w[17]-a[0*sna_+4]*w[15];
          w[ 2] = a[0*sna_+0]*w[27]-a[0*sna_+1]*w[21]+a[0*sna_+2]*w[18]-a[0*sna_+5]*w[15];
          w[ 3] = a[0*sna_+0]*w[28]-a[0*sna_+1]*w[22]+a[0*sna_+3]*w[17]-a[0*sna_+4]*w[16];
          w[ 4] = a[0*sna_+0]*w[29]-a[0*sna_+1]*w[23]+a[0*sna_+3]*w[18]-a[0*sna_+5]*w[16];
          w[ 5] = a[0*sna_+0]*w[30]-a[0*sna_+1]*w[24]+a[0*sna_+4]*w[18]-a[0*sna_+5]*w[17];
          w[ 6] = a[0*sna_+0]*w[31]-a[0*sna_+2]*w[22]+a[0*sna_+3]*w[20]-a[0*sna_+4]*w[19];
          w[ 7] = a[0*sna_+0]*w[32]-a[0*sna_+2]*w[23]+a[0*sna_+3]*w[21]-a[0*sna_+5]*w[19];
          w[ 8] = a[0*sna_+0]*w[33]-a[0*sna_+2]*w[24]+a[0*sna_+4]*w[21]-a[0*sna_+5]*w[20];
          w[ 9] = a[0*sna_+0]*w[34]-a[0*sna_+3]*w[24]+a[0*sna_+4]*w[23]-a[0*sna_+5]*w[22];
          w[10] = a[0*sna_+1]*w[31]-a[0*sna_+2]*w[28]+a[0*sna_+3]*w[26]-a[0*sna_+4]*w[25];
          w[11] = a[0*sna_+1]*w[32]-a[0*sna_+2]*w[29]+a[0*sna_+3]*w[27]-a[0*sna_+5]*w[25];
          w[12] = a[0*sna_+1]*w[33]-a[0*sna_+2]*w[30]+a[0*sna_+4]*w[27]-a[0*sna_+5]*w[26];
          w[13] = a[0*sna_+1]*w[34]-a[0*sna_+3]*w[30]+a[0*sna_+4]*w[29]-a[0*sna_+5]*w[28];
          w[14] = a[0*sna_+2]*w[34]-a[0*sna_+3]*w[33]+a[0*sna_+4]*w[32]-a[0*sna_+5]*w[31];
          b[0*snb_+4] = d*( a[5*sna_+1]*w[14]-a[5*sna_+2]*w[13]+a[5*sna_+3]*w[12]-a[5*sna_+4]*w[11]+a[5*sna_+5]*w[10]);
          b[1*snb_+4] = d*(-a[5*sna_+0]*w[14]+a[5*sna_+2]*w[9]-a[5*sna_+3]*w[8]+a[5*sna_+4]*w[7]-a[5*sna_+5]*w[6]);
          b[2*snb_+4] = d*( a[5*sna_+0]*w[13]-a[5*sna_+1]*w[9]+a[5*sna_+3]*w[5]-a[5*sna_+4]*w[4]+a[5*sna_+5]*w[3]);
          b[3*snb_+4] = d*(-a[5*sna_+0]*w[12]+a[5*sna_+1]*w[8]-a[5*sna_+2]*w[5]+a[5*sna_+4]*w[2]-a[5*sna_+5]*w[1]);
          b[4*snb_+4] = d*( a[5*sna_+0]*w[11]-a[5*sna_+1]*w[7]+a[5*sna_+2]*w[4]-a[5*sna_+3]*w[2]+a[5*sna_+5]*w[0]);
          b[5*snb_+4] = d*(-a[5*sna_+0]*w[10]+a[5*sna_+1]*w[6]-a[5*sna_+2]*w[3]+a[5*sna_+3]*w[1]-a[5*sna_+4]*w[0]);
          b[0*snb_+5] = d*(-a[4*sna_+1]*w[14]+a[4*sna_+2]*w[13]-a[4*sna_+3]*w[12]+a[4*sna_+4]*w[11]-a[4*sna_+5]*w[10]);
          b[1*snb_+5] = d*( a[4*sna_+0]*w[14]-a[4*sna_+2]*w[9]+a[4*sna_+3]*w[8]-a[4*sna_+4]*w[7]+a[4*sna_+5]*w[6]);
          b[2*snb_+5] = d*(-a[4*sna_+0]*w[13]+a[4*sna_+1]*w[9]-a[4*sna_+3]*w[5]+a[4*sna_+4]*w[4]-a[4*sna_+5]*w[3]);
          b[3*snb_+5] = d*( a[4*sna_+0]*w[12]-a[4*sna_+1]*w[8]+a[4*sna_+2]*w[5]-a[4*sna_+4]*w[2]+a[4*sna_+5]*w[1]);
          b[4*snb_+5] = d*(-a[4*sna_+0]*w[11]+a[4*sna_+1]*w[7]-a[4*sna_+2]*w[4]+a[4*sna_+3]*w[2]-a[4*sna_+5]*w[0]);
          b[5*snb_+5] = d*( a[4*sna_+0]*w[10]-a[4*sna_+1]*w[6]+a[4*sna_+2]*w[3]-a[4*sna_+3]*w[1]+a[4*sna_+4]*w[0]);
          return det;
        }

        #ifdef __CUDACC__
        template<typename T_, typename ThreadGroup_>
        CUDA_DEVICE static __forceinline__ void grouped_compute(const ThreadGroup_& tg, T_* b, const T_* a, const T_&)
        {
          constexpr int n_ = 6;
          // copy matrix a to b
          for (int i = tg.thread_rank(); i < n_*n_; i += tg.num_threads())
          {
            const int row = i/n_;
            const int col = i%n_;
            b[row*snb_+col] = a[row*sna_+col];
          }

          // create shared pivot array
          __shared__ int p[n_];
          for (int i = tg.thread_rank(); i < n_; i += tg.num_threads())
          {
            p[i] = 0;
          }
          tg.sync();
          // call grouped invert from first warp subgroup
          if(tg.thread_rank() < 32)
          {
            auto a_g = cg::coalesced_threads();

            CudaMath::cuda_grouped_invert_matrix(a_g, n_, snb_, b, p);
          }
          tg.sync();

        }
        #endif
      };

      // generic square matrix inversion:
      template<int n_, int sna_, int snb_>
      struct InverseHelper<n_, n_, sna_, snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static T_ compute(T_* b, const T_* a)
        {
          // copy matrix a to b
          for (int i(0); i < n_; ++i)
          {
            for (int j(0); j < n_; ++j)
            {
              b[i*snb_+j] = a[i*sna_+j];
            }
          }

          // create pivot array
          int p[n_];

          // perform matrix inversion
          #ifndef __CUDACC__
          return Math::invert_matrix(n_, snb_, b, p);
          #else
          return CudaMath::cuda_invert_matrix(n_, snb_, b, p);
          #endif
        }

        #ifdef __CUDACC__
        template<typename T_, typename ThreadGroup_>
        CUDA_DEVICE static __forceinline__ void grouped_compute(const ThreadGroup_& tg, T_* b, const T_* a, const T_&)
        {
          // copy matrix a to b
          for (int i = tg.thread_rank(); i < n_*n_; i += tg.num_threads())
          {
            const int row = i/n_;
            const int col = i%n_;
            b[row*snb_+col] = a[row*sna_+col];
          }

          // create shared pivot array
          __shared__ int p[n_];
          for (int i = tg.thread_rank(); i < n_; i += tg.num_threads())
          {
            p[i] = 0;
          }
          tg.sync();
          // call grouped invert from first warp subgroup
          if(tg.thread_rank() < 32)
          {
            auto a_g = cg::coalesced_threads();

            CudaMath::cuda_grouped_invert_matrix(a_g, n_, snb_, b, p);
          }
          tg.sync();

        }
        #endif
      };

      template<int m_, int n_, int sna_, int snb_>
      struct CofactorHelper
      {
        template<typename T_>
        CUDA_HOST_DEVICE static void compute(T_* b, const T_* a)
        {
          static_assert(m_ == n_, "cofactor computation is only available for square matrices!");

          // compute inverse
          const T_ det = Intern::InverseHelper<m_,n_,sna_,snb_>::compute(b, a);

          for(int i(0); i < m_; ++i)
          {
            for(int j(0); j <= i; ++j)
            {
              b[i*n_+j] *= det;
            }
            for(int j(i+1); j < n_; ++j)
            {
              std::swap(b[i*n_+j], b[j*m_+i]);
              b[i*n_+j] *= det;
            }
          }
        }
      };

      template<int sna_, int snb_>
      struct CofactorHelper<1,1, sna_, snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static void compute(T_* b, const T_*)
        {
          b[0] = T_(1);
        }
      };

      template<int sna_, int snb_>
      struct CofactorHelper<2,2, sna_, snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static void compute(T_* b, const T_* a)
        {
          b[0*snb_+0] =  a[1*sna_+1];
          b[0*snb_+1] = -a[0*sna_+1];
          b[1*snb_+0] = -a[1*sna_+0];
          b[1*snb_+1] =  a[0*sna_+0];
        }
      };

      template<int sna_, int snb_>
      struct CofactorHelper<3,3, sna_, snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static void compute(T_* b, const T_* a)
        {
          b[0*snb_+0] = a[1*sna_+1]*a[2*sna_+2] - a[1*sna_+2]*a[2*sna_+1];
          b[1*snb_+0] = a[1*sna_+2]*a[2*sna_+0] - a[1*sna_+0]*a[2*sna_+2];
          b[2*snb_+0] = a[1*sna_+0]*a[2*sna_+1] - a[1*sna_+1]*a[2*sna_+0];
          b[0*snb_+1] = a[0*sna_+2]*a[2*sna_+1] - a[0*sna_+1]*a[2*sna_+2];
          b[1*snb_+1] = a[0*sna_+0]*a[2*sna_+2] - a[0*sna_+2]*a[2*sna_+0];
          b[2*snb_+1] = a[0*sna_+1]*a[2*sna_+0] - a[0*sna_+0]*a[2*sna_+1];
          b[0*snb_+2] = a[0*sna_+1]*a[1*sna_+2] - a[0*sna_+2]*a[1*sna_+1];
          b[1*snb_+2] = a[0*sna_+2]*a[1*sna_+0] - a[0*sna_+0]*a[1*sna_+2];
          b[2*snb_+2] = a[0*sna_+0]*a[1*sna_+1] - a[0*sna_+1]*a[1*sna_+0];
        }
      };

      template<int sna_, int snb_>
      struct CofactorHelper<4,4, sna_, snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static void compute(T_* b, const T_* a)
        {
          T_ w[6];
          w[0] = a[2*sna_+0]*a[3*sna_+1]-a[2*sna_+1]*a[3*sna_+0];
          w[1] = a[2*sna_+0]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+0];
          w[2] = a[2*sna_+0]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+0];
          w[3] = a[2*sna_+1]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+1];
          w[4] = a[2*sna_+1]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+1];
          w[5] = a[2*sna_+2]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+2];

          b[0*snb_+0] = a[1*sna_+1]*w[5]-a[1*sna_+2]*w[4]+a[1*sna_+3]*w[3];
          b[1*snb_+0] =-a[1*sna_+0]*w[5]+a[1*sna_+2]*w[2]-a[1*sna_+3]*w[1];
          b[2*snb_+0] = a[1*sna_+0]*w[4]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          b[3*snb_+0] =-a[1*sna_+0]*w[3]+a[1*sna_+1]*w[1]-a[1*sna_+2]*w[0];

          b[0*snb_+1] =-a[0*sna_+1]*w[5]+a[0*sna_+2]*w[4]-a[0*sna_+3]*w[3];
          b[1*snb_+1] = a[0*sna_+0]*w[5]-a[0*sna_+2]*w[2]+a[0*sna_+3]*w[1];
          b[2*snb_+1] =-a[0*sna_+0]*w[4]+a[0*sna_+1]*w[2]-a[0*sna_+3]*w[0];
          b[3*snb_+1] = a[0*sna_+0]*w[3]-a[0*sna_+1]*w[1]+a[0*sna_+2]*w[0];

          w[0] = a[0*sna_+0]*a[1*sna_+1]-a[0*sna_+1]*a[1*sna_+0];
          w[1] = a[0*sna_+0]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+0];
          w[2] = a[0*sna_+0]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+0];
          w[3] = a[0*sna_+1]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+1];
          w[4] = a[0*sna_+1]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+1];
          w[5] = a[0*sna_+2]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+2];

          b[0*snb_+2] = a[3*sna_+1]*w[5]-a[3*sna_+2]*w[4]+a[3*sna_+3]*w[3];
          b[1*snb_+2] =-a[3*sna_+0]*w[5]+a[3*sna_+2]*w[2]-a[3*sna_+3]*w[1];
          b[2*snb_+2] = a[3*sna_+0]*w[4]-a[3*sna_+1]*w[2]+a[3*sna_+3]*w[0];
          b[3*snb_+2] =-a[3*sna_+0]*w[3]+a[3*sna_+1]*w[1]-a[3*sna_+2]*w[0];

          b[0*snb_+3] =-a[2*sna_+1]*w[5]+a[2*sna_+2]*w[4]-a[2*sna_+3]*w[3];
          b[1*snb_+3] = a[2*sna_+0]*w[5]-a[2*sna_+2]*w[2]+a[2*sna_+3]*w[1];
          b[2*snb_+3] =-a[2*sna_+0]*w[4]+a[2*sna_+1]*w[2]-a[2*sna_+3]*w[0];
          b[3*snb_+3] = a[2*sna_+0]*w[3]-a[2*sna_+1]*w[1]+a[2*sna_+2]*w[0];
        }
      };

      template<int sna_, int snb_>
      struct CofactorHelper<5,5, sna_, snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static void compute(T_* b, const T_* a)
        {
          T_ w[20];
          w[ 0] = a[3*sna_+0]*a[4*sna_+1]-a[3*sna_+1]*a[4*sna_+0];
          w[ 1] = a[3*sna_+0]*a[4*sna_+2]-a[3*sna_+2]*a[4*sna_+0];
          w[ 2] = a[3*sna_+0]*a[4*sna_+3]-a[3*sna_+3]*a[4*sna_+0];
          w[ 3] = a[3*sna_+0]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+0];
          w[ 4] = a[3*sna_+1]*a[4*sna_+2]-a[3*sna_+2]*a[4*sna_+1];
          w[ 5] = a[3*sna_+1]*a[4*sna_+3]-a[3*sna_+3]*a[4*sna_+1];
          w[ 6] = a[3*sna_+1]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+1];
          w[ 7] = a[3*sna_+2]*a[4*sna_+3]-a[3*sna_+3]*a[4*sna_+2];
          w[ 8] = a[3*sna_+2]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+2];
          w[ 9] = a[3*sna_+3]*a[4*sna_+4]-a[3*sna_+4]*a[4*sna_+3];
          w[10] = a[2*sna_+0]*w[4]-a[2*sna_+1]*w[1]+a[2*sna_+2]*w[0];
          w[11] = a[2*sna_+0]*w[5]-a[2*sna_+1]*w[2]+a[2*sna_+3]*w[0];
          w[12] = a[2*sna_+0]*w[6]-a[2*sna_+1]*w[3]+a[2*sna_+4]*w[0];
          w[13] = a[2*sna_+0]*w[7]-a[2*sna_+2]*w[2]+a[2*sna_+3]*w[1];
          w[14] = a[2*sna_+0]*w[8]-a[2*sna_+2]*w[3]+a[2*sna_+4]*w[1];
          w[15] = a[2*sna_+0]*w[9]-a[2*sna_+3]*w[3]+a[2*sna_+4]*w[2];
          w[16] = a[2*sna_+1]*w[7]-a[2*sna_+2]*w[5]+a[2*sna_+3]*w[4];
          w[17] = a[2*sna_+1]*w[8]-a[2*sna_+2]*w[6]+a[2*sna_+4]*w[4];
          w[18] = a[2*sna_+1]*w[9]-a[2*sna_+3]*w[6]+a[2*sna_+4]*w[5];
          w[19] = a[2*sna_+2]*w[9]-a[2*sna_+3]*w[8]+a[2*sna_+4]*w[7];

          b[0*snb_+0] = a[1*sna_+1]*w[19]-a[1*sna_+2]*w[18]+a[1*sna_+3]*w[17]-a[1*sna_+4]*w[16];
          b[1*snb_+0] =-a[1*sna_+0]*w[19]+a[1*sna_+2]*w[15]-a[1*sna_+3]*w[14]+a[1*sna_+4]*w[13];
          b[2*snb_+0] = a[1*sna_+0]*w[18]-a[1*sna_+1]*w[15]+a[1*sna_+3]*w[12]-a[1*sna_+4]*w[11];
          b[3*snb_+0] =-a[1*sna_+0]*w[17]+a[1*sna_+1]*w[14]-a[1*sna_+2]*w[12]+a[1*sna_+4]*w[10];
          b[4*snb_+0] = a[1*sna_+0]*w[16]-a[1*sna_+1]*w[13]+a[1*sna_+2]*w[11]-a[1*sna_+3]*w[10];

          b[0*snb_+1] =-a[0*sna_+1]*w[19]+a[0*sna_+2]*w[18]-a[0*sna_+3]*w[17]+a[0*sna_+4]*w[16];
          b[1*snb_+1] = a[0*sna_+0]*w[19]-a[0*sna_+2]*w[15]+a[0*sna_+3]*w[14]-a[0*sna_+4]*w[13];
          b[2*snb_+1] =-a[0*sna_+0]*w[18]+a[0*sna_+1]*w[15]-a[0*sna_+3]*w[12]+a[0*sna_+4]*w[11];
          b[3*snb_+1] = a[0*sna_+0]*w[17]-a[0*sna_+1]*w[14]+a[0*sna_+2]*w[12]-a[0*sna_+4]*w[10];
          b[4*snb_+1] =-a[0*sna_+0]*w[16]+a[0*sna_+1]*w[13]-a[0*sna_+2]*w[11]+a[0*sna_+3]*w[10];

          w[10] = a[1*sna_+0]*w[4]-a[1*sna_+1]*w[1]+a[1*sna_+2]*w[0];
          w[11] = a[1*sna_+0]*w[5]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          w[12] = a[1*sna_+0]*w[6]-a[1*sna_+1]*w[3]+a[1*sna_+4]*w[0];
          w[13] = a[1*sna_+0]*w[7]-a[1*sna_+2]*w[2]+a[1*sna_+3]*w[1];
          w[14] = a[1*sna_+0]*w[8]-a[1*sna_+2]*w[3]+a[1*sna_+4]*w[1];
          w[15] = a[1*sna_+0]*w[9]-a[1*sna_+3]*w[3]+a[1*sna_+4]*w[2];
          w[16] = a[1*sna_+1]*w[7]-a[1*sna_+2]*w[5]+a[1*sna_+3]*w[4];
          w[17] = a[1*sna_+1]*w[8]-a[1*sna_+2]*w[6]+a[1*sna_+4]*w[4];
          w[18] = a[1*sna_+1]*w[9]-a[1*sna_+3]*w[6]+a[1*sna_+4]*w[5];
          w[19] = a[1*sna_+2]*w[9]-a[1*sna_+3]*w[8]+a[1*sna_+4]*w[7];

          b[0*snb_+2] = a[0*sna_+1]*w[19]-a[0*sna_+2]*w[18]+a[0*sna_+3]*w[17]-a[0*sna_+4]*w[16];
          b[1*snb_+2] =-a[0*sna_+0]*w[19]+a[0*sna_+2]*w[15]-a[0*sna_+3]*w[14]+a[0*sna_+4]*w[13];
          b[2*snb_+2] = a[0*sna_+0]*w[18]-a[0*sna_+1]*w[15]+a[0*sna_+3]*w[12]-a[0*sna_+4]*w[11];
          b[3*snb_+2] =-a[0*sna_+0]*w[17]+a[0*sna_+1]*w[14]-a[0*sna_+2]*w[12]+a[0*sna_+4]*w[10];
          b[4*snb_+2] = a[0*sna_+0]*w[16]-a[0*sna_+1]*w[13]+a[0*sna_+2]*w[11]-a[0*sna_+3]*w[10];

          w[ 0] = a[0*sna_+0]*a[1*sna_+1]-a[0*sna_+1]*a[1*sna_+0];
          w[ 1] = a[0*sna_+0]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+0];
          w[ 2] = a[0*sna_+0]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+0];
          w[ 3] = a[0*sna_+0]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+0];
          w[ 4] = a[0*sna_+1]*a[1*sna_+2]-a[0*sna_+2]*a[1*sna_+1];
          w[ 5] = a[0*sna_+1]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+1];
          w[ 6] = a[0*sna_+1]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+1];
          w[ 7] = a[0*sna_+2]*a[1*sna_+3]-a[0*sna_+3]*a[1*sna_+2];
          w[ 8] = a[0*sna_+2]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+2];
          w[ 9] = a[0*sna_+3]*a[1*sna_+4]-a[0*sna_+4]*a[1*sna_+3];
          w[10] = a[2*sna_+0]*w[4]-a[2*sna_+1]*w[1]+a[2*sna_+2]*w[0];
          w[11] = a[2*sna_+0]*w[5]-a[2*sna_+1]*w[2]+a[2*sna_+3]*w[0];
          w[12] = a[2*sna_+0]*w[6]-a[2*sna_+1]*w[3]+a[2*sna_+4]*w[0];
          w[13] = a[2*sna_+0]*w[7]-a[2*sna_+2]*w[2]+a[2*sna_+3]*w[1];
          w[14] = a[2*sna_+0]*w[8]-a[2*sna_+2]*w[3]+a[2*sna_+4]*w[1];
          w[15] = a[2*sna_+0]*w[9]-a[2*sna_+3]*w[3]+a[2*sna_+4]*w[2];
          w[16] = a[2*sna_+1]*w[7]-a[2*sna_+2]*w[5]+a[2*sna_+3]*w[4];
          w[17] = a[2*sna_+1]*w[8]-a[2*sna_+2]*w[6]+a[2*sna_+4]*w[4];
          w[18] = a[2*sna_+1]*w[9]-a[2*sna_+3]*w[6]+a[2*sna_+4]*w[5];
          w[19] = a[2*sna_+2]*w[9]-a[2*sna_+3]*w[8]+a[2*sna_+4]*w[7];

          b[0*snb_+3] =  a[4*sna_+1]*w[19]-a[4*sna_+2]*w[18]+a[4*sna_+3]*w[17]-a[4*sna_+4]*w[16];
          b[1*snb_+3] = -a[4*sna_+0]*w[19]+a[4*sna_+2]*w[15]-a[4*sna_+3]*w[14]+a[4*sna_+4]*w[13];
          b[2*snb_+3] =  a[4*sna_+0]*w[18]-a[4*sna_+1]*w[15]+a[4*sna_+3]*w[12]-a[4*sna_+4]*w[11];
          b[3*snb_+3] = -a[4*sna_+0]*w[17]+a[4*sna_+1]*w[14]-a[4*sna_+2]*w[12]+a[4*sna_+4]*w[10];
          b[4*snb_+3] =  a[4*sna_+0]*w[16]-a[4*sna_+1]*w[13]+a[4*sna_+2]*w[11]-a[4*sna_+3]*w[10];

          b[0*snb_+4] = -a[3*sna_+1]*w[19]+a[3*sna_+2]*w[18]-a[3*sna_+3]*w[17]+a[3*sna_+4]*w[16];
          b[1*snb_+4] =  a[3*sna_+0]*w[19]-a[3*sna_+2]*w[15]+a[3*sna_+3]*w[14]-a[3*sna_+4]*w[13];
          b[2*snb_+4] = -a[3*sna_+0]*w[18]+a[3*sna_+1]*w[15]-a[3*sna_+3]*w[12]+a[3*sna_+4]*w[11];
          b[3*snb_+4] =  a[3*sna_+0]*w[17]-a[3*sna_+1]*w[14]+a[3*sna_+2]*w[12]-a[3*sna_+4]*w[10];
          b[4*snb_+4] = -a[3*sna_+0]*w[16]+a[3*sna_+1]*w[13]-a[3*sna_+2]*w[11]+a[3*sna_+3]*w[10];
        }
      };

      template<int sna_, int snb_>
      struct CofactorHelper<6,6, sna_, snb_>
      {
        template<typename T_>
        CUDA_HOST_DEVICE static void compute(T_* b, const T_* a)
        {
          T_ w[35];
          w[ 0] = a[4*sna_+0]*a[5*sna_+1]-a[4*sna_+1]*a[5*sna_+0];
          w[ 1] = a[4*sna_+0]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+0];
          w[ 2] = a[4*sna_+0]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+0];
          w[ 3] = a[4*sna_+0]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+0];
          w[ 4] = a[4*sna_+0]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+0];
          w[ 5] = a[4*sna_+1]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+1];
          w[ 6] = a[4*sna_+1]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+1];
          w[ 7] = a[4*sna_+1]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+1];
          w[ 8] = a[4*sna_+1]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+1];
          w[ 9] = a[4*sna_+2]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+2];
          w[10] = a[4*sna_+2]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+2];
          w[11] = a[4*sna_+2]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+2];
          w[12] = a[4*sna_+3]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+3];
          w[13] = a[4*sna_+3]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+3];
          w[14] = a[4*sna_+4]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+4];
          w[15] = a[3*sna_+0]*w[5]-a[3*sna_+1]*w[1]+a[3*sna_+2]*w[0];
          w[16] = a[3*sna_+0]*w[6]-a[3*sna_+1]*w[2]+a[3*sna_+3]*w[0];
          w[17] = a[3*sna_+0]*w[7]-a[3*sna_+1]*w[3]+a[3*sna_+4]*w[0];
          w[18] = a[3*sna_+0]*w[8]-a[3*sna_+1]*w[4]+a[3*sna_+5]*w[0];
          w[19] = a[3*sna_+0]*w[9]-a[3*sna_+2]*w[2]+a[3*sna_+3]*w[1];
          w[20] = a[3*sna_+0]*w[10]-a[3*sna_+2]*w[3]+a[3*sna_+4]*w[1];
          w[21] = a[3*sna_+0]*w[11]-a[3*sna_+2]*w[4]+a[3*sna_+5]*w[1];
          w[22] = a[3*sna_+0]*w[12]-a[3*sna_+3]*w[3]+a[3*sna_+4]*w[2];
          w[23] = a[3*sna_+0]*w[13]-a[3*sna_+3]*w[4]+a[3*sna_+5]*w[2];
          w[24] = a[3*sna_+0]*w[14]-a[3*sna_+4]*w[4]+a[3*sna_+5]*w[3];
          w[25] = a[3*sna_+1]*w[9]-a[3*sna_+2]*w[6]+a[3*sna_+3]*w[5];
          w[26] = a[3*sna_+1]*w[10]-a[3*sna_+2]*w[7]+a[3*sna_+4]*w[5];
          w[27] = a[3*sna_+1]*w[11]-a[3*sna_+2]*w[8]+a[3*sna_+5]*w[5];
          w[28] = a[3*sna_+1]*w[12]-a[3*sna_+3]*w[7]+a[3*sna_+4]*w[6];
          w[29] = a[3*sna_+1]*w[13]-a[3*sna_+3]*w[8]+a[3*sna_+5]*w[6];
          w[30] = a[3*sna_+1]*w[14]-a[3*sna_+4]*w[8]+a[3*sna_+5]*w[7];
          w[31] = a[3*sna_+2]*w[12]-a[3*sna_+3]*w[10]+a[3*sna_+4]*w[9];
          w[32] = a[3*sna_+2]*w[13]-a[3*sna_+3]*w[11]+a[3*sna_+5]*w[9];
          w[33] = a[3*sna_+2]*w[14]-a[3*sna_+4]*w[11]+a[3*sna_+5]*w[10];
          w[34] = a[3*sna_+3]*w[14]-a[3*sna_+4]*w[13]+a[3*sna_+5]*w[12];

          w[ 0] = a[2*sna_+0]*w[25]-a[2*sna_+1]*w[19]+a[2*sna_+2]*w[16]-a[2*sna_+3]*w[15];
          w[ 1] = a[2*sna_+0]*w[26]-a[2*sna_+1]*w[20]+a[2*sna_+2]*w[17]-a[2*sna_+4]*w[15];
          w[ 2] = a[2*sna_+0]*w[27]-a[2*sna_+1]*w[21]+a[2*sna_+2]*w[18]-a[2*sna_+5]*w[15];
          w[ 3] = a[2*sna_+0]*w[28]-a[2*sna_+1]*w[22]+a[2*sna_+3]*w[17]-a[2*sna_+4]*w[16];
          w[ 4] = a[2*sna_+0]*w[29]-a[2*sna_+1]*w[23]+a[2*sna_+3]*w[18]-a[2*sna_+5]*w[16];
          w[ 5] = a[2*sna_+0]*w[30]-a[2*sna_+1]*w[24]+a[2*sna_+4]*w[18]-a[2*sna_+5]*w[17];
          w[ 6] = a[2*sna_+0]*w[31]-a[2*sna_+2]*w[22]+a[2*sna_+3]*w[20]-a[2*sna_+4]*w[19];
          w[ 7] = a[2*sna_+0]*w[32]-a[2*sna_+2]*w[23]+a[2*sna_+3]*w[21]-a[2*sna_+5]*w[19];
          w[ 8] = a[2*sna_+0]*w[33]-a[2*sna_+2]*w[24]+a[2*sna_+4]*w[21]-a[2*sna_+5]*w[20];
          w[ 9] = a[2*sna_+0]*w[34]-a[2*sna_+3]*w[24]+a[2*sna_+4]*w[23]-a[2*sna_+5]*w[22];
          w[10] = a[2*sna_+1]*w[31]-a[2*sna_+2]*w[28]+a[2*sna_+3]*w[26]-a[2*sna_+4]*w[25];
          w[11] = a[2*sna_+1]*w[32]-a[2*sna_+2]*w[29]+a[2*sna_+3]*w[27]-a[2*sna_+5]*w[25];
          w[12] = a[2*sna_+1]*w[33]-a[2*sna_+2]*w[30]+a[2*sna_+4]*w[27]-a[2*sna_+5]*w[26];
          w[13] = a[2*sna_+1]*w[34]-a[2*sna_+3]*w[30]+a[2*sna_+4]*w[29]-a[2*sna_+5]*w[28];
          w[14] = a[2*sna_+2]*w[34]-a[2*sna_+3]*w[33]+a[2*sna_+4]*w[32]-a[2*sna_+5]*w[31];

          b[0*snb_+0] =  a[1*sna_+1]*w[14]-a[1*sna_+2]*w[13]+a[1*sna_+3]*w[12]-a[1*sna_+4]*w[11]+a[1*sna_+5]*w[10];
          b[1*snb_+0] = -a[1*sna_+0]*w[14]+a[1*sna_+2]*w[9]-a[1*sna_+3]*w[8]+a[1*sna_+4]*w[7]-a[1*sna_+5]*w[6];
          b[2*snb_+0] =  a[1*sna_+0]*w[13]-a[1*sna_+1]*w[9]+a[1*sna_+3]*w[5]-a[1*sna_+4]*w[4]+a[1*sna_+5]*w[3];
          b[3*snb_+0] = -a[1*sna_+0]*w[12]+a[1*sna_+1]*w[8]-a[1*sna_+2]*w[5]+a[1*sna_+4]*w[2]-a[1*sna_+5]*w[1];
          b[4*snb_+0] =  a[1*sna_+0]*w[11]-a[1*sna_+1]*w[7]+a[1*sna_+2]*w[4]-a[1*sna_+3]*w[2]+a[1*sna_+5]*w[0];
          b[5*snb_+0] = -a[1*sna_+0]*w[10]+a[1*sna_+1]*w[6]-a[1*sna_+2]*w[3]+a[1*sna_+3]*w[1]-a[1*sna_+4]*w[0];

          b[0*snb_+1] =-a[0*sna_+1]*w[14]+a[0*sna_+2]*w[13]-a[0*sna_+3]*w[12]+a[0*sna_+4]*w[11]-a[0*sna_+5]*w[10];
          b[1*snb_+1] = a[0*sna_+0]*w[14]-a[0*sna_+2]*w[9]+a[0*sna_+3]*w[8]-a[0*sna_+4]*w[7]+a[0*sna_+5]*w[6];
          b[2*snb_+1] =-a[0*sna_+0]*w[13]+a[0*sna_+1]*w[9]-a[0*sna_+3]*w[5]+a[0*sna_+4]*w[4]-a[0*sna_+5]*w[3];
          b[3*snb_+1] = a[0*sna_+0]*w[12]-a[0*sna_+1]*w[8]+a[0*sna_+2]*w[5]-a[0*sna_+4]*w[2]+a[0*sna_+5]*w[1];
          b[4*snb_+1] =-a[0*sna_+0]*w[11]+a[0*sna_+1]*w[7]-a[0*sna_+2]*w[4]+a[0*sna_+3]*w[2]-a[0*sna_+5]*w[0];
          b[5*snb_+1] = a[0*sna_+0]*w[10]-a[0*sna_+1]*w[6]+a[0*sna_+2]*w[3]-a[0*sna_+3]*w[1]+a[0*sna_+4]*w[0];

          w[ 0] = a[4*sna_+0]*a[5*sna_+1]-a[4*sna_+1]*a[5*sna_+0];
          w[ 1] = a[4*sna_+0]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+0];
          w[ 2] = a[4*sna_+0]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+0];
          w[ 3] = a[4*sna_+0]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+0];
          w[ 4] = a[4*sna_+0]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+0];
          w[ 5] = a[4*sna_+1]*a[5*sna_+2]-a[4*sna_+2]*a[5*sna_+1];
          w[ 6] = a[4*sna_+1]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+1];
          w[ 7] = a[4*sna_+1]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+1];
          w[ 8] = a[4*sna_+1]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+1];
          w[ 9] = a[4*sna_+2]*a[5*sna_+3]-a[4*sna_+3]*a[5*sna_+2];
          w[10] = a[4*sna_+2]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+2];
          w[11] = a[4*sna_+2]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+2];
          w[12] = a[4*sna_+3]*a[5*sna_+4]-a[4*sna_+4]*a[5*sna_+3];
          w[13] = a[4*sna_+3]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+3];
          w[14] = a[4*sna_+4]*a[5*sna_+5]-a[4*sna_+5]*a[5*sna_+4];
          w[15] = a[1*sna_+0]*w[5]-a[1*sna_+1]*w[1]+a[1*sna_+2]*w[0];
          w[16] = a[1*sna_+0]*w[6]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          w[17] = a[1*sna_+0]*w[7]-a[1*sna_+1]*w[3]+a[1*sna_+4]*w[0];
          w[18] = a[1*sna_+0]*w[8]-a[1*sna_+1]*w[4]+a[1*sna_+5]*w[0];
          w[19] = a[1*sna_+0]*w[9]-a[1*sna_+2]*w[2]+a[1*sna_+3]*w[1];
          w[20] = a[1*sna_+0]*w[10]-a[1*sna_+2]*w[3]+a[1*sna_+4]*w[1];
          w[21] = a[1*sna_+0]*w[11]-a[1*sna_+2]*w[4]+a[1*sna_+5]*w[1];
          w[22] = a[1*sna_+0]*w[12]-a[1*sna_+3]*w[3]+a[1*sna_+4]*w[2];
          w[23] = a[1*sna_+0]*w[13]-a[1*sna_+3]*w[4]+a[1*sna_+5]*w[2];
          w[24] = a[1*sna_+0]*w[14]-a[1*sna_+4]*w[4]+a[1*sna_+5]*w[3];
          w[25] = a[1*sna_+1]*w[9]-a[1*sna_+2]*w[6]+a[1*sna_+3]*w[5];
          w[26] = a[1*sna_+1]*w[10]-a[1*sna_+2]*w[7]+a[1*sna_+4]*w[5];
          w[27] = a[1*sna_+1]*w[11]-a[1*sna_+2]*w[8]+a[1*sna_+5]*w[5];
          w[28] = a[1*sna_+1]*w[12]-a[1*sna_+3]*w[7]+a[1*sna_+4]*w[6];
          w[29] = a[1*sna_+1]*w[13]-a[1*sna_+3]*w[8]+a[1*sna_+5]*w[6];
          w[30] = a[1*sna_+1]*w[14]-a[1*sna_+4]*w[8]+a[1*sna_+5]*w[7];
          w[31] = a[1*sna_+2]*w[12]-a[1*sna_+3]*w[10]+a[1*sna_+4]*w[9];
          w[32] = a[1*sna_+2]*w[13]-a[1*sna_+3]*w[11]+a[1*sna_+5]*w[9];
          w[33] = a[1*sna_+2]*w[14]-a[1*sna_+4]*w[11]+a[1*sna_+5]*w[10];
          w[34] = a[1*sna_+3]*w[14]-a[1*sna_+4]*w[13]+a[1*sna_+5]*w[12];
          w[ 0] = a[0*sna_+0]*w[25]-a[0*sna_+1]*w[19]+a[0*sna_+2]*w[16]-a[0*sna_+3]*w[15];
          w[ 1] = a[0*sna_+0]*w[26]-a[0*sna_+1]*w[20]+a[0*sna_+2]*w[17]-a[0*sna_+4]*w[15];
          w[ 2] = a[0*sna_+0]*w[27]-a[0*sna_+1]*w[21]+a[0*sna_+2]*w[18]-a[0*sna_+5]*w[15];
          w[ 3] = a[0*sna_+0]*w[28]-a[0*sna_+1]*w[22]+a[0*sna_+3]*w[17]-a[0*sna_+4]*w[16];
          w[ 4] = a[0*sna_+0]*w[29]-a[0*sna_+1]*w[23]+a[0*sna_+3]*w[18]-a[0*sna_+5]*w[16];
          w[ 5] = a[0*sna_+0]*w[30]-a[0*sna_+1]*w[24]+a[0*sna_+4]*w[18]-a[0*sna_+5]*w[17];
          w[ 6] = a[0*sna_+0]*w[31]-a[0*sna_+2]*w[22]+a[0*sna_+3]*w[20]-a[0*sna_+4]*w[19];
          w[ 7] = a[0*sna_+0]*w[32]-a[0*sna_+2]*w[23]+a[0*sna_+3]*w[21]-a[0*sna_+5]*w[19];
          w[ 8] = a[0*sna_+0]*w[33]-a[0*sna_+2]*w[24]+a[0*sna_+4]*w[21]-a[0*sna_+5]*w[20];
          w[ 9] = a[0*sna_+0]*w[34]-a[0*sna_+3]*w[24]+a[0*sna_+4]*w[23]-a[0*sna_+5]*w[22];
          w[10] = a[0*sna_+1]*w[31]-a[0*sna_+2]*w[28]+a[0*sna_+3]*w[26]-a[0*sna_+4]*w[25];
          w[11] = a[0*sna_+1]*w[32]-a[0*sna_+2]*w[29]+a[0*sna_+3]*w[27]-a[0*sna_+5]*w[25];
          w[12] = a[0*sna_+1]*w[33]-a[0*sna_+2]*w[30]+a[0*sna_+4]*w[27]-a[0*sna_+5]*w[26];
          w[13] = a[0*sna_+1]*w[34]-a[0*sna_+3]*w[30]+a[0*sna_+4]*w[29]-a[0*sna_+5]*w[28];
          w[14] = a[0*sna_+2]*w[34]-a[0*sna_+3]*w[33]+a[0*sna_+4]*w[32]-a[0*sna_+5]*w[31];

          b[0*snb_+2] = a[3*sna_+1]*w[14]-a[3*sna_+2]*w[13]+a[3*sna_+3]*w[12]-a[3*sna_+4]*w[11]+a[3*sna_+5]*w[10];
          b[1*snb_+2] =-a[3*sna_+0]*w[14]+a[3*sna_+2]*w[9]-a[3*sna_+3]*w[8]+a[3*sna_+4]*w[7]-a[3*sna_+5]*w[6];
          b[2*snb_+2] = a[3*sna_+0]*w[13]-a[3*sna_+1]*w[9]+a[3*sna_+3]*w[5]-a[3*sna_+4]*w[4]+a[3*sna_+5]*w[3];
          b[3*snb_+2] =-a[3*sna_+0]*w[12]+a[3*sna_+1]*w[8]-a[3*sna_+2]*w[5]+a[3*sna_+4]*w[2]-a[3*sna_+5]*w[1];
          b[4*snb_+2] = a[3*sna_+0]*w[11]-a[3*sna_+1]*w[7]+a[3*sna_+2]*w[4]-a[3*sna_+3]*w[2]+a[3*sna_+5]*w[0];
          b[5*snb_+2] =-a[3*sna_+0]*w[10]+a[3*sna_+1]*w[6]-a[3*sna_+2]*w[3]+a[3*sna_+3]*w[1]-a[3*sna_+4]*w[0];

          b[0*snb_+3] =-a[2*sna_+1]*w[14]+a[2*sna_+2]*w[13]-a[2*sna_+3]*w[12]+a[2*sna_+4]*w[11]-a[2*sna_+5]*w[10];
          b[1*snb_+3] = a[2*sna_+0]*w[14]-a[2*sna_+2]*w[9]+a[2*sna_+3]*w[8]-a[2*sna_+4]*w[7]+a[2*sna_+5]*w[6];
          b[2*snb_+3] =-a[2*sna_+0]*w[13]+a[2*sna_+1]*w[9]-a[2*sna_+3]*w[5]+a[2*sna_+4]*w[4]-a[2*sna_+5]*w[3];
          b[3*snb_+3] = a[2*sna_+0]*w[12]-a[2*sna_+1]*w[8]+a[2*sna_+2]*w[5]-a[2*sna_+4]*w[2]+a[2*sna_+5]*w[1];
          b[4*snb_+3] =-a[2*sna_+0]*w[11]+a[2*sna_+1]*w[7]-a[2*sna_+2]*w[4]+a[2*sna_+3]*w[2]-a[2*sna_+5]*w[0];
          b[5*snb_+3] = a[2*sna_+0]*w[10]-a[2*sna_+1]*w[6]+a[2*sna_+2]*w[3]-a[2*sna_+3]*w[1]+a[2*sna_+4]*w[0];

          w[ 0] = a[2*sna_+0]*a[3*sna_+1]-a[2*sna_+1]*a[3*sna_+0];
          w[ 1] = a[2*sna_+0]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+0];
          w[ 2] = a[2*sna_+0]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+0];
          w[ 3] = a[2*sna_+0]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+0];
          w[ 4] = a[2*sna_+0]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+0];
          w[ 5] = a[2*sna_+1]*a[3*sna_+2]-a[2*sna_+2]*a[3*sna_+1];
          w[ 6] = a[2*sna_+1]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+1];
          w[ 7] = a[2*sna_+1]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+1];
          w[ 8] = a[2*sna_+1]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+1];
          w[ 9] = a[2*sna_+2]*a[3*sna_+3]-a[2*sna_+3]*a[3*sna_+2];
          w[10] = a[2*sna_+2]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+2];
          w[11] = a[2*sna_+2]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+2];
          w[12] = a[2*sna_+3]*a[3*sna_+4]-a[2*sna_+4]*a[3*sna_+3];
          w[13] = a[2*sna_+3]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+3];
          w[14] = a[2*sna_+4]*a[3*sna_+5]-a[2*sna_+5]*a[3*sna_+4];
          w[15] = a[1*sna_+0]*w[5]-a[1*sna_+1]*w[1]+a[1*sna_+2]*w[0];
          w[16] = a[1*sna_+0]*w[6]-a[1*sna_+1]*w[2]+a[1*sna_+3]*w[0];
          w[17] = a[1*sna_+0]*w[7]-a[1*sna_+1]*w[3]+a[1*sna_+4]*w[0];
          w[18] = a[1*sna_+0]*w[8]-a[1*sna_+1]*w[4]+a[1*sna_+5]*w[0];
          w[19] = a[1*sna_+0]*w[9]-a[1*sna_+2]*w[2]+a[1*sna_+3]*w[1];
          w[20] = a[1*sna_+0]*w[10]-a[1*sna_+2]*w[3]+a[1*sna_+4]*w[1];
          w[21] = a[1*sna_+0]*w[11]-a[1*sna_+2]*w[4]+a[1*sna_+5]*w[1];
          w[22] = a[1*sna_+0]*w[12]-a[1*sna_+3]*w[3]+a[1*sna_+4]*w[2];
          w[23] = a[1*sna_+0]*w[13]-a[1*sna_+3]*w[4]+a[1*sna_+5]*w[2];
          w[24] = a[1*sna_+0]*w[14]-a[1*sna_+4]*w[4]+a[1*sna_+5]*w[3];
          w[25] = a[1*sna_+1]*w[9]-a[1*sna_+2]*w[6]+a[1*sna_+3]*w[5];
          w[26] = a[1*sna_+1]*w[10]-a[1*sna_+2]*w[7]+a[1*sna_+4]*w[5];
          w[27] = a[1*sna_+1]*w[11]-a[1*sna_+2]*w[8]+a[1*sna_+5]*w[5];
          w[28] = a[1*sna_+1]*w[12]-a[1*sna_+3]*w[7]+a[1*sna_+4]*w[6];
          w[29] = a[1*sna_+1]*w[13]-a[1*sna_+3]*w[8]+a[1*sna_+5]*w[6];
          w[30] = a[1*sna_+1]*w[14]-a[1*sna_+4]*w[8]+a[1*sna_+5]*w[7];
          w[31] = a[1*sna_+2]*w[12]-a[1*sna_+3]*w[10]+a[1*sna_+4]*w[9];
          w[32] = a[1*sna_+2]*w[13]-a[1*sna_+3]*w[11]+a[1*sna_+5]*w[9];
          w[33] = a[1*sna_+2]*w[14]-a[1*sna_+4]*w[11]+a[1*sna_+5]*w[10];
          w[34] = a[1*sna_+3]*w[14]-a[1*sna_+4]*w[13]+a[1*sna_+5]*w[12];

          w[ 0] = a[0*sna_+0]*w[25]-a[0*sna_+1]*w[19]+a[0*sna_+2]*w[16]-a[0*sna_+3]*w[15];
          w[ 1] = a[0*sna_+0]*w[26]-a[0*sna_+1]*w[20]+a[0*sna_+2]*w[17]-a[0*sna_+4]*w[15];
          w[ 2] = a[0*sna_+0]*w[27]-a[0*sna_+1]*w[21]+a[0*sna_+2]*w[18]-a[0*sna_+5]*w[15];
          w[ 3] = a[0*sna_+0]*w[28]-a[0*sna_+1]*w[22]+a[0*sna_+3]*w[17]-a[0*sna_+4]*w[16];
          w[ 4] = a[0*sna_+0]*w[29]-a[0*sna_+1]*w[23]+a[0*sna_+3]*w[18]-a[0*sna_+5]*w[16];
          w[ 5] = a[0*sna_+0]*w[30]-a[0*sna_+1]*w[24]+a[0*sna_+4]*w[18]-a[0*sna_+5]*w[17];
          w[ 6] = a[0*sna_+0]*w[31]-a[0*sna_+2]*w[22]+a[0*sna_+3]*w[20]-a[0*sna_+4]*w[19];
          w[ 7] = a[0*sna_+0]*w[32]-a[0*sna_+2]*w[23]+a[0*sna_+3]*w[21]-a[0*sna_+5]*w[19];
          w[ 8] = a[0*sna_+0]*w[33]-a[0*sna_+2]*w[24]+a[0*sna_+4]*w[21]-a[0*sna_+5]*w[20];
          w[ 9] = a[0*sna_+0]*w[34]-a[0*sna_+3]*w[24]+a[0*sna_+4]*w[23]-a[0*sna_+5]*w[22];
          w[10] = a[0*sna_+1]*w[31]-a[0*sna_+2]*w[28]+a[0*sna_+3]*w[26]-a[0*sna_+4]*w[25];
          w[11] = a[0*sna_+1]*w[32]-a[0*sna_+2]*w[29]+a[0*sna_+3]*w[27]-a[0*sna_+5]*w[25];
          w[12] = a[0*sna_+1]*w[33]-a[0*sna_+2]*w[30]+a[0*sna_+4]*w[27]-a[0*sna_+5]*w[26];
          w[13] = a[0*sna_+1]*w[34]-a[0*sna_+3]*w[30]+a[0*sna_+4]*w[29]-a[0*sna_+5]*w[28];
          w[14] = a[0*sna_+2]*w[34]-a[0*sna_+3]*w[33]+a[0*sna_+4]*w[32]-a[0*sna_+5]*w[31];

          b[0*snb_+4] = a[5*sna_+1]*w[14]-a[5*sna_+2]*w[13]+a[5*sna_+3]*w[12]-a[5*sna_+4]*w[11]+a[5*sna_+5]*w[10];
          b[1*snb_+4] =-a[5*sna_+0]*w[14]+a[5*sna_+2]*w[9]-a[5*sna_+3]*w[8]+a[5*sna_+4]*w[7]-a[5*sna_+5]*w[6];
          b[2*snb_+4] = a[5*sna_+0]*w[13]-a[5*sna_+1]*w[9]+a[5*sna_+3]*w[5]-a[5*sna_+4]*w[4]+a[5*sna_+5]*w[3];
          b[3*snb_+4] =-a[5*sna_+0]*w[12]+a[5*sna_+1]*w[8]-a[5*sna_+2]*w[5]+a[5*sna_+4]*w[2]-a[5*sna_+5]*w[1];
          b[4*snb_+4] = a[5*sna_+0]*w[11]-a[5*sna_+1]*w[7]+a[5*sna_+2]*w[4]-a[5*sna_+3]*w[2]+a[5*sna_+5]*w[0];
          b[5*snb_+4] =-a[5*sna_+0]*w[10]+a[5*sna_+1]*w[6]-a[5*sna_+2]*w[3]+a[5*sna_+3]*w[1]-a[5*sna_+4]*w[0];

          b[0*snb_+5] =-a[4*sna_+1]*w[14]+a[4*sna_+2]*w[13]-a[4*sna_+3]*w[12]+a[4*sna_+4]*w[11]-a[4*sna_+5]*w[10];
          b[1*snb_+5] = a[4*sna_+0]*w[14]-a[4*sna_+2]*w[9]+a[4*sna_+3]*w[8]-a[4*sna_+4]*w[7]+a[4*sna_+5]*w[6];
          b[2*snb_+5] =-a[4*sna_+0]*w[13]+a[4*sna_+1]*w[9]-a[4*sna_+3]*w[5]+a[4*sna_+4]*w[4]-a[4*sna_+5]*w[3];
          b[3*snb_+5] = a[4*sna_+0]*w[12]-a[4*sna_+1]*w[8]+a[4*sna_+2]*w[5]-a[4*sna_+4]*w[2]+a[4*sna_+5]*w[1];
          b[4*snb_+5] =-a[4*sna_+0]*w[11]+a[4*sna_+1]*w[7]-a[4*sna_+2]*w[4]+a[4*sna_+3]*w[2]-a[4*sna_+5]*w[0];
          b[5*snb_+5] = a[4*sna_+0]*w[10]-a[4*sna_+1]*w[6]+a[4*sna_+2]*w[3]-a[4*sna_+3]*w[1]+a[4*sna_+4]*w[0];
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Tiny
} // namespace FEAT

#endif // KERNEL_UTIL_TINY_ALGEBRA_HPP

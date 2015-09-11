#pragma once
#ifndef KERNEL_UTIL_TINY_ALGEBRA_HPP
#define KERNEL_UTIL_TINY_ALGEBRA_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>

namespace FEAST
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
     * Technically, the Tensor3 class template realises a l-tuple of m-by-n matrices.
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
       * This is the end of the recursion where the type is something other than a FEAST::Tiny::... object.
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
      template<typename T_, int m_, int n_, int sm_, int sn_>

      // Same for Matrix
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

      // forward declarations
      template<int m_, int n_>
      struct DetHelper;

      template<int m_, int n_>
      struct VolHelper;

      template<int m_, int n_>
      struct InverseHelper;
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
      Vector()
      {
      }

      /// \brief value-assignment constructor
      explicit Vector(DataType value)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = value;
        }
      }

      /// copy constructor
      template<int sx_>
      Vector(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = x.v[i];
        }
      }

      /// value-assignment operator
      Vector& operator=(DataType value)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = value;
        }
        return *this;
      }

      /// copy-assignment operator
      template<int sx_>
      Vector& operator=(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = x.v[i];
        }
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
      T_& operator()(int i)
      {
        ASSERT((i >= 0) && (i < n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator()() */
      const T_& operator()(int i) const
      {
        ASSERT((i >= 0) && (i < n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator()() */
      T_& operator[](int i)
      {
        ASSERT((i >= 0) && (i < n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator[]() */
      const T_& operator[](int i) const
      {
        ASSERT((i >= 0) && (i < n_), "index i out-of-bounds");
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
      Vector<T_, nn_, s_>& size_cast()
      {
        static_assert((nn_ > 0) && (nn_ <= s_), "invalid cast length");
        return reinterpret_cast<Vector<T_, nn_, s_>&>(*this);
      }

      /** \copydoc size_cast() */
      template<int nn_>
      const Vector<T_, nn_, s_>& size_cast() const
      {
        static_assert((nn_ > 0) && (nn_ <= s_), "invalid cast length");
        return reinterpret_cast<const Vector<T_, nn_, s_>&>(*this);
      }

      /// scalar-multiply operator
      Vector& operator*=(T_ alpha)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] *= alpha;
        }
        return *this;
      }

      /// element-wise-multiply operator
      template <int sx_>
      Vector& operator*=(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] *= x.v[i];
        }
        return *this;
      }

      /// vector-add operator
      template<int sx_>
      Vector& operator+=(const Vector<T_, n_, sx_>& x)
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] += x.v[i];
        }
        return *this;
      }

      /// vector-subtract operator
      template<int sx_>
      Vector& operator-=(const Vector<T_, n_, sx_>& x)
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
      void format(DataType alpha = DataType(0))
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = alpha;
        }
      }

      /**
       * \brief Normalises the vector.
       */
      void normalise()
      {
        DataType norm2(this->norm_euclid());
        ASSERT(norm2 > Math::eps<DataType>(), "Trying to normalise a zero vector!");
        norm2 = DataType(1)/norm2;
        (*this) *= norm2;
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
       * The (right) multiplicant vector for the product.
       *
       * \returns \p *this
       */
      template<int m_, int sma_, int sna_, int sx_>
      Vector& set_mat_vec_mult(const Matrix<T_, n_, m_, sma_, sna_>& a, const Vector<T_, m_, sx_>& x)
      {
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
       * \f[ y\top \leftarrow x^\top \cdot A \Longleftrightarrow y \leftarrow A^\top\cdot x \f]
       *
       * \param[in] x
       * The (left) multiplicant vector for the product.
       *
       * \param[in] a
       * The matrix for the product.
       *
       * \returns \p *this
       */
      template<int m_, int sma_, int sna_, int sx_>
      Vector& set_vec_mat_mult(const Vector<T_, m_, sx_>& x, const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
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
       * \brief Computes the euclid norm of the vector.
       *
       * \returns
       * The euclid norm of the vector.
       */
      ValueType norm_euclid() const
      {
        DataType r(DataType(0));
        for(int i(0); i < n_; ++i)
          r += Math::sqr(v[i]);
        return Math::sqrt(r);
      }
    }; // class Vector

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
    inline T_ dot(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      T_ r(0);

      /// \compilerhack Intel C++ 14 loop vectorisation bug
#if defined(FEAST_COMPILER_INTEL) && (FEAST_COMPILER_INTEL < 1500)
      for(int i(0); (i+1) < n_+1; ++i)
#else
      for(int i(0); i < n_; ++i)
#endif
      {
        r += a.v[i] * b.v[i];
      }
      return r;
    }
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
    inline Vector<T_, n_> operator*(typename Vector<T_, n_>::DataType alpha, const Vector<T_, n_, s_>& x)
    {
      return Vector<T_, n_>(x) *= alpha;
    }

    /// scalar-right-multiplication operator
    template<typename T_, int n_, int s_>
    inline Vector<T_, n_> operator*(const Vector<T_, n_, s_>& x, typename Vector<T_, n_>::DataType alpha)
    {
      return Vector<T_, n_>(x) *= alpha;
    }

    /// vector element-wise-product operator
    template<typename T_, int n_, int sa_, int sb_>
    inline Vector<T_, n_> component_product(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      return Vector<T_, n_>(a) *= b;
    }

    /// vector addition operator
    template<typename T_, int n_, int sa_, int sb_>
    inline Vector<T_, n_> operator+(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      return Vector<T_, n_>(a) += b;
    }

    /// vector subtraction operator
    template<typename T_, int n_, int sa_, int sb_>
    inline Vector<T_, n_> operator-(const Vector<T_, n_, sa_>& a, const Vector<T_, n_, sb_>& b)
    {
      return Vector<T_, n_>(a) -= b;
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
      Matrix()
      {
      }

      /// value-assignment constructor
      explicit Matrix(DataType value)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] = value;
        }
      }

      /// copy constructor
      template<int sma_, int sna_>
      Matrix(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] = a.v[i];
        }
      }

      /// value-assignment operator
      Matrix& operator=(DataType value)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] = value;
        }
        return *this;
      }

      /// assignment operator
      template<int sma_, int sna_>
      Matrix& operator=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] = a.v[i];
        }
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
      T_& operator()(int i, int j)
      {
        ASSERT( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERT( (j >= 0) && (j < n_), "index j out-of-bounds");
        return v[i][j];
      }

      /** \copydoc operator()() */
      const T_& operator()(int i, int j) const
      {
        ASSERT( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERT( (j >= 0) && (j < n_), "index j out-of-bounds");
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
      RowType& operator[](int i)
      {
        ASSERT( (i >= 0) && (i <m_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator[]() */
      const RowType& operator[](int i) const
      {
        ASSERT( (i >= 0) && (i <m_), "index i out-of-bounds");
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
      Matrix<T_, mm_, nn_, sm_, sn_>& size_cast()
      {
        static_assert((mm_ > 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ > 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<Matrix<T_, mm_, nn_, sm_, sn_>&>(*this);
      }

      /** \copydoc size_cast() */
      template<int mm_, int nn_>
      const Matrix<T_, mm_, nn_, sm_, sn_>& size_cast() const
      {
        static_assert((mm_ > 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ > 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<const Matrix<T_, mm_, nn_, sm_, sn_>&>(*this);
      }

      /// scalar-right-multiply-by operator
      Matrix& operator*=(DataType alpha)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] *= alpha;
        }
        return *this;
      }

      /// matrix component-wise addition operator
      template<int sma_, int sna_>
      Matrix& operator+=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          v[i] += a.v[i];
        }
        return *this;
      }

      /// matrix component-wise subtraction operator
      template<int sma_, int sna_>
      Matrix& operator-=(const Matrix<T_, m_, n_, sma_, sna_>& a)
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
      void format(DataType alpha = DataType(0))
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
      DataType hessian_sqr_norm() const
      {
        DataType r(0);
        for(int i(0); i < m_; ++i)
        {
          r += Math::sqr(v[i][i]);
          for(int j(0); j < n_; ++j)
          {
            r += Math::sqr(v[i][j]);
          }
        }
        return r / DataType(2);
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
      DataType trace() const
      {
        int k = (m_ < n_ ? m_ : n_);
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
      DataType det() const
      {
        // Note: We need to specify the template arguments for the 'compute' function explicitly as the
        //       Intel C++ compiler cannot deduct the arguments automatically due to a compiler bug.
        typedef DataType Tv[sm_][sn_];
        return Intern::DetHelper<m_, n_>::template compute<DataType, sm_, sn_>(*reinterpret_cast<const Tv*>(this));
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
      DataType vol() const
      {
        // Note: We need to specify the template arguments for the 'compute' function explicitly as the
        //       Intel C++ compiler cannot deduct the arguments automatically due to a compiler bug.
        typedef DataType Tv[sm_][sn_];
        return Intern::VolHelper<m_, n_>::template compute<DataType, sm_, sn_>(*reinterpret_cast<const Tv*>(this));
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
      Matrix& set_inverse(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        // Note: We need to specify the template arguments for the 'compute' function explicitly as the
        //       Intel C++ compiler cannot deduct the arguments automatically due to a compiler bug.
        typedef DataType Tv[sm_][sn_];
        typedef DataType Ta[sma_][sna_];
        Intern::InverseHelper<m_,n_>::template compute<DataType, sm_, sn_, sma_, sna_>(
          *reinterpret_cast<Tv*>(this), *reinterpret_cast<const Ta*>(&a));
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
      Matrix& set_transpose(const Matrix<T_, n_, m_, sma_, sna_>& a)
      {
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
       * \brief Computes the scalar product of two vectors with this matrix.
       *
       * This function returns
       * \f[ x^\top\cdot A\cdot y = \sum_{i=0}^{m-1}\sum_{j=0}^{n-1} x_i\cdot A_{ij}\cdot y_j\f]
       *
       * \param[in] x
       * The left muliplicant vector of size \p m_.
       *
       * \param[in] y
       * The right multiplicant vector of size \p n_.
       *
       * \returns
       * The scalar product of \p x and \p y with this matrix.
       */
      template<int snx_, int sny_>
      DataType scalar_product(const Vector<T_, m_, snx_>& x, const Vector<T_, n_, sny_>& y) const
      {
        DataType r(DataType(0));
        for(int i(0); i < m_; ++i)
        {
          r += x[i] * dot(v[i], y);
        }
        return r;
      }

      /**
       * \brief Adds the algebraic matrix-product of two other matrices onto this matrix.
       *
       * Let \e C denote \c this matrix, and let \e A denote the left m-by-l matrix and \e B the right l-by-n matrix,
       * then this function computes:
       * \f[ C\leftarrow C + \alpha A\cdot B \f]
       *
       * \param[in] a
       * The left m-by-l multiplicant matrix.
       *
       * \param[in] b
       * The right l-by-n multiplicant matrix.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int l_, int sma_, int sna_, int smb_, int snb_>
      Matrix& add_mat_mat_mult(
        const Matrix<T_, m_, l_, sma_, sna_>& a,
        const Matrix<T_, l_, n_, smb_, snb_>& b,
        DataType alpha = DataType(1))
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            DataType r(0);
            for(int k(0); k < l_; ++k)
            {
              r += a.v[i][k] * b.v[k][j];
            }
            operator()(i,j) += alpha * r;
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
       * The left m-by-l multiplicant matrix.
       *
       * \param[in] b
       * The right l-by-n multiplicant matrix.
       *
       * \returns \p *this
       */
      template<int l_, int sma_, int sna_, int smb_, int snb_>
      Matrix& set_mat_mat_mult(const Matrix<T_, m_, l_, sma_, sna_>& a, const Matrix<T_, l_, n_, smb_, snb_>& b)
      {
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
       * The inner k-by-l multiplicant matrix.
       *
       * \param[in] b
       * The left k-by-m multiplicant matrix.
       *
       * \param[in] d
       * The right l-by-n multiplicant matrix.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int k_, int l_, int sma_, int sna_, int smb_, int snb_, int smd_, int snd_>
      Matrix& add_double_mat_mult(
        const Matrix<T_, k_, l_, sma_, sna_>& a,
        const Matrix<T_, k_, m_, smb_, snb_>& b,
        const Matrix<T_, l_, n_, smd_, snd_>& d,
        DataType alpha = DataType(1))
      {
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
            operator()(i,j) += alpha * r;
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
       * The inner k-by-l multiplicant matrix.
       *
       * \param[in] b
       * The left k-by-m multiplicant matrix.
       *
       * \param[in] d
       * The right l-by-n multiplicant matrix.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int k_, int l_, int sma_, int sna_, int smb_, int snb_, int smd_, int snd_>
      Matrix& set_double_mat_mult(
        const Matrix<T_, k_, l_, sma_, sna_>& a,
        const Matrix<T_, k_, m_, smb_, snb_>& b,
        const Matrix<T_, l_, n_, smd_, snd_>& d,
        T_ alpha = T_(1))
      {
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
       * The l-size vector that serves as a left multiplicant.
       *
       * \param[in] t
       * The l-by-m-by-n tensor that serves as a right multiplicant.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int l_, int snv_, int slt_, int smt_, int snt_>
      Matrix& add_vec_tensor_mult(
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
            operator()(i,j) += alpha * r;
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
       * The l-size vector that serves as a left multiplicant.
       *
       * \param[in] t
       * The l-by-m-by-n tensor that serves as a right multiplicant.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int l_, int snv_, int slt_, int smt_, int snt_>
      Matrix& set_vec_tensor_mult(
        const Vector<T_, l_, snv_>& x,
        const Tensor3<T_, l_, m_, n_, slt_, smt_, snt_>& t,
        DataType alpha = DataType(1))
      {
        format();
        return add_vec_tensor_mult(x, t, alpha);
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
     * A vector that is orthogonal to the m_-1 columns of the input matrix, but not normalised.
     *
     * \note So far, this is only implemented for m_ = 2,3.
     */
    template<typename T_, int m_, int sm_, int sn_>
    Vector<T_, m_> orthogonal(const Matrix<T_, m_, m_-1, sm_, sn_>& tau);
#endif

    /// \cond internal
    template<typename T_, int sm_, int sn_>
    Vector<T_, 2> orthogonal(const Matrix<T_, 2, 1, sm_, sn_>& tau)
    {
      Vector<T_, 2, sm_> nu(T_(0));

      // 2d "cross" product. The sign has to be on the second component so the input is rotated in negative direction
      nu[0] =  tau[1][0];
      nu[1] = -tau[0][0];

      return nu;
    }

    template<typename T_, int sm_, int sn_>
    Vector<T_, 3> orthogonal(const Matrix<T_, 3, 2, sm_, sn_>& tau)
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
    inline Vector<T_, m_> operator*(const Matrix<T_, m_, n_, sm_, sn_>& a, const Vector<T_, n_, sx_>& x)
    {
      return Vector<T_, m_>().set_mat_vec_mult(a, x);
    }

    /// vector-matrix-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_, int sx_>
    inline Vector<T_, n_> operator*(const Vector<T_, m_, sx_>& x, const Matrix<T_, m_, n_, sm_, sn_>& a)
    {
      return Vector<T_, n_>().set_vec_mat_mult(x, a);
    }

    /// scalar-left-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_>
    inline Matrix<T_, m_, n_> operator*(typename Matrix<T_, m_, n_>::DataType alpha, const Matrix<T_, m_, n_, sm_, sn_>& a)
    {
      return Matrix<T_, m_, n_>(a) *= alpha;
    }

    /// scalar-right-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_>
    inline Matrix<T_, m_, n_, sm_, sn_> operator*(const Matrix<T_, m_, n_, sm_, sn_>& a, typename Matrix<T_, m_, n_>::DataType alpha)
    {
      return Matrix<T_, m_, n_>(a) *= alpha;
    }

    /// algebraic matrix-matrix-multiply operator
    template<typename T_, int m_, int n_, int l_, int sma_, int sna_, int smb_, int snb_>
    inline Matrix<T_, m_, n_> operator*(const Matrix<T_, m_, l_, sma_, sna_>& a, const Matrix<T_, l_, n_, smb_, snb_>& b)
    {
      return Matrix<T_, m_, n_>().set_mat_mat_mult(a, b);
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
    inline T_ dot(const Matrix<T_, m_, n_, sma_, sna_>& a, const Matrix<T_, m_, n_, smb_, snb_>& b)
    {
      T_ r(0);
      for(int i(0); i < m_; ++i)
      {
        for(int j(0); j < n_; ++j)
          r += a(i,j) * b(i,j);
      }
      return r;
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
      Tensor3()
      {
      }

      /// value-assignment constructor
      explicit Tensor3(DataType value)
      {
        for(int i(0); i < l_; ++i)
          v[i] = value;
      }

      /// copy-constructor
      template<int sla_, int sma_, int sna_>
      Tensor3(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] = a.v[i];
      }

      /// value-assignment operator
      Tensor3& operator=(DataType value)
      {
        for(int i(0); i < l_; ++i)
          v[i] = value;
        return *this;
      }

      /// copy-assignment operator
      template<int sla_, int sma_, int sna_>
      Tensor3& operator=(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] = a.v[i];
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
      T_& operator()(int h, int i, int j)
      {
        ASSERT( (h >= 0) && (h < l_), "index h out-of-bounds");
        ASSERT( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERT( (j >= 0) && (j < n_), "index j out-of-bounds");
        return v[h](i,j);
      }

      /** \copydoc operator()() */
      const T_& operator()(int h, int i, int j) const
      {
        ASSERT( (h >= 0) && (h < l_), "index h out-of-bounds");
        ASSERT( (i >= 0) && (i < m_), "index i out-of-bounds");
        ASSERT( (j >= 0) && (j < n_), "index j out-of-bounds");
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
      PlaneType& operator[](int h)
      {
        ASSERT( (h >= 0) && (h < l_), "index h out-of-bounds");
        return v[h];
      }

      /** \copydoc operator[]() */
      const PlaneType& operator[](int h) const
      {
        ASSERT( (h >= 0) && (h < l_), "index h out-of-bounds");
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
      Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>& size_cast()
      {
        static_assert((ll_ >= 0) && (ll_ <= sl_), "invalid cast tube count");
        static_assert((mm_ >= 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ >= 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>&>(*this);
      }

      /** \copydoc size_cast() */
      template<int ll_, int mm_, int nn_>
      const Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>& size_cast() const
      {
        static_assert((ll_ >= 0) && (ll_ <= sl_), "invalid cast tube count");
        static_assert((mm_ >= 0) && (mm_ <= sm_), "invalid cast row count");
        static_assert((nn_ >= 0) && (nn_ <= sn_), "invalid cast column count");
        return reinterpret_cast<Tensor3<T_, ll_, mm_, nn_, sl_, sm_, sn_>&>(*this);
      }

      /// scalar right-multiply-by operator
      Tensor3& operator*=(DataType alpha)
      {
        for(int i(0); i < l_; ++i)
          v[i] *= alpha;
        return *this;
      }

      /// tensor component-wise addition operator
      template<int sla_, int sma_, int sna_>
      Tensor3& operator+=(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] += a.v[i];
        return *this;
      }

      /// tensor component-wise subtraction operator
      template<int sla_, int sma_, int sna_>
      Tensor3& operator-=(const Tensor3<T_, l_, m_, n_, sla_, sma_, sna_>& a)
      {
        for(int i(0); i < l_; ++i)
          v[i] -= a.v[i];
        return *this;
      }

      /// formats the tensor
      void format(DataType alpha = DataType(0))
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
       * The l-by-k matrix that serves as the left multiplicant.
       *
       * \param[in] t
       * The k-by-m-by-n tensor that serves as the right multiplicant.
       *
       * \param[in] alpha
       * A scaling factor for the product.
       *
       * \returns \p *this
       */
      template<int k_, int sma_, int sna_, int slt_, int smt_, int snt_>
      Tensor3& add_mat_tensor_mult(
        const Matrix<T_, l_, k_, sma_, sna_>& a,
        const Tensor3<T_, k_, m_, n_, slt_, smt_, snt_>& t,
        DataType alpha = DataType(1))
      {
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
       * The l-by-m'-by-n' tensor that serves as the inner multiplicant.
       *
       * \param[in] b
       * The m'-by-m matrix that serves as the left multiplicant.
       *
       * \param[in] d
       * The n'-by-n matrix that serves as the right multiplicant.
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
      Tensor3& add_double_mat_mult(
        const Tensor3<T_, lt_, mt_, nt_, slt_, smt_, snt_>& t,
        const Matrix<T_, nt_, n_, smb_, snb_>& b,
        const Matrix<T_, mt_, m_, smd_, snd_>& d,
        DataType alpha = DataType(1))
      {
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
    }; // class Tensor3<...>

    /// scalar left-multiply operator
    template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
    inline Tensor3<T_, l_, m_, n_, sl_, sm_, sn_> operator*(
      typename Tensor3<T_, l_, m_, n_>::DataType alpha, const Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>& a)
    {
      return Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>(a) *= alpha;
    }

    /// scalar right-multiply operator
    template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
    inline Tensor3<T_, l_, m_, n_, sl_, sm_, sn_> operator*(
      const Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>& a, typename Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>::DataType alpha)
    {
      return Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>(a) *= alpha;
    }

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    // Internal helpers implementation
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /// \cond internal
    namespace Intern
    {
      template<>
      struct DetHelper<1,1>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          return a[0][0];
        }
      };

      template<>
      struct DetHelper<2,2>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          return a[0][0]*a[1][1] - a[0][1]*a[1][0];
        }
      };

      template<>
      struct DetHelper<3,3>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          return  a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
                + a[0][1]*(a[1][2]*a[2][0] - a[1][0]*a[2][2])
                + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);
        }
      };

      template<>
      struct DetHelper<4,4>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // 2x2 determinants of rows 3-4
          T_ w[6] =
          {
            a[2][0]*a[3][1] - a[2][1]*a[3][0],
            a[2][0]*a[3][2] - a[2][2]*a[3][0],
            a[2][0]*a[3][3] - a[2][3]*a[3][0],
            a[2][1]*a[3][2] - a[2][2]*a[3][1],
            a[2][1]*a[3][3] - a[2][3]*a[3][1],
            a[2][2]*a[3][3] - a[2][3]*a[3][2]
          };

          return
            + a[0][0] * (a[1][1]*w[5] - a[1][2]*w[4] + a[1][3]*w[3])
            - a[0][1] * (a[1][0]*w[5] - a[1][2]*w[2] + a[1][3]*w[1])
            + a[0][2] * (a[1][0]*w[4] - a[1][1]*w[2] + a[1][3]*w[0])
            - a[0][3] * (a[1][0]*w[3] - a[1][1]*w[1] + a[1][2]*w[0]);
        }
      };

      template<>
      struct DetHelper<5,5>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // 2x2 determinants of rows 4-5
          T_ v[10] =
          {
            a[3][0]*a[4][1] - a[3][1]*a[4][0],
            a[3][0]*a[4][2] - a[3][2]*a[4][0],
            a[3][0]*a[4][3] - a[3][3]*a[4][0],
            a[3][0]*a[4][4] - a[3][4]*a[4][0],
            a[3][1]*a[4][2] - a[3][2]*a[4][1],
            a[3][1]*a[4][3] - a[3][3]*a[4][1],
            a[3][1]*a[4][4] - a[3][4]*a[4][1],
            a[3][2]*a[4][3] - a[3][3]*a[4][2],
            a[3][2]*a[4][4] - a[3][4]*a[4][2],
            a[3][3]*a[4][4] - a[3][4]*a[4][3]
          };
          // 3x3 determinants of rows 3-4-5
          T_ w[10] =
          {
            a[2][0]*v[4] - a[2][1]*v[1] + a[2][2]*v[0],
            a[2][0]*v[5] - a[2][1]*v[2] + a[2][3]*v[0],
            a[2][0]*v[6] - a[2][1]*v[3] + a[2][4]*v[0],
            a[2][0]*v[7] - a[2][2]*v[2] + a[2][3]*v[1],
            a[2][0]*v[8] - a[2][2]*v[3] + a[2][4]*v[1],
            a[2][0]*v[9] - a[2][3]*v[3] + a[2][4]*v[2],
            a[2][1]*v[7] - a[2][2]*v[5] + a[2][3]*v[4],
            a[2][1]*v[8] - a[2][2]*v[6] + a[2][4]*v[4],
            a[2][1]*v[9] - a[2][3]*v[6] + a[2][4]*v[5],
            a[2][2]*v[9] - a[2][3]*v[8] + a[2][4]*v[7]
          };

          return
            + a[0][0]*(a[1][1]*w[9] - a[1][2]*w[8] + a[1][3]*w[7] - a[1][4]*w[6])
            - a[0][1]*(a[1][0]*w[9] - a[1][2]*w[5] + a[1][3]*w[4] - a[1][4]*w[3])
            + a[0][2]*(a[1][0]*w[8] - a[1][1]*w[5] + a[1][3]*w[2] - a[1][4]*w[1])
            - a[0][3]*(a[1][0]*w[7] - a[1][1]*w[4] + a[1][2]*w[2] - a[1][4]*w[0])
            + a[0][4]*(a[1][0]*w[6] - a[1][1]*w[3] + a[1][2]*w[1] - a[1][3]*w[0]);
        }
      };

      template<>
      struct DetHelper<6,6>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // 2x2 determinants of rows 5-6
          T_ v[15] =
          {
            a[4][0]*a[5][1] - a[4][1]*a[5][0],
            a[4][0]*a[5][2] - a[4][2]*a[5][0],
            a[4][0]*a[5][3] - a[4][3]*a[5][0],
            a[4][0]*a[5][4] - a[4][4]*a[5][0],
            a[4][0]*a[5][5] - a[4][5]*a[5][0],
            a[4][1]*a[5][2] - a[4][2]*a[5][1],
            a[4][1]*a[5][3] - a[4][3]*a[5][1],
            a[4][1]*a[5][4] - a[4][4]*a[5][1],
            a[4][1]*a[5][5] - a[4][5]*a[5][1],
            a[4][2]*a[5][3] - a[4][3]*a[5][2],
            a[4][2]*a[5][4] - a[4][4]*a[5][2],
            a[4][2]*a[5][5] - a[4][5]*a[5][2],
            a[4][3]*a[5][4] - a[4][4]*a[5][3],
            a[4][3]*a[5][5] - a[4][5]*a[5][3],
            a[4][4]*a[5][5] - a[4][5]*a[5][4]
          };
          // 3x3 determinants of rows 4-5-6
          T_ w[20] =
          {
            a[3][0]*v[ 5] - a[3][1]*v[ 1] + a[3][2]*v[ 0],
            a[3][0]*v[ 6] - a[3][1]*v[ 2] + a[3][3]*v[ 0],
            a[3][0]*v[ 7] - a[3][1]*v[ 3] + a[3][4]*v[ 0],
            a[3][0]*v[ 8] - a[3][1]*v[ 4] + a[3][5]*v[ 0],
            a[3][0]*v[ 9] - a[3][2]*v[ 2] + a[3][3]*v[ 1],
            a[3][0]*v[10] - a[3][2]*v[ 3] + a[3][4]*v[ 1],
            a[3][0]*v[11] - a[3][2]*v[ 4] + a[3][5]*v[ 1],
            a[3][0]*v[12] - a[3][3]*v[ 3] + a[3][4]*v[ 2],
            a[3][0]*v[13] - a[3][3]*v[ 4] + a[3][5]*v[ 2],
            a[3][0]*v[14] - a[3][4]*v[ 4] + a[3][5]*v[ 3],
            a[3][1]*v[ 9] - a[3][2]*v[ 6] + a[3][3]*v[ 5],
            a[3][1]*v[10] - a[3][2]*v[ 7] + a[3][4]*v[ 5],
            a[3][1]*v[11] - a[3][2]*v[ 8] + a[3][5]*v[ 5],
            a[3][1]*v[12] - a[3][3]*v[ 7] + a[3][4]*v[ 6],
            a[3][1]*v[13] - a[3][3]*v[ 8] + a[3][5]*v[ 6],
            a[3][1]*v[14] - a[3][4]*v[ 8] + a[3][5]*v[ 7],
            a[3][2]*v[12] - a[3][3]*v[10] + a[3][4]*v[ 9],
            a[3][2]*v[13] - a[3][3]*v[11] + a[3][5]*v[ 9],
            a[3][2]*v[14] - a[3][4]*v[11] + a[3][5]*v[10],
            a[3][3]*v[14] - a[3][4]*v[13] + a[3][5]*v[12]
          };
          // 4x4 determinants of rows 3-4-5-6
          v[ 0] = a[2][0]*w[10] - a[2][1]*w[ 4] + a[2][2]*w[ 1] - a[2][3]*w[ 0];
          v[ 1] = a[2][0]*w[11] - a[2][1]*w[ 5] + a[2][2]*w[ 2] - a[2][4]*w[ 0];
          v[ 2] = a[2][0]*w[12] - a[2][1]*w[ 6] + a[2][2]*w[ 3] - a[2][5]*w[ 0];
          v[ 3] = a[2][0]*w[13] - a[2][1]*w[ 7] + a[2][3]*w[ 2] - a[2][4]*w[ 1];
          v[ 4] = a[2][0]*w[14] - a[2][1]*w[ 8] + a[2][3]*w[ 3] - a[2][5]*w[ 1];
          v[ 5] = a[2][0]*w[15] - a[2][1]*w[ 9] + a[2][4]*w[ 3] - a[2][5]*w[ 2];
          v[ 6] = a[2][0]*w[16] - a[2][2]*w[ 7] + a[2][3]*w[ 5] - a[2][4]*w[ 4];
          v[ 7] = a[2][0]*w[17] - a[2][2]*w[ 8] + a[2][3]*w[ 6] - a[2][5]*w[ 4];
          v[ 8] = a[2][0]*w[18] - a[2][2]*w[ 9] + a[2][4]*w[ 6] - a[2][5]*w[ 5];
          v[ 9] = a[2][0]*w[19] - a[2][3]*w[ 9] + a[2][4]*w[ 8] - a[2][5]*w[ 7];
          v[10] = a[2][1]*w[16] - a[2][2]*w[13] + a[2][3]*w[11] - a[2][4]*w[10];
          v[11] = a[2][1]*w[17] - a[2][2]*w[14] + a[2][3]*w[12] - a[2][5]*w[10];
          v[12] = a[2][1]*w[18] - a[2][2]*w[15] + a[2][4]*w[12] - a[2][5]*w[11];
          v[13] = a[2][1]*w[19] - a[2][3]*w[15] + a[2][4]*w[14] - a[2][5]*w[13];
          v[14] = a[2][2]*w[19] - a[2][3]*w[18] + a[2][4]*w[17] - a[2][5]*w[16];

          return
            + a[0][0]*(a[1][1]*v[14] - a[1][2]*v[13] + a[1][3]*v[12] - a[1][4]*v[11] + a[1][5]*v[10])
            - a[0][1]*(a[1][0]*v[14] - a[1][2]*v[ 9] + a[1][3]*v[ 8] - a[1][4]*v[ 7] + a[1][5]*v[ 6])
            + a[0][2]*(a[1][0]*v[13] - a[1][1]*v[ 9] + a[1][3]*v[ 5] - a[1][4]*v[ 4] + a[1][5]*v[ 3])
            - a[0][3]*(a[1][0]*v[12] - a[1][1]*v[ 8] + a[1][2]*v[ 5] - a[1][4]*v[ 2] + a[1][5]*v[ 1])
            + a[0][4]*(a[1][0]*v[11] - a[1][1]*v[ 7] + a[1][2]*v[ 4] - a[1][3]*v[ 2] + a[1][5]*v[ 0])
            - a[0][5]*(a[1][0]*v[10] - a[1][1]*v[ 6] + a[1][2]*v[ 3] - a[1][3]*v[ 1] + a[1][4]*v[ 0]);
        }
      };

      template<int m_, int n_>
      struct VolHelper
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // generic fallback implementation: compute b := a^T * a and return sqrt(det(b))
          T_ b[n_][n_];
          for(int i(0); i < n_; ++i)
          {
            for(int j(0); j < n_; ++j)
            {
              b[i][j] = T_(0);
              for(int k(0); k < m_; ++k)
              {
                b[i][j] += a[k][i]*a[k][j];
              }
            }
          }
          return Math::sqrt(DetHelper<n_,n_>::template compute<T_,n_,n_>(b));
        }
      };

      template<int n_>
      struct VolHelper<n_,n_>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // square matrix special case: vol(a) = abs(det(a))
          return Math::abs(DetHelper<n_,n_>::template compute<T_,sm_,sn_>(a));
        }
      };

      template<>
      struct VolHelper<2,1>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // This is the euclid norm of the only matrix column.
          return Math::sqrt(Math::sqr(a[0][0]) + Math::sqr(a[1][0]));
        }
      };

      template<>
      struct VolHelper<3,1>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // This is the euclid norm of the only matrix column.
          return Math::sqrt(Math::sqr(a[0][0]) + Math::sqr(a[1][0]) + Math::sqr(a[2][0]));
        }
      };

      template<>
      struct VolHelper<3,2>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // This is the euclid norm of the 3D cross product of the two matrix columns.
          return Math::sqrt(
            Math::sqr(a[1][0]*a[2][1] - a[2][0]*a[1][1]) +
            Math::sqr(a[2][0]*a[0][1] - a[0][0]*a[2][1]) +
            Math::sqr(a[0][0]*a[1][1] - a[1][0]*a[0][1]));
        }
      };

      template<>
      struct InverseHelper<1,1>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          b[0][0] = T_(1) / a[0][0];
        }
      };

      template<>
      struct InverseHelper<2,2>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          T_ d = T_(1) / (a[0][0]*a[1][1] - a[0][1]*a[1][0]);
          b[0][0] =  d*a[1][1];
          b[0][1] = -d*a[0][1];
          b[1][0] = -d*a[1][0];
          b[1][1] =  d*a[0][0];
        }
      };

      template<>
      struct InverseHelper<3,3>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          b[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
          b[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
          b[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
          T_ d = T_(1) / (a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0]);
          b[0][0] *= d;
          b[1][0] *= d;
          b[2][0] *= d;
          b[0][1] = d*(a[0][2]*a[2][1] - a[0][1]*a[2][2]);
          b[1][1] = d*(a[0][0]*a[2][2] - a[0][2]*a[2][0]);
          b[2][1] = d*(a[0][1]*a[2][0] - a[0][0]*a[2][1]);
          b[0][2] = d*(a[0][1]*a[1][2] - a[0][2]*a[1][1]);
          b[1][2] = d*(a[0][2]*a[1][0] - a[0][0]*a[1][2]);
          b[2][2] = d*(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
        }
      };

      template<>
      struct InverseHelper<4,4>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          T_ w[6];
          w[0] = a[2][0]*a[3][1]-a[2][1]*a[3][0];
          w[1] = a[2][0]*a[3][2]-a[2][2]*a[3][0];
          w[2] = a[2][0]*a[3][3]-a[2][3]*a[3][0];
          w[3] = a[2][1]*a[3][2]-a[2][2]*a[3][1];
          w[4] = a[2][1]*a[3][3]-a[2][3]*a[3][1];
          w[5] = a[2][2]*a[3][3]-a[2][3]*a[3][2];
          b[0][0] = a[1][1]*w[5]-a[1][2]*w[4]+a[1][3]*w[3];
          b[1][0] =-a[1][0]*w[5]+a[1][2]*w[2]-a[1][3]*w[1];
          b[2][0] = a[1][0]*w[4]-a[1][1]*w[2]+a[1][3]*w[0];
          b[3][0] =-a[1][0]*w[3]+a[1][1]*w[1]-a[1][2]*w[0];
          T_ d = T_(1) / (a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0]+a[0][3]*b[3][0]);
          b[0][0] *= d;
          b[1][0] *= d;
          b[2][0] *= d;
          b[3][0] *= d;
          b[0][1] = d*(-a[0][1]*w[5]+a[0][2]*w[4]-a[0][3]*w[3]);
          b[1][1] = d*( a[0][0]*w[5]-a[0][2]*w[2]+a[0][3]*w[1]);
          b[2][1] = d*(-a[0][0]*w[4]+a[0][1]*w[2]-a[0][3]*w[0]);
          b[3][1] = d*( a[0][0]*w[3]-a[0][1]*w[1]+a[0][2]*w[0]);
          w[0] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
          w[1] = a[0][0]*a[1][2]-a[0][2]*a[1][0];
          w[2] = a[0][0]*a[1][3]-a[0][3]*a[1][0];
          w[3] = a[0][1]*a[1][2]-a[0][2]*a[1][1];
          w[4] = a[0][1]*a[1][3]-a[0][3]*a[1][1];
          w[5] = a[0][2]*a[1][3]-a[0][3]*a[1][2];
          b[0][2] = d*( a[3][1]*w[5]-a[3][2]*w[4]+a[3][3]*w[3]);
          b[1][2] = d*(-a[3][0]*w[5]+a[3][2]*w[2]-a[3][3]*w[1]);
          b[2][2] = d*( a[3][0]*w[4]-a[3][1]*w[2]+a[3][3]*w[0]);
          b[3][2] = d*(-a[3][0]*w[3]+a[3][1]*w[1]-a[3][2]*w[0]);
          b[0][3] = d*(-a[2][1]*w[5]+a[2][2]*w[4]-a[2][3]*w[3]);
          b[1][3] = d*( a[2][0]*w[5]-a[2][2]*w[2]+a[2][3]*w[1]);
          b[2][3] = d*(-a[2][0]*w[4]+a[2][1]*w[2]-a[2][3]*w[0]);
          b[3][3] = d*( a[2][0]*w[3]-a[2][1]*w[1]+a[2][2]*w[0]);
        }
      };

      template<>
      struct InverseHelper<5,5>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          T_ w[20];
          w[ 0] = a[3][0]*a[4][1]-a[3][1]*a[4][0];
          w[ 1] = a[3][0]*a[4][2]-a[3][2]*a[4][0];
          w[ 2] = a[3][0]*a[4][3]-a[3][3]*a[4][0];
          w[ 3] = a[3][0]*a[4][4]-a[3][4]*a[4][0];
          w[ 4] = a[3][1]*a[4][2]-a[3][2]*a[4][1];
          w[ 5] = a[3][1]*a[4][3]-a[3][3]*a[4][1];
          w[ 6] = a[3][1]*a[4][4]-a[3][4]*a[4][1];
          w[ 7] = a[3][2]*a[4][3]-a[3][3]*a[4][2];
          w[ 8] = a[3][2]*a[4][4]-a[3][4]*a[4][2];
          w[ 9] = a[3][3]*a[4][4]-a[3][4]*a[4][3];
          w[10] = a[2][0]*w[4]-a[2][1]*w[1]+a[2][2]*w[0];
          w[11] = a[2][0]*w[5]-a[2][1]*w[2]+a[2][3]*w[0];
          w[12] = a[2][0]*w[6]-a[2][1]*w[3]+a[2][4]*w[0];
          w[13] = a[2][0]*w[7]-a[2][2]*w[2]+a[2][3]*w[1];
          w[14] = a[2][0]*w[8]-a[2][2]*w[3]+a[2][4]*w[1];
          w[15] = a[2][0]*w[9]-a[2][3]*w[3]+a[2][4]*w[2];
          w[16] = a[2][1]*w[7]-a[2][2]*w[5]+a[2][3]*w[4];
          w[17] = a[2][1]*w[8]-a[2][2]*w[6]+a[2][4]*w[4];
          w[18] = a[2][1]*w[9]-a[2][3]*w[6]+a[2][4]*w[5];
          w[19] = a[2][2]*w[9]-a[2][3]*w[8]+a[2][4]*w[7];
          b[0][0] = a[1][1]*w[19]-a[1][2]*w[18]+a[1][3]*w[17]-a[1][4]*w[16];
          b[1][0] =-a[1][0]*w[19]+a[1][2]*w[15]-a[1][3]*w[14]+a[1][4]*w[13];
          b[2][0] = a[1][0]*w[18]-a[1][1]*w[15]+a[1][3]*w[12]-a[1][4]*w[11];
          b[3][0] =-a[1][0]*w[17]+a[1][1]*w[14]-a[1][2]*w[12]+a[1][4]*w[10];
          b[4][0] = a[1][0]*w[16]-a[1][1]*w[13]+a[1][2]*w[11]-a[1][3]*w[10];
          T_ d = T_(1) / (a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0]+a[0][3]*b[3][0]+a[0][4]*b[4][0]);
          b[0][0] *= d;
          b[1][0] *= d;
          b[2][0] *= d;
          b[3][0] *= d;
          b[4][0] *= d;
          b[0][1] = d*(-a[0][1]*w[19]+a[0][2]*w[18]-a[0][3]*w[17]+a[0][4]*w[16]);
          b[1][1] = d*( a[0][0]*w[19]-a[0][2]*w[15]+a[0][3]*w[14]-a[0][4]*w[13]);
          b[2][1] = d*(-a[0][0]*w[18]+a[0][1]*w[15]-a[0][3]*w[12]+a[0][4]*w[11]);
          b[3][1] = d*( a[0][0]*w[17]-a[0][1]*w[14]+a[0][2]*w[12]-a[0][4]*w[10]);
          b[4][1] = d*(-a[0][0]*w[16]+a[0][1]*w[13]-a[0][2]*w[11]+a[0][3]*w[10]);
          w[10] = a[1][0]*w[4]-a[1][1]*w[1]+a[1][2]*w[0];
          w[11] = a[1][0]*w[5]-a[1][1]*w[2]+a[1][3]*w[0];
          w[12] = a[1][0]*w[6]-a[1][1]*w[3]+a[1][4]*w[0];
          w[13] = a[1][0]*w[7]-a[1][2]*w[2]+a[1][3]*w[1];
          w[14] = a[1][0]*w[8]-a[1][2]*w[3]+a[1][4]*w[1];
          w[15] = a[1][0]*w[9]-a[1][3]*w[3]+a[1][4]*w[2];
          w[16] = a[1][1]*w[7]-a[1][2]*w[5]+a[1][3]*w[4];
          w[17] = a[1][1]*w[8]-a[1][2]*w[6]+a[1][4]*w[4];
          w[18] = a[1][1]*w[9]-a[1][3]*w[6]+a[1][4]*w[5];
          w[19] = a[1][2]*w[9]-a[1][3]*w[8]+a[1][4]*w[7];
          b[0][2] = d*( a[0][1]*w[19]-a[0][2]*w[18]+a[0][3]*w[17]-a[0][4]*w[16]);
          b[1][2] = d*(-a[0][0]*w[19]+a[0][2]*w[15]-a[0][3]*w[14]+a[0][4]*w[13]);
          b[2][2] = d*( a[0][0]*w[18]-a[0][1]*w[15]+a[0][3]*w[12]-a[0][4]*w[11]);
          b[3][2] = d*(-a[0][0]*w[17]+a[0][1]*w[14]-a[0][2]*w[12]+a[0][4]*w[10]);
          b[4][2] = d*( a[0][0]*w[16]-a[0][1]*w[13]+a[0][2]*w[11]-a[0][3]*w[10]);
          w[ 0] = a[0][0]*a[1][1]-a[0][1]*a[1][0];
          w[ 1] = a[0][0]*a[1][2]-a[0][2]*a[1][0];
          w[ 2] = a[0][0]*a[1][3]-a[0][3]*a[1][0];
          w[ 3] = a[0][0]*a[1][4]-a[0][4]*a[1][0];
          w[ 4] = a[0][1]*a[1][2]-a[0][2]*a[1][1];
          w[ 5] = a[0][1]*a[1][3]-a[0][3]*a[1][1];
          w[ 6] = a[0][1]*a[1][4]-a[0][4]*a[1][1];
          w[ 7] = a[0][2]*a[1][3]-a[0][3]*a[1][2];
          w[ 8] = a[0][2]*a[1][4]-a[0][4]*a[1][2];
          w[ 9] = a[0][3]*a[1][4]-a[0][4]*a[1][3];
          w[10] = a[2][0]*w[4]-a[2][1]*w[1]+a[2][2]*w[0];
          w[11] = a[2][0]*w[5]-a[2][1]*w[2]+a[2][3]*w[0];
          w[12] = a[2][0]*w[6]-a[2][1]*w[3]+a[2][4]*w[0];
          w[13] = a[2][0]*w[7]-a[2][2]*w[2]+a[2][3]*w[1];
          w[14] = a[2][0]*w[8]-a[2][2]*w[3]+a[2][4]*w[1];
          w[15] = a[2][0]*w[9]-a[2][3]*w[3]+a[2][4]*w[2];
          w[16] = a[2][1]*w[7]-a[2][2]*w[5]+a[2][3]*w[4];
          w[17] = a[2][1]*w[8]-a[2][2]*w[6]+a[2][4]*w[4];
          w[18] = a[2][1]*w[9]-a[2][3]*w[6]+a[2][4]*w[5];
          w[19] = a[2][2]*w[9]-a[2][3]*w[8]+a[2][4]*w[7];
          b[0][3] = d*( a[4][1]*w[19]-a[4][2]*w[18]+a[4][3]*w[17]-a[4][4]*w[16]);
          b[1][3] = d*(-a[4][0]*w[19]+a[4][2]*w[15]-a[4][3]*w[14]+a[4][4]*w[13]);
          b[2][3] = d*( a[4][0]*w[18]-a[4][1]*w[15]+a[4][3]*w[12]-a[4][4]*w[11]);
          b[3][3] = d*(-a[4][0]*w[17]+a[4][1]*w[14]-a[4][2]*w[12]+a[4][4]*w[10]);
          b[4][3] = d*( a[4][0]*w[16]-a[4][1]*w[13]+a[4][2]*w[11]-a[4][3]*w[10]);
          b[0][4] = d*(-a[3][1]*w[19]+a[3][2]*w[18]-a[3][3]*w[17]+a[3][4]*w[16]);
          b[1][4] = d*( a[3][0]*w[19]-a[3][2]*w[15]+a[3][3]*w[14]-a[3][4]*w[13]);
          b[2][4] = d*(-a[3][0]*w[18]+a[3][1]*w[15]-a[3][3]*w[12]+a[3][4]*w[11]);
          b[3][4] = d*( a[3][0]*w[17]-a[3][1]*w[14]+a[3][2]*w[12]-a[3][4]*w[10]);
          b[4][4] = d*(-a[3][0]*w[16]+a[3][1]*w[13]-a[3][2]*w[11]+a[3][3]*w[10]);
        }
      };

      template<>
      struct InverseHelper<6,6>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          T_ w[35];
          w[ 0] = a[4][0]*a[5][1]-a[4][1]*a[5][0];
          w[ 1] = a[4][0]*a[5][2]-a[4][2]*a[5][0];
          w[ 2] = a[4][0]*a[5][3]-a[4][3]*a[5][0];
          w[ 3] = a[4][0]*a[5][4]-a[4][4]*a[5][0];
          w[ 4] = a[4][0]*a[5][5]-a[4][5]*a[5][0];
          w[ 5] = a[4][1]*a[5][2]-a[4][2]*a[5][1];
          w[ 6] = a[4][1]*a[5][3]-a[4][3]*a[5][1];
          w[ 7] = a[4][1]*a[5][4]-a[4][4]*a[5][1];
          w[ 8] = a[4][1]*a[5][5]-a[4][5]*a[5][1];
          w[ 9] = a[4][2]*a[5][3]-a[4][3]*a[5][2];
          w[10] = a[4][2]*a[5][4]-a[4][4]*a[5][2];
          w[11] = a[4][2]*a[5][5]-a[4][5]*a[5][2];
          w[12] = a[4][3]*a[5][4]-a[4][4]*a[5][3];
          w[13] = a[4][3]*a[5][5]-a[4][5]*a[5][3];
          w[14] = a[4][4]*a[5][5]-a[4][5]*a[5][4];
          w[15] = a[3][0]*w[5]-a[3][1]*w[1]+a[3][2]*w[0];
          w[16] = a[3][0]*w[6]-a[3][1]*w[2]+a[3][3]*w[0];
          w[17] = a[3][0]*w[7]-a[3][1]*w[3]+a[3][4]*w[0];
          w[18] = a[3][0]*w[8]-a[3][1]*w[4]+a[3][5]*w[0];
          w[19] = a[3][0]*w[9]-a[3][2]*w[2]+a[3][3]*w[1];
          w[20] = a[3][0]*w[10]-a[3][2]*w[3]+a[3][4]*w[1];
          w[21] = a[3][0]*w[11]-a[3][2]*w[4]+a[3][5]*w[1];
          w[22] = a[3][0]*w[12]-a[3][3]*w[3]+a[3][4]*w[2];
          w[23] = a[3][0]*w[13]-a[3][3]*w[4]+a[3][5]*w[2];
          w[24] = a[3][0]*w[14]-a[3][4]*w[4]+a[3][5]*w[3];
          w[25] = a[3][1]*w[9]-a[3][2]*w[6]+a[3][3]*w[5];
          w[26] = a[3][1]*w[10]-a[3][2]*w[7]+a[3][4]*w[5];
          w[27] = a[3][1]*w[11]-a[3][2]*w[8]+a[3][5]*w[5];
          w[28] = a[3][1]*w[12]-a[3][3]*w[7]+a[3][4]*w[6];
          w[29] = a[3][1]*w[13]-a[3][3]*w[8]+a[3][5]*w[6];
          w[30] = a[3][1]*w[14]-a[3][4]*w[8]+a[3][5]*w[7];
          w[31] = a[3][2]*w[12]-a[3][3]*w[10]+a[3][4]*w[9];
          w[32] = a[3][2]*w[13]-a[3][3]*w[11]+a[3][5]*w[9];
          w[33] = a[3][2]*w[14]-a[3][4]*w[11]+a[3][5]*w[10];
          w[34] = a[3][3]*w[14]-a[3][4]*w[13]+a[3][5]*w[12];
          w[ 0] = a[2][0]*w[25]-a[2][1]*w[19]+a[2][2]*w[16]-a[2][3]*w[15];
          w[ 1] = a[2][0]*w[26]-a[2][1]*w[20]+a[2][2]*w[17]-a[2][4]*w[15];
          w[ 2] = a[2][0]*w[27]-a[2][1]*w[21]+a[2][2]*w[18]-a[2][5]*w[15];
          w[ 3] = a[2][0]*w[28]-a[2][1]*w[22]+a[2][3]*w[17]-a[2][4]*w[16];
          w[ 4] = a[2][0]*w[29]-a[2][1]*w[23]+a[2][3]*w[18]-a[2][5]*w[16];
          w[ 5] = a[2][0]*w[30]-a[2][1]*w[24]+a[2][4]*w[18]-a[2][5]*w[17];
          w[ 6] = a[2][0]*w[31]-a[2][2]*w[22]+a[2][3]*w[20]-a[2][4]*w[19];
          w[ 7] = a[2][0]*w[32]-a[2][2]*w[23]+a[2][3]*w[21]-a[2][5]*w[19];
          w[ 8] = a[2][0]*w[33]-a[2][2]*w[24]+a[2][4]*w[21]-a[2][5]*w[20];
          w[ 9] = a[2][0]*w[34]-a[2][3]*w[24]+a[2][4]*w[23]-a[2][5]*w[22];
          w[10] = a[2][1]*w[31]-a[2][2]*w[28]+a[2][3]*w[26]-a[2][4]*w[25];
          w[11] = a[2][1]*w[32]-a[2][2]*w[29]+a[2][3]*w[27]-a[2][5]*w[25];
          w[12] = a[2][1]*w[33]-a[2][2]*w[30]+a[2][4]*w[27]-a[2][5]*w[26];
          w[13] = a[2][1]*w[34]-a[2][3]*w[30]+a[2][4]*w[29]-a[2][5]*w[28];
          w[14] = a[2][2]*w[34]-a[2][3]*w[33]+a[2][4]*w[32]-a[2][5]*w[31];
          b[0][0] =  a[1][1]*w[14]-a[1][2]*w[13]+a[1][3]*w[12]-a[1][4]*w[11]+a[1][5]*w[10];
          b[1][0] = -a[1][0]*w[14]+a[1][2]*w[9]-a[1][3]*w[8]+a[1][4]*w[7]-a[1][5]*w[6];
          b[2][0] =  a[1][0]*w[13]-a[1][1]*w[9]+a[1][3]*w[5]-a[1][4]*w[4]+a[1][5]*w[3];
          b[3][0] = -a[1][0]*w[12]+a[1][1]*w[8]-a[1][2]*w[5]+a[1][4]*w[2]-a[1][5]*w[1];
          b[4][0] =  a[1][0]*w[11]-a[1][1]*w[7]+a[1][2]*w[4]-a[1][3]*w[2]+a[1][5]*w[0];
          b[5][0] = -a[1][0]*w[10]+a[1][1]*w[6]-a[1][2]*w[3]+a[1][3]*w[1]-a[1][4]*w[0];
          T_ d = T_(1) / (a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0]
                        + a[0][3]*b[3][0] + a[0][4]*b[4][0] + a[0][5]*b[5][0]);
          b[0][0] *= d;
          b[1][0] *= d;
          b[2][0] *= d;
          b[3][0] *= d;
          b[4][0] *= d;
          b[5][0] *= d;
          b[0][1] = d*(-a[0][1]*w[14]+a[0][2]*w[13]-a[0][3]*w[12]+a[0][4]*w[11]-a[0][5]*w[10]);
          b[1][1] = d*( a[0][0]*w[14]-a[0][2]*w[9]+a[0][3]*w[8]-a[0][4]*w[7]+a[0][5]*w[6]);
          b[2][1] = d*(-a[0][0]*w[13]+a[0][1]*w[9]-a[0][3]*w[5]+a[0][4]*w[4]-a[0][5]*w[3]);
          b[3][1] = d*( a[0][0]*w[12]-a[0][1]*w[8]+a[0][2]*w[5]-a[0][4]*w[2]+a[0][5]*w[1]);
          b[4][1] = d*(-a[0][0]*w[11]+a[0][1]*w[7]-a[0][2]*w[4]+a[0][3]*w[2]-a[0][5]*w[0]);
          b[5][1] = d*( a[0][0]*w[10]-a[0][1]*w[6]+a[0][2]*w[3]-a[0][3]*w[1]+a[0][4]*w[0]);
          w[ 0] = a[4][0]*a[5][1]-a[4][1]*a[5][0];
          w[ 1] = a[4][0]*a[5][2]-a[4][2]*a[5][0];
          w[ 2] = a[4][0]*a[5][3]-a[4][3]*a[5][0];
          w[ 3] = a[4][0]*a[5][4]-a[4][4]*a[5][0];
          w[ 4] = a[4][0]*a[5][5]-a[4][5]*a[5][0];
          w[ 5] = a[4][1]*a[5][2]-a[4][2]*a[5][1];
          w[ 6] = a[4][1]*a[5][3]-a[4][3]*a[5][1];
          w[ 7] = a[4][1]*a[5][4]-a[4][4]*a[5][1];
          w[ 8] = a[4][1]*a[5][5]-a[4][5]*a[5][1];
          w[ 9] = a[4][2]*a[5][3]-a[4][3]*a[5][2];
          w[10] = a[4][2]*a[5][4]-a[4][4]*a[5][2];
          w[11] = a[4][2]*a[5][5]-a[4][5]*a[5][2];
          w[12] = a[4][3]*a[5][4]-a[4][4]*a[5][3];
          w[13] = a[4][3]*a[5][5]-a[4][5]*a[5][3];
          w[14] = a[4][4]*a[5][5]-a[4][5]*a[5][4];
          w[15] = a[1][0]*w[5]-a[1][1]*w[1]+a[1][2]*w[0];
          w[16] = a[1][0]*w[6]-a[1][1]*w[2]+a[1][3]*w[0];
          w[17] = a[1][0]*w[7]-a[1][1]*w[3]+a[1][4]*w[0];
          w[18] = a[1][0]*w[8]-a[1][1]*w[4]+a[1][5]*w[0];
          w[19] = a[1][0]*w[9]-a[1][2]*w[2]+a[1][3]*w[1];
          w[20] = a[1][0]*w[10]-a[1][2]*w[3]+a[1][4]*w[1];
          w[21] = a[1][0]*w[11]-a[1][2]*w[4]+a[1][5]*w[1];
          w[22] = a[1][0]*w[12]-a[1][3]*w[3]+a[1][4]*w[2];
          w[23] = a[1][0]*w[13]-a[1][3]*w[4]+a[1][5]*w[2];
          w[24] = a[1][0]*w[14]-a[1][4]*w[4]+a[1][5]*w[3];
          w[25] = a[1][1]*w[9]-a[1][2]*w[6]+a[1][3]*w[5];
          w[26] = a[1][1]*w[10]-a[1][2]*w[7]+a[1][4]*w[5];
          w[27] = a[1][1]*w[11]-a[1][2]*w[8]+a[1][5]*w[5];
          w[28] = a[1][1]*w[12]-a[1][3]*w[7]+a[1][4]*w[6];
          w[29] = a[1][1]*w[13]-a[1][3]*w[8]+a[1][5]*w[6];
          w[30] = a[1][1]*w[14]-a[1][4]*w[8]+a[1][5]*w[7];
          w[31] = a[1][2]*w[12]-a[1][3]*w[10]+a[1][4]*w[9];
          w[32] = a[1][2]*w[13]-a[1][3]*w[11]+a[1][5]*w[9];
          w[33] = a[1][2]*w[14]-a[1][4]*w[11]+a[1][5]*w[10];
          w[34] = a[1][3]*w[14]-a[1][4]*w[13]+a[1][5]*w[12];
          w[ 0] = a[0][0]*w[25]-a[0][1]*w[19]+a[0][2]*w[16]-a[0][3]*w[15];
          w[ 1] = a[0][0]*w[26]-a[0][1]*w[20]+a[0][2]*w[17]-a[0][4]*w[15];
          w[ 2] = a[0][0]*w[27]-a[0][1]*w[21]+a[0][2]*w[18]-a[0][5]*w[15];
          w[ 3] = a[0][0]*w[28]-a[0][1]*w[22]+a[0][3]*w[17]-a[0][4]*w[16];
          w[ 4] = a[0][0]*w[29]-a[0][1]*w[23]+a[0][3]*w[18]-a[0][5]*w[16];
          w[ 5] = a[0][0]*w[30]-a[0][1]*w[24]+a[0][4]*w[18]-a[0][5]*w[17];
          w[ 6] = a[0][0]*w[31]-a[0][2]*w[22]+a[0][3]*w[20]-a[0][4]*w[19];
          w[ 7] = a[0][0]*w[32]-a[0][2]*w[23]+a[0][3]*w[21]-a[0][5]*w[19];
          w[ 8] = a[0][0]*w[33]-a[0][2]*w[24]+a[0][4]*w[21]-a[0][5]*w[20];
          w[ 9] = a[0][0]*w[34]-a[0][3]*w[24]+a[0][4]*w[23]-a[0][5]*w[22];
          w[10] = a[0][1]*w[31]-a[0][2]*w[28]+a[0][3]*w[26]-a[0][4]*w[25];
          w[11] = a[0][1]*w[32]-a[0][2]*w[29]+a[0][3]*w[27]-a[0][5]*w[25];
          w[12] = a[0][1]*w[33]-a[0][2]*w[30]+a[0][4]*w[27]-a[0][5]*w[26];
          w[13] = a[0][1]*w[34]-a[0][3]*w[30]+a[0][4]*w[29]-a[0][5]*w[28];
          w[14] = a[0][2]*w[34]-a[0][3]*w[33]+a[0][4]*w[32]-a[0][5]*w[31];
          b[0][2] = d*( a[3][1]*w[14]-a[3][2]*w[13]+a[3][3]*w[12]-a[3][4]*w[11]+a[3][5]*w[10]);
          b[1][2] = d*(-a[3][0]*w[14]+a[3][2]*w[9]-a[3][3]*w[8]+a[3][4]*w[7]-a[3][5]*w[6]);
          b[2][2] = d*( a[3][0]*w[13]-a[3][1]*w[9]+a[3][3]*w[5]-a[3][4]*w[4]+a[3][5]*w[3]);
          b[3][2] = d*(-a[3][0]*w[12]+a[3][1]*w[8]-a[3][2]*w[5]+a[3][4]*w[2]-a[3][5]*w[1]);
          b[4][2] = d*( a[3][0]*w[11]-a[3][1]*w[7]+a[3][2]*w[4]-a[3][3]*w[2]+a[3][5]*w[0]);
          b[5][2] = d*(-a[3][0]*w[10]+a[3][1]*w[6]-a[3][2]*w[3]+a[3][3]*w[1]-a[3][4]*w[0]);
          b[0][3] = d*(-a[2][1]*w[14]+a[2][2]*w[13]-a[2][3]*w[12]+a[2][4]*w[11]-a[2][5]*w[10]);
          b[1][3] = d*( a[2][0]*w[14]-a[2][2]*w[9]+a[2][3]*w[8]-a[2][4]*w[7]+a[2][5]*w[6]);
          b[2][3] = d*(-a[2][0]*w[13]+a[2][1]*w[9]-a[2][3]*w[5]+a[2][4]*w[4]-a[2][5]*w[3]);
          b[3][3] = d*( a[2][0]*w[12]-a[2][1]*w[8]+a[2][2]*w[5]-a[2][4]*w[2]+a[2][5]*w[1]);
          b[4][3] = d*(-a[2][0]*w[11]+a[2][1]*w[7]-a[2][2]*w[4]+a[2][3]*w[2]-a[2][5]*w[0]);
          b[5][3] = d*( a[2][0]*w[10]-a[2][1]*w[6]+a[2][2]*w[3]-a[2][3]*w[1]+a[2][4]*w[0]);
          w[ 0] = a[2][0]*a[3][1]-a[2][1]*a[3][0];
          w[ 1] = a[2][0]*a[3][2]-a[2][2]*a[3][0];
          w[ 2] = a[2][0]*a[3][3]-a[2][3]*a[3][0];
          w[ 3] = a[2][0]*a[3][4]-a[2][4]*a[3][0];
          w[ 4] = a[2][0]*a[3][5]-a[2][5]*a[3][0];
          w[ 5] = a[2][1]*a[3][2]-a[2][2]*a[3][1];
          w[ 6] = a[2][1]*a[3][3]-a[2][3]*a[3][1];
          w[ 7] = a[2][1]*a[3][4]-a[2][4]*a[3][1];
          w[ 8] = a[2][1]*a[3][5]-a[2][5]*a[3][1];
          w[ 9] = a[2][2]*a[3][3]-a[2][3]*a[3][2];
          w[10] = a[2][2]*a[3][4]-a[2][4]*a[3][2];
          w[11] = a[2][2]*a[3][5]-a[2][5]*a[3][2];
          w[12] = a[2][3]*a[3][4]-a[2][4]*a[3][3];
          w[13] = a[2][3]*a[3][5]-a[2][5]*a[3][3];
          w[14] = a[2][4]*a[3][5]-a[2][5]*a[3][4];
          w[15] = a[1][0]*w[5]-a[1][1]*w[1]+a[1][2]*w[0];
          w[16] = a[1][0]*w[6]-a[1][1]*w[2]+a[1][3]*w[0];
          w[17] = a[1][0]*w[7]-a[1][1]*w[3]+a[1][4]*w[0];
          w[18] = a[1][0]*w[8]-a[1][1]*w[4]+a[1][5]*w[0];
          w[19] = a[1][0]*w[9]-a[1][2]*w[2]+a[1][3]*w[1];
          w[20] = a[1][0]*w[10]-a[1][2]*w[3]+a[1][4]*w[1];
          w[21] = a[1][0]*w[11]-a[1][2]*w[4]+a[1][5]*w[1];
          w[22] = a[1][0]*w[12]-a[1][3]*w[3]+a[1][4]*w[2];
          w[23] = a[1][0]*w[13]-a[1][3]*w[4]+a[1][5]*w[2];
          w[24] = a[1][0]*w[14]-a[1][4]*w[4]+a[1][5]*w[3];
          w[25] = a[1][1]*w[9]-a[1][2]*w[6]+a[1][3]*w[5];
          w[26] = a[1][1]*w[10]-a[1][2]*w[7]+a[1][4]*w[5];
          w[27] = a[1][1]*w[11]-a[1][2]*w[8]+a[1][5]*w[5];
          w[28] = a[1][1]*w[12]-a[1][3]*w[7]+a[1][4]*w[6];
          w[29] = a[1][1]*w[13]-a[1][3]*w[8]+a[1][5]*w[6];
          w[30] = a[1][1]*w[14]-a[1][4]*w[8]+a[1][5]*w[7];
          w[31] = a[1][2]*w[12]-a[1][3]*w[10]+a[1][4]*w[9];
          w[32] = a[1][2]*w[13]-a[1][3]*w[11]+a[1][5]*w[9];
          w[33] = a[1][2]*w[14]-a[1][4]*w[11]+a[1][5]*w[10];
          w[34] = a[1][3]*w[14]-a[1][4]*w[13]+a[1][5]*w[12];
          w[ 0] = a[0][0]*w[25]-a[0][1]*w[19]+a[0][2]*w[16]-a[0][3]*w[15];
          w[ 1] = a[0][0]*w[26]-a[0][1]*w[20]+a[0][2]*w[17]-a[0][4]*w[15];
          w[ 2] = a[0][0]*w[27]-a[0][1]*w[21]+a[0][2]*w[18]-a[0][5]*w[15];
          w[ 3] = a[0][0]*w[28]-a[0][1]*w[22]+a[0][3]*w[17]-a[0][4]*w[16];
          w[ 4] = a[0][0]*w[29]-a[0][1]*w[23]+a[0][3]*w[18]-a[0][5]*w[16];
          w[ 5] = a[0][0]*w[30]-a[0][1]*w[24]+a[0][4]*w[18]-a[0][5]*w[17];
          w[ 6] = a[0][0]*w[31]-a[0][2]*w[22]+a[0][3]*w[20]-a[0][4]*w[19];
          w[ 7] = a[0][0]*w[32]-a[0][2]*w[23]+a[0][3]*w[21]-a[0][5]*w[19];
          w[ 8] = a[0][0]*w[33]-a[0][2]*w[24]+a[0][4]*w[21]-a[0][5]*w[20];
          w[ 9] = a[0][0]*w[34]-a[0][3]*w[24]+a[0][4]*w[23]-a[0][5]*w[22];
          w[10] = a[0][1]*w[31]-a[0][2]*w[28]+a[0][3]*w[26]-a[0][4]*w[25];
          w[11] = a[0][1]*w[32]-a[0][2]*w[29]+a[0][3]*w[27]-a[0][5]*w[25];
          w[12] = a[0][1]*w[33]-a[0][2]*w[30]+a[0][4]*w[27]-a[0][5]*w[26];
          w[13] = a[0][1]*w[34]-a[0][3]*w[30]+a[0][4]*w[29]-a[0][5]*w[28];
          w[14] = a[0][2]*w[34]-a[0][3]*w[33]+a[0][4]*w[32]-a[0][5]*w[31];
          b[0][4] = d*( a[5][1]*w[14]-a[5][2]*w[13]+a[5][3]*w[12]-a[5][4]*w[11]+a[5][5]*w[10]);
          b[1][4] = d*(-a[5][0]*w[14]+a[5][2]*w[9]-a[5][3]*w[8]+a[5][4]*w[7]-a[5][5]*w[6]);
          b[2][4] = d*( a[5][0]*w[13]-a[5][1]*w[9]+a[5][3]*w[5]-a[5][4]*w[4]+a[5][5]*w[3]);
          b[3][4] = d*(-a[5][0]*w[12]+a[5][1]*w[8]-a[5][2]*w[5]+a[5][4]*w[2]-a[5][5]*w[1]);
          b[4][4] = d*( a[5][0]*w[11]-a[5][1]*w[7]+a[5][2]*w[4]-a[5][3]*w[2]+a[5][5]*w[0]);
          b[5][4] = d*(-a[5][0]*w[10]+a[5][1]*w[6]-a[5][2]*w[3]+a[5][3]*w[1]-a[5][4]*w[0]);
          b[0][5] = d*(-a[4][1]*w[14]+a[4][2]*w[13]-a[4][3]*w[12]+a[4][4]*w[11]-a[4][5]*w[10]);
          b[1][5] = d*( a[4][0]*w[14]-a[4][2]*w[9]+a[4][3]*w[8]-a[4][4]*w[7]+a[4][5]*w[6]);
          b[2][5] = d*(-a[4][0]*w[13]+a[4][1]*w[9]-a[4][3]*w[5]+a[4][4]*w[4]-a[4][5]*w[3]);
          b[3][5] = d*( a[4][0]*w[12]-a[4][1]*w[8]+a[4][2]*w[5]-a[4][4]*w[2]+a[4][5]*w[1]);
          b[4][5] = d*(-a[4][0]*w[11]+a[4][1]*w[7]-a[4][2]*w[4]+a[4][3]*w[2]-a[4][5]*w[0]);
          b[5][5] = d*( a[4][0]*w[10]-a[4][1]*w[6]+a[4][2]*w[3]-a[4][3]*w[1]+a[4][4]*w[0]);
        }
      };

      // generic square matrix inversion:
      // \author Christoph Lohmann
      template<int n_>
      struct InverseHelper<n_, n_>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          // copy matrix a to b
          for (int i(0); i < n_; ++i)
          {
            for (int j(0); j < n_; ++j)
            {
              b[i][j] = a[i][j];
            }
          }

          // create pivot array
          int p[3*n_];

          // perform matrix inversion
          Math::invert_matrix(n_, snb_, &b[0][0], p);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Tiny
} // namespace FEAST

#endif // KERNEL_UTIL_TINY_ALGEBRA_HPP

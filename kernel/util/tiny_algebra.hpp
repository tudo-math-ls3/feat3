#pragma once
#ifndef KERNEL_UTIL_TINY_ALGEBRA_HPP
#define KERNEL_UTIL_TINY_ALGEBRA_HPP 1

// includes, FEAST
#include <kernel/util/assertion.hpp>

// includes, STL
#include <cmath>

namespace FEAST
{
  /// Tiny namespace
  namespace Tiny
  {
    /**
     * \brief Tiny Vector class template
     *
     * This class template implements a vector whose datatype and size is given at compile-time.
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
     * This class template implements a matrix whose datatype and sizes are given at compile-time.
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

    /// \cond internal
    namespace Intern
    {
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
      /// dummy enum
      enum
      {
        /// the length of the vector
        n = n_,
        /// the stride of the vector
        s = s_
      };

      /// the data type of the vector
      typedef T_ DataType;

      /// actual vector data
      T_ v[s_];

      /// default constructor
      Vector()
      {
      }

      /// \brief value-assignment constructor
      explicit Vector(T_ value)
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
      Vector& operator=(T_ value)
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
      T_& operator()(Index i)
      {
        ASSERT(i < Index(n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator()() */
      const T_& operator()(Index i) const
      {
        ASSERT(i < Index(n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator()() */
      T_& operator[](Index i)
      {
        ASSERT(i < Index(n_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator()() */
      const T_& operator[](Index i) const
      {
        ASSERT(i < Index(n_), "index i out-of-bounds");
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
        return this;
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
       * \brief Clears the vector.
       *
       * \param[in] alpha
       * The value that the vector is to be set to.
       */
      void clear(T_ alpha = T_(0))
      {
        for(int i(0); i < n_; ++i)
        {
          v[i] = alpha;
        }
      }

      /**
       * \brief Sets this vector to the result of a matrix-vector product.
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
        for(int i(0); i < m_; ++i)
        {
          v[i] = T_(0);
          for(int j(0); j < n_; ++j)
          {
            v[i] += a.v[i][j] * x.v[j];
          }
        }
        return *this;
      }

      /**
       * \brief Sets this vector to the result of a vector-matrix product.
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
    }; // class Vector

    /**
     * \brief Computes the dot-product of two vectors.
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
      for(int i(0); i < n_; ++i)
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
    inline Vector<T_, n_> operator*(T_ alpha, const Vector<T_, n_, s_>& x)
    {
      return Vector<T_, n_>(x) *= alpha;
    }

    /// scalar-right-multiplication operator
    template<typename T_, int n_, int s_>
    inline Vector<T_, n_> operator*(const Vector<T_, n_, s_>& x, T_ alpha)
    {
      return Vector<T_, n_>(x) *= alpha;
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
      static_assert(n_ > 0, "invalid column length");
      static_assert(sm_ >= m_, "invalid row stride");
      static_assert(sn_ >= n_, "invalid column stride");

    public:
      /// dummy enum
      enum
      {
        /// the row count of the matrix
        m = m_,
        /// the column count of the matrix
        n = n_,
        /// the row stride of the matrix
        sm = sm_,
        /// the column stride of the matrix
        sn = sn_
      };

      /// the data type of the matrix
      typedef T_ DataType;

      /// actual matrix data
      T_ v[sm_][sn_];

      /// default constructor
      Matrix()
      {
      }

      /// value-assignment constructor
      explicit Matrix(T_ value)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j]= value;
          }
        }
      }

      /// copy constructor
      template<int sma_, int sna_>
      Matrix(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] = a.v[i][j];
          }
        }
      }

      /// value-assignment operator
      Matrix& operator=(T_ value)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j]= value;
          }
        }
        return *this;
      }

      /// assignment operator
      template<int sma_, int sna_>
      Matrix& operator=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] = a.v[i][j];
          }
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
      T_& operator()(Index i, Index j)
      {
        ASSERT(i < Index(m_), "index i out-of-bounds");
        ASSERT(j < Index(n_), "index j out-of-bounds");
        return v[i][j];
      }

      /** \copydoc operator()() */
      const T_& operator()(Index i, Index j) const
      {
        ASSERT(i < Index(m_), "index i out-of-bounds");
        ASSERT(j < Index(n_), "index j out-of-bounds");
        return v[i][j];
      }

      /**
       * \brief Row-Access operator.
       *
       * \param[in] i
       * The index of the row that is to be returned.
       *
       * \returns
       * A (const) pointer to the <c>i</c>-th row of the matrix.S
       */
      T_* operator[](Index i)
      {
        ASSERT(i < Index(m_), "index i out-of-bounds");
        return v[i];
      }

      /** \copydoc operator[]() */
      const T_* operator[](Index i) const
      {
        ASSERT(i < Index(m_), "index i out-of-bounds");
        return v[i];
      }

      /**
       * \brief Size-cast function.
       *
       * This function casts this matrix's reference to another size.
       *
       * \tparam nn_, mm_
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
      Matrix& operator*=(T_ alpha)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] *= alpha;
          }
        }
        return *this;
      }

      /// matrix component-wise addition operator
      template<int sma_, int sna_>
      Matrix& operator+=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] += a.v[i][j];
          }
        }
        return *this;
      }

      /// matrix component-wise subtraction operator
      template<int sma_, int sna_>
      Matrix& operator-=(const Matrix<T_, m_, n_, sma_, sna_>& a)
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] -= a.v[i][j];
          }
        }
        return *this;
      }

      /**
       * \brief Clears the matrix.
       *
       * \param[in] alpha
       * The value that the matrix is to be set to.
       */
      void clear(T_ alpha = T_(0))
      {
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] = alpha;
          }
        }
      }

      /**
       * \brief Returns the determinant of the matrix.
       *
       * \warning This function only works for \p m_ = \p n_ and will intentionally fail to compile in any other case.
       *
       * \returns The determinant of the matrix.
       */
      T_ det() const
      {
        // Note: We need to specify the template arguments for the 'compute' function explicitly as the
        //       Intel C++ compiler cannot deduct the arguments automatically due to a compiler bug.
        return Intern::DetHelper<m_, n_>::template compute<T_, sm_, sn_>(v);
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
      T_ vol() const
      {
        // Note: We need to specify the template arguments for the 'compute' function explicitly as the
        //       Intel C++ compiler cannot deduct the arguments automatically due to a compiler bug.
        return Intern::VolHelper<m_, n_>::template compute<T_, sm_, sn_>(v);
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
        Intern::InverseHelper<m_,n_>::template compute<T_, sm_, sn_, sma_, sna_>(v, a.v);
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
       * \brief Sets this matrix to the algebraic matrix-product of two other matrices.
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
        for(int i(0); i < m_; ++i)
        {
          for(int j(0); j < n_; ++j)
          {
            v[i][j] = T_(0);
            for(int k(0); k < l_; ++k)
            {
              v[i][j] += a.v[i][k] * b.v[k][j];
            }
          }
        }
        return *this;
      }
    }; // class Matrix

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
    inline Matrix<T_, m_, n_> operator*(T_ alpha, const Matrix<T_, m_, n_, sm_, sn_>& a)
    {
      return Matrix<T_, m_, n_>(a) *= alpha;
    }

    /// scalar-right-multiply operator
    template<typename T_, int m_, int n_, int sm_, int sn_>
    inline Matrix<T_, m_, n_, sm_, sn_> operator*(const Matrix<T_, m_, n_, sm_, sn_>& a, T_ alpha)
    {
      return Matrix<T_, m_, n_>(a) *= alpha;
    }

    /// algebraic matrix-matrix-multiply operator
    template<typename T_, int m_, int n_, int l_, int sma_, int sna_, int smb_, int snb_>
    inline Matrix<T_, m_, n_> operator*(const Matrix<T_, m_, l_, sma_, sna_>& a, const Matrix<T_, l_, n_, smb_, snb_>& b)
    {
      return Matrix<T_, m_, n_>().set_mat_mat_mult(a, b);
    }

    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */
    // Internal helpers implementation
    /* ************************************************************************************************************* */
    /* ************************************************************************************************************* */

    /// \cond internal
    namespace Intern
    {
      template<typename T_>
      inline T_ sqr(T_ x)
      {
        return x*x;
      }

      template<typename T_>
      inline T_ cub(T_ x)
      {
        return x*x*x;
      }

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
          return std::sqrt(DetHelper<n_,n_>::template compute<T_,n_,n_>(b));
        }
      };

      template<int n_>
      struct VolHelper<n_,n_>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // square matrix special case: vol(a) = abs(det(a))
          return std::abs(DetHelper<n_,n_>::template compute<T_,sm_,sn_>(a));
        }
      };

      template<>
      struct VolHelper<2,1>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // This is the euclid norm of the only matrix column.
          return std::sqrt(sqr(a[0][0]) + sqr(a[1][0]));
        }
      };

      template<>
      struct VolHelper<3,1>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // This is the euclid norm of the only matrix column.
          return std::sqrt(sqr(a[0][0]) + sqr(a[1][0]) + sqr(a[2][0]));
        }
      };

      template<>
      struct VolHelper<3,2>
      {
        template<typename T_, int sm_, int sn_>
        static T_ compute(const T_ (&a)[sm_][sn_])
        {
          // This is the euclid norm of the 3D cross product of the two matrix columns.
          return std::sqrt(
            sqr(a[1][0]*a[2][1] - a[2][0]*a[1][1]) +
            sqr(a[2][0]*a[0][1] - a[0][0]*a[2][1]) +
            sqr(a[0][0]*a[1][1] - a[1][0]*a[0][1]));
        }
      };

      // generic square matrix inversion:
      // currently disabled to avoid dependency on linear_algebra.hpp
      /*template<int n_>
      struct InverseHelper<n_, n_>
      {
        template<typename T_, int smb_, int snb_, int sma_, int sna_>
        static void compute(T_ (&b)[smb_][snb_], const T_ (&a)[sma_][sna_])
        {
          T_ lu[n_][n_];
          int p[n_];
          LinAlg::mat_copy<false>(n_, n_, n_, &lu[0][0], sna_, &a[0][0]);
          LinAlg::mat_identity(n_, snb_, &b[0][0]);
          LinAlg::mat_factorise(n_, n_, n_, &lu[0][0], p);
          LinAlg::mat_solve_mat<false>(n_, n_, snb_, &b[0][0], n_, &lu[0][0], p);
        }
      };*/

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
    } // namespace Intern
    /// \endcond
  } // namespace Tiny
} // namespace FEAST

#endif // KERNEL_UTIL_TINY_ALGEBRA_HPP

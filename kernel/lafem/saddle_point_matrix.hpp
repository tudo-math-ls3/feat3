#pragma once
#ifndef KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP
#define KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP 1

// includes, FEAST
#include <kernel/lafem/tuple_vector.hpp>

// includes, system
#include <type_traits>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Saddle-Point matrix element helper class template
     *
     * This class template is a helper to realise the "at" member function of the SaddlePointMatrix class template.
     *
     * \author Peter Zajac
     */
    template<
      typename MatrixA_,
      typename MatrixB_,
      typename MatrixD_,
      Index i_,
      Index j_>
    struct SaddlePointMatrixElement;

    /**
     * \brief Saddle-Point matrix meta class template
     *
     * This class template implements a meta container for saddle-point matrices:
     * \verbatim
       / A B \
       \ D 0 /
       \endverbatim
     * where
     *  - \e A is an n-by-n (meta) matrix
     *  - \e B is an n-by-m (meta) matrix
     *  - \e D is an m-by-n (meta) matrix
     *
     * \tparam MatrixA_
     * The type of the upper-left (square) matrix block \e A.
     *
     * \tparam MatrixB_
     * The type of the upper-right (rectangular) matrix block \e B.
     *
     * \tparam MatrixD_
     * The type of the lower-left (rectangular) matrix block \e D.
     *
     * \author Peter Zajac
     */
    template<
      typename MatrixA_,
      typename MatrixB_ = MatrixA_,
      typename MatrixD_ = MatrixB_>
    class SaddlePointMatrix
    {
    public:
      /// type of sub-matrix A
      typedef MatrixA_ MatrixTypeA;
      /// type of sub-matrix B
      typedef MatrixB_ MatrixTypeB;
      /// type of sub-matrix D
      typedef MatrixD_ MatrixTypeD;

      // ensure that all matrices have the same mem- and data-types
      static_assert(std::is_same<typename MatrixA_::MemType, typename MatrixB_::MemType>::value,
        "A and B have different mem-types");
      static_assert(std::is_same<typename MatrixA_::MemType, typename MatrixD_::MemType>::value,
        "A and D have different mem-types");
      static_assert(std::is_same<typename MatrixA_::DataType, typename MatrixB_::DataType>::value,
        "A and B have different data-types");
      static_assert(std::is_same<typename MatrixA_::DataType, typename MatrixD_::DataType>::value,
        "A and D have different data-types");

      // ensure that the compatible vector types are the same
      static_assert(std::is_same<typename MatrixA_::VectorTypeL, typename MatrixB_::VectorTypeL>::value,
        "A and B have different compatible L-vectors");
      static_assert(std::is_same<typename MatrixA_::VectorTypeR, typename MatrixD_::VectorTypeR>::value,
        "A and D have different compatible R-vectors");

      /// memory type
      typedef typename MatrixTypeA::MemType MemType;
      /// data type
      typedef typename MatrixTypeA::DataType DataType;

      /// Compatible L-vector type
      typedef TupleVector<typename MatrixTypeA::VectorTypeL, typename MatrixTypeD::VectorTypeL> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename MatrixTypeA::VectorTypeR, typename MatrixTypeB::VectorTypeR> VectorTypeR;

      /// dummy enum
      enum
      {
        /// number of row blocks (vertical size)
        num_row_blocks = 2,
        /// number of column blocks (horizontal size)
        num_col_blocks = 2
      };

    protected:
      /// sub-matrix A
      MatrixTypeA _matrix_a;
      /// sub-matrix B
      MatrixTypeB _matrix_b;
      /// sub-matrix C
      MatrixTypeD _matrix_d;

    public:
      /// default CTOR
      SaddlePointMatrix()
      {
      }

      /**
       * \brief Sub-Matrix move constructor
       *
       * \param[in] matrix_a, matrix_b, matrix_d
       * The three sub-matrices which are to be inserted.
       */
      explicit SaddlePointMatrix(
        MatrixTypeA&& matrix_a,
        MatrixTypeB&& matrix_b,
        MatrixTypeD&& matrix_d) :
        _matrix_a(matrix_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d)
      {
        ASSERT(_matrix_a.rows() == _matrix_b.rows(), "row count mismatch: A.rows != B.rows");
        ASSERT(_matrix_a.colums() == _matrix_d.colums(), "row count mismatch: A.cols != D.cols");
      }

      /// move ctor
      SaddlePointMatrix(SaddlePointMatrix&& other) :
        _matrix_a(std::move(other._matrix_a)),
        _matrix_b(std::move(other._matrix_b)),
        _matrix_d(std::move(other._matrix_d))
      {
      }

      /// move-assign operator
      SaddlePointMatrix& operator=(SaddlePointMatrix&& other)
      {
        if(this == &other)
          return *this;
        _matrix_a = std::move(other._matrix_a);
        _matrix_b = std::move(other._matrix_b);
        _matrix_d = std::move(other._matrix_d);
        return *this;
      }

      /// virtual destructor
      virtual ~SaddlePointMatrix()
      {
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      SaddlePointMatrix clone() const
      {
        return SaddlePointMatrix(_matrix_a.clone(), _matrix_b.clone(), _matrix_d.clone());
      }

      /// Returns the sub-matrix block A.
      MatrixTypeA& block_a()
      {
        return _matrix_a;
      }

      /// Returns the sub-matrix block A.
      const MatrixTypeA& block_a() const
      {
        return _matrix_a;
      }

      /// Returns the sub-matrix block B.
      MatrixTypeB& block_b()
      {
        return _matrix_b;
      }

      /// Returns the sub-matrix block B.
      const MatrixTypeB& block_b() const
      {
        return _matrix_b;
      }

      /// Returns the sub-matrix block D.
      MatrixTypeD& block_d()
      {
        return _matrix_d;
      }

      /// Returns the sub-matrix block D.
      const MatrixTypeD& block_d() const
      {
        return _matrix_d;
      }

      /**
       * \brief Returns a sub-matrix block.
       *
       * \tparam i_
       * The row index of the sub-matrix block that is to be returned.
       *
       * \tparam j_
       * The column index of the sub-matrix block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-matrix at position <em>(i_,j_)</em>.
       */
      template<Index i_, Index j_>
      typename SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::Type& at()
      {
        static_assert((i_ < Index(1)) || (j_ < Index(1)), "sub-matrix block (1,1) does not exist in SaddlePointMatrix");
        return SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::get(*this);
      }

      /** \copydoc at() */
      template<Index i_, Index j_>
      const typename SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::Type& at() const
      {
        static_assert((i_ < Index(1)) || (j_ < Index(1)), "sub-matrix block (1,1) does not exist in SaddlePointMatrix");
        return SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::get(*this);
      }

      /// Returns the total number of rows in this matrix.
      Index rows() const
      {
        return _matrix_a.rows() + _matrix_d.rows();
      }

      /// Returns the total number of columns in this matrix.
      Index columns() const
      {
        return _matrix_a.columns() + _matrix_b.columns();
      }

      /// Returns the total number of non-zeros in this matrix.
      Index used_elements() const
      {
        return _matrix_a.used_elements() + _matrix_b.used_elements() + _matrix_d.used_elements();
      }

      /**
       * \brief Applies this matrix onto a vector.
       *
       * This function performs
       *  \f[r \leftarrow this\cdot x \f]
       *
       * \param[out] r
       * The vector the receives the result.
       *
       * \param[in] x
       * The multiplicant vector.
       */
      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x)
      {
        block_a().template apply<Algo_>(r.template at<0>(), x.template at<0>());
        block_b().template apply<Algo_>(r.template at<0>(), x.template at<1>(), r.template at<0>(), DataType(1));
        block_d().template apply<Algo_>(r.template at<1>(), x.template at<0>());
      }

      /**
       * \brief Applies this matrix onto a vector.
       *
       * This function performs
       *  \f[r \leftarrow y + \alpha\cdot this\cdot x \f]
       *
       * \param[out] r
       * The vector the receives the result.
       *
       * \param[in] x
       * The multiplicant vector.
       *
       * \param[in] y
       * The summand vector
       * \param[in] alpha A scalar to scale the product with.
       */
      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1))
      {
        block_a().template apply<Algo_>(r.template at<0>(), x.template at<0>(), y.template at<0>(), alpha);
        block_b().template apply<Algo_>(r.template at<0>(), x.template at<1>(), r.template at<0>(), alpha);
        block_d().template apply<Algo_>(r.template at<1>(), x.template at<0>(), y.template at<1>(), alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(block_a().create_vector_l(), block_d().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(block_a().create_vector_r(), block_b().create_vector_r());
      }
    }; // class SaddlePointMatrix<...>

    /// \cond internal
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
    struct SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, Index(0), Index(0)>
    {
      typedef MatrixA_ Type;
      static Type& get(SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix)
      {
        return matrix.block_a();
      }
      static const Type& get(const SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix)
      {
        return matrix.block_a();
      }
    };

    template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
    struct SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, Index(0), Index(1)>
    {
      typedef MatrixB_ Type;
      static Type& get(SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix)
      {
        return matrix.block_b();
      }
      static const Type& get(const SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix)
      {
        return matrix.block_b();
      }
    };

    template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
    struct SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, Index(1), Index(0)>
    {
      typedef MatrixD_ Type;
      static Type& get(SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix)
      {
        return matrix.block_d();
      }
      static const Type& get(const SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix)
      {
        return matrix.block_d();
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP

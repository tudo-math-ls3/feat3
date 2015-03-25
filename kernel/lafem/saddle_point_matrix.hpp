#pragma once
#ifndef KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP
#define KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP 1

// includes, FEAST
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/dense_vector.hpp>

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
      /// index type
      typedef typename MatrixTypeA::IndexType IndexType;

      /// Compatible L-vector type
      typedef TupleVector<typename MatrixTypeA::VectorTypeL, typename MatrixTypeD::VectorTypeL> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename MatrixTypeA::VectorTypeR, typename MatrixTypeB::VectorTypeR> VectorTypeR;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class SaddlePointMatrix<typename MatrixA_::template ContainerType<Mem2_, DT2_, IT2_>,
                                                    typename MatrixB_::template ContainerType<Mem2_, DT2_, IT2_>,
                                                    typename MatrixD_::template ContainerType<Mem2_, DT2_, IT2_> >;

      /// number of row blocks (vertical size)
      static constexpr int num_row_blocks = 2;
      /// number of column blocks (horizontal size)
      static constexpr int num_col_blocks = 2;

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
        _matrix_a(std::move(matrix_a)),
        _matrix_b(std::move(matrix_b)),
        _matrix_d(std::move(matrix_d))
      {
        ASSERT(_matrix_a.rows() == _matrix_b.rows(), "row count mismatch: A.rows != B.rows");
        ASSERT(_matrix_a.columns() == _matrix_d.columns(), "row count mismatch: A.cols != D.cols");
      }

      /// move ctor
      SaddlePointMatrix(SaddlePointMatrix&& other) :
        _matrix_a(std::move(other._matrix_a)),
        _matrix_b(std::move(other._matrix_b)),
        _matrix_d(std::move(other._matrix_d))
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a saddle-point-matrix based on the source file.
       */
      explicit SaddlePointMatrix(FileMode mode, String filename)
      {
        CONTEXT("When creating SaddlePointMatrix");

        read_from(mode, filename);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
      {
        CONTEXT("When reading in SaddlePointMatrix");

        String directory;
        auto found = filename.rfind("/");
        if (found != std::string::npos)
        {
          directory = filename.substr(0, found + 1);
        }

        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);

        String line;
        std::getline(file, line);
        if (line.find("%%MatrixMarket saddlepointmatrix coordinate real general") == String::npos)
          throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a complatible file");

        do {
          if (file.eof())
            throw InternalError(__func__, __FILE__, __LINE__, "Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        MatrixA_ tmp_a(mode, directory + line);
        _matrix_a = std::move(tmp_a);

        std::getline(file, line);
        line.trim_me();
        MatrixB_ tmp_b(mode, directory + line);
        _matrix_b = std::move(tmp_b);

        std::getline(file, line);
        line.trim_me();
        MatrixD_ tmp_d(mode, directory + line);
        _matrix_d = std::move(tmp_d);

        file.close();
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
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        CONTEXT("When writing out SaddlePointMatrix");

        std::ofstream file(filename.c_str(), std::ofstream::out);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);

        String suffix, directory;
        auto found = filename.rfind(".");
        if (found != std::string::npos)
        {
          suffix = filename.substr(found);
          filename.erase(found);
        }
        found = filename.rfind("/");
        if (found != std::string::npos)
        {
          directory = filename.substr(0, found + 1);
          filename.erase(0, found + 1);
        }

        file << "%%MatrixMarket saddlepointmatrix coordinate real general" << std::endl;
        file << filename << "_a" << suffix << std::endl;
        file << filename << "_b" << suffix << std::endl;
        file << filename << "_d" << suffix << std::endl;

        file.close();

        _matrix_a.write_out(mode, directory + filename + "_a" + suffix);
        _matrix_b.write_out(mode, directory + filename + "_b" + suffix);
        _matrix_d.write_out(mode, directory + filename + "_d" + suffix);
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

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("SaddePointMatrix<") + MatrixTypeA::name() + "," +
          MatrixTypeB::name() + "," +
          MatrixTypeD::name() + ">";
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
#ifdef FEAST_COMPILER_MICROSOFT
      template<typename Algo_, typename VectorL_, typename VectorR_>
      void apply(VectorL_& r, const VectorR_& x) const
#else
      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x) const
#endif
      {
        block_a().template apply<Algo_>(r.template at<0>(), x.template at<0>());
        block_b().template apply<Algo_>(r.template at<0>(), x.template at<1>(), r.template at<0>(), DataType(1));
        block_d().template apply<Algo_>(r.template at<1>(), x.template at<0>());
      }

      template<typename Algo_>
      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, block_a().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, block_d().rows(), block_a().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, block_a().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, block_b().columns(), block_a().columns());

        block_a().template apply<Algo_>(r_first, x_first);
        block_b().template apply<Algo_>(r_first, x_rest, r_first, DataType(1));
        block_d().template apply<Algo_>(r_rest, x_first);
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
#ifdef FEAST_COMPILER_MICROSOFT
      template<typename Algo_, typename VectorL_, typename VectorR_>
      void apply(VectorL_& r, const VectorR_& x, const VectorL_& y, DataType alpha = DataType(1)) const
#else
      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
#endif
      {
        block_a().template apply<Algo_>(r.template at<0>(), x.template at<0>(), y.template at<0>(), alpha);
        block_b().template apply<Algo_>(r.template at<0>(), x.template at<1>(), r.template at<0>(), alpha);
        block_d().template apply<Algo_>(r.template at<1>(), x.template at<0>(), y.template at<1>(), alpha);
      }

      template<typename Algo_>
      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x,
                 const DenseVector<MemType, DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, block_a().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, block_d().rows(), block_a().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, block_a().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, block_b().columns(), block_a().columns());

        DenseVector<MemType, DataType, IndexType> y_first(y, block_a().rows(), 0);
        DenseVector<MemType, DataType, IndexType> y_rest(y, block_d().rows(), block_a().rows());

        block_a().template apply<Algo_>(r_first, x_first, y_first, alpha);
        block_b().template apply<Algo_>(r_first, x_rest, r_first, alpha);
        block_d().template apply<Algo_>(r_rest, x_first, y_rest, alpha);
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

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        const Index arows(this->block_a().rows());

        if (row < arows)
        {
          return this->block_a().get_length_of_line(row) + this->block_b().get_length_of_line(row);
        }
        else
        {
          return this->block_d().get_length_of_line(row - arows);
        }
      }

      /// \cond internal
      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, typename MatrixA_::DataType * const pval_set, typename MatrixA_::IndexType * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const Index arows(this->block_a().rows());

        if (row < arows)
        {
          const Index length_of_a(this->block_a().get_length_of_line(row));

          this->block_a().set_line(row, pval_set, pcol_set, col_start, stride);
          this->block_b().set_line(row, pval_set + stride * length_of_a, pcol_set + stride * length_of_a, col_start + this->block_a().columns(), stride);
        }
        else
        {
          this->block_d().set_line(row - arows, pval_set, pcol_set, col_start, stride);
        }
      }
      /// \endcond

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const typename SaddlePointMatrix::template ContainerType<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting SaddlePointMatrix");

        this->block_a().convert(other.block_a());
        this->block_b().convert(other.block_b());
        this->block_d().convert(other.block_d());
      }

      /**
       * \brief SaddlePointMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_>
      friend bool operator== (const SaddlePointMatrix & a, const ContainerType<Mem2_> & b)
      {
        CONTEXT("When comparing SaddlePointMatrices");

        return (a.name() == b.name()) && (a.block_a() == b.block_a())
          && (a.block_b() == b.block_b()) && (a.block_d() == b.block_d());
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

#pragma once
#ifndef KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP
#define KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP 1

// includes, FEAT
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/container.hpp>

// includes, system
#include <type_traits>
#include <fstream>

namespace FEAT
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
      int i_,
      int j_>
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

      /// this typedef lets you create a matrix container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

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
        XASSERT(_matrix_a.rows() == _matrix_b.rows());
        XASSERT(_matrix_a.columns() == _matrix_d.columns());
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
      SaddlePointMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return SaddlePointMatrix(_matrix_a.clone(mode), _matrix_b.clone(mode), _matrix_d.clone(mode));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _matrix_a.bytes() + _matrix_b.bytes() + _matrix_d.bytes();
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
      template<int i_, int j_>
      typename SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::Type& at()
      {
        static_assert((0 <= i_) && (0 <= j_), "invalid sub-matrix index");
        static_assert((i_ <= 1) && (j_ <= 1), "invalid sub-matrix index");
        static_assert((i_ < 1) || (j_ < 1), "sub-matrix block (1,1) does not exist in SaddlePointMatrix");
        return SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::get(*this);
      }

      /** \copydoc at() */
      template<int i_, int j_>
      const typename SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::Type& at() const
      {
        static_assert((0 <= i_) && (0 <= j_), "invalid sub-matrix index");
        static_assert((i_ <= 1) && (j_ <= 1), "invalid sub-matrix index");
        static_assert((i_ < 1) || (j_ < 1), "sub-matrix block (1,1) does not exist in SaddlePointMatrix");
        return SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, i_, j_>::get(*this);
      }

      /**
       * \brief Returns the total number of rows in this matrix.
       *
       * \returns Matrix row count if perspective_ = false.
       * \returns Raw matrix row count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        return _matrix_a.template rows<perspective_>() + _matrix_d.template rows<perspective_>();
      }

      /**
       * \brief Returns the total number of columns in this matrix.
       *
       * \returns Matrix column count if perspective_ = false.
       * \returns Raw matrix column count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        return _matrix_a.template columns<perspective_>() + _matrix_b.template columns<perspective_>();
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \returns Matrix non zero element count if perspective_ = false.
       * \returns Raw matrix non zero element count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        return _matrix_a.template used_elements<perspective_>() + _matrix_b.template used_elements<perspective_>() + _matrix_d.template used_elements<perspective_>();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("SaddePointMatrix<") + MatrixTypeA::name() + "," +
          MatrixTypeB::name() + "," +
          MatrixTypeD::name() + ">";
      }

      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return rows<perspective_>() * columns<perspective_>();
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] value
       * The value to which the matrix' entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        block_a().format(value);
        block_b().format(value);
        block_d().format(value);
      }

      /**
       * \brief Free all allocated arrays
       */
      void clear()
      {
        block_a().clear();
        block_b().clear();
        block_d().clear();
      }

      /// extract main diagonal vector from matrix
      void extract_diag(VectorTypeL& diag) const
      {
        _matrix_a.extract_diag(diag.first());
        diag.rest().format();
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
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        block_a().apply(r.template at<0>(), x.template at<0>());
        block_b().apply(r.template at<0>(), x.template at<1>(), r.template at<0>(), DataType(1));
        block_d().apply(r.template at<1>(), x.template at<0>());
      }

      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, block_a().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, block_d().rows(), block_a().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, block_a().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, block_b().columns(), block_a().columns());

        block_a().apply(r_first, x_first);
        block_b().apply(r_first, x_rest, r_first, DataType(1));
        block_d().apply(r_rest, x_first);
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
       *
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        block_a().apply(r.template at<0>(), x.template at<0>(), y.template at<0>(), alpha);
        block_b().apply(r.template at<0>(), x.template at<1>(), r.template at<0>(), alpha);
        block_d().apply(r.template at<1>(), x.template at<0>(), y.template at<1>(), alpha);
      }

      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x,
                 const DenseVector<MemType, DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, block_a().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, block_d().rows(), block_a().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, block_a().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, block_b().columns(), block_a().columns());

        DenseVector<MemType, DataType, IndexType> y_first(y, block_a().rows(), 0);
        DenseVector<MemType, DataType, IndexType> y_rest(y, block_d().rows(), block_a().rows());

        block_a().apply(r_first, x_first, y_first, alpha);
        block_b().apply(r_first, x_rest, r_first, alpha);
        block_d().apply(r_rest, x_first, y_rest, alpha);
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
        const Index arows(this->block_a().template rows<Perspective::pod>());

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
        const Index arows(this->block_a().template rows<Perspective::pod>());

        if (row < arows)
        {
          const Index length_of_a(this->block_a().get_length_of_line(row));

          this->block_a().set_line(row, pval_set, pcol_set, col_start, stride);
          this->block_b().set_line(row, pval_set + stride * length_of_a, pcol_set + stride * length_of_a, col_start + this->block_a().template columns<Perspective::pod>(), stride);
        }
        else
        {
          this->block_d().set_line(row - arows, pval_set, pcol_set, col_start, stride);
        }
      }

      void set_line_reverse(const Index row, typename MatrixA_::DataType * const pval_set, const Index stride = 1)
      {
        const Index arows(this->block_a().template rows<Perspective::pod>());

        if (row < arows)
        {
          const Index length_of_a(this->block_a().get_length_of_line(row));

          this->block_a().set_line_reverse(row, pval_set, stride);
          this->block_b().set_line_reverse(row, pval_set + stride * length_of_a, stride);
        }
        else
        {
          this->block_d().set_line_reverse(row - arows, pval_set, stride);
        }
      }
      /// \endcond

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename MatrixA2_, typename MatrixB2_, typename MatrixD2_>
      void convert(const SaddlePointMatrix<MatrixA2_, MatrixB2_, MatrixD2_> & other)
      {
        this->block_a().convert(other.block_a());
        this->block_b().convert(other.block_b());
        this->block_d().convert(other.block_d());
      }

      template <typename MatrixA2_, typename MatrixB2_, typename MatrixD2_>
      void convert_reverse(SaddlePointMatrix<MatrixA2_, MatrixB2_, MatrixD2_> & other) const
      {
        this->block_a().convert_reverse(other.block_a());
        this->block_b().convert_reverse(other.block_b());
        this->block_d().convert_reverse(other.block_d());
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
        return (a.name() == b.name()) && (a.block_a() == b.block_a())
          && (a.block_b() == b.block_b()) && (a.block_d() == b.block_d());
      }
    }; // class SaddlePointMatrix<...>

    /// \cond internal
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
    struct SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, 0, 0>
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
    struct SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, 0, 1>
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
    struct SaddlePointMatrixElement<MatrixA_, MatrixB_, MatrixD_, 1, 0>
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
} // namespace FEAT

#endif // KERNEL_LAFEM_SADDLE_POINT_MATRIX_HPP

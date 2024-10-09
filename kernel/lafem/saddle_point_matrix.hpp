// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
     * This class template is a helper to realize the "at" member function of the SaddlePointMatrix class template.
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

      // ensure that all matrices have the same and data-types
      static_assert(std::is_same<typename MatrixA_::DataType, typename MatrixB_::DataType>::value,
                    "A and B have different data-types");
      static_assert(std::is_same<typename MatrixA_::DataType, typename MatrixD_::DataType>::value,
                    "A and D have different data-types");

      // ensure that the compatible vector types are the same
      static_assert(std::is_same<typename MatrixA_::VectorTypeL, typename MatrixB_::VectorTypeL>::value,
                    "A and B have different compatible L-vectors");
      static_assert(std::is_same<typename MatrixA_::VectorTypeR, typename MatrixD_::VectorTypeR>::value,
                    "A and D have different compatible R-vectors");

      /// data type
      typedef typename MatrixTypeA::DataType DataType;
      /// index type
      typedef typename MatrixTypeA::IndexType IndexType;

      /// Compatible L-vector type
      typedef TupleVector<typename MatrixTypeA::VectorTypeL, typename MatrixTypeD::VectorTypeL> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename MatrixTypeA::VectorTypeR, typename MatrixTypeB::VectorTypeR> VectorTypeR;
      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = SaddlePointMatrix<typename MatrixA_::template ContainerType<DT2_, IT2_>,
                                              typename MatrixB_::template ContainerType<DT2_, IT2_>,
                                              typename MatrixD_::template ContainerType<DT2_, IT2_> >;

      /// this typedef lets you create a matrix container with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

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
      explicit SaddlePointMatrix(FileMode mode, const String& filename)
      {
        read_from(mode, filename);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, const String& filename)
      {
        String directory;
        auto found = filename.rfind("/");
        if (found != std::string::npos)
        {
          directory = filename.substr(0, found + 1);
        }

        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);

        String line;
        std::getline(file, line);
        if (line.find("%%MatrixMarket saddlepointmatrix coordinate real general") == String::npos)
          XABORTM("Input-file is not a complatible file");

        do {
          if (file.eof())
            XABORTM("Wrong Input-file");
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
          XABORTM("Unable to open Matrix file " + filename);

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

        file << "%%MatrixMarket saddlepointmatrix coordinate real general" << "\n";
        file << filename << "_a" << suffix << "\n";
        file << filename << "_b" << suffix << "\n";
        file << filename << "_d" << suffix << "\n";

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

      void lump_rows(VectorTypeL& lump, bool /*sync*/ = true) const
      {
        _matrix_a.lump_rows(lump.first());
        lump.rest().format();
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

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        DenseVector<DataType, IndexType> r_first(r, block_a().rows(), 0);
        DenseVector<DataType, IndexType> r_rest(r, block_d().rows(), block_a().rows());

        DenseVector<DataType, IndexType> x_first(x, block_a().columns(), 0);
        DenseVector<DataType, IndexType> x_rest(x, block_b().columns(), block_a().columns());

        block_a().apply(r_first, x_first);
        block_b().apply(r_first, x_rest, r_first, DataType(1));
        block_d().apply(r_rest, x_first);
      }

      /**
      * \brief Applies this matrix onto a vector.
      *
      * This function performs
      *  \f[r \leftarrow this^\top \cdot x \f]
      *
      * \param[out] r
      * The vector the receives the result.
      *
      * \param[in] x
      * The multiplicant vector.
      */
      void apply_transposed(VectorTypeR& r, const VectorTypeL& x) const
      {
        block_a().apply_transposed(r.template at<0>(), x.template at<0>());
        block_d().apply_transposed(r.template at<0>(), x.template at<1>(), r.template at<0>(), DataType(1));
        block_b().apply_transposed(r.template at<1>(), x.template at<0>());
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");

        DenseVector<DataType, IndexType> r_first(r, block_a().columns(), 0);
        DenseVector<DataType, IndexType> r_rest(r, block_b().columns(), block_a().columns());

        DenseVector<DataType, IndexType> x_first(x, block_a().rows(), 0);
        DenseVector<DataType, IndexType> x_rest(x, block_d().rows(), block_a().rows());

        block_a().apply_transposed(r_first, x_first);
        block_d().apply_transposed(r_first, x_rest, r_first, DataType(1));
        block_b().applytransposed(r_rest, x_first);
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

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
                 const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

        DenseVector<DataType, IndexType> r_first(r, block_a().rows(), 0);
        DenseVector<DataType, IndexType> r_rest(r, block_d().rows(), block_a().rows());

        DenseVector<DataType, IndexType> x_first(x, block_a().columns(), 0);
        DenseVector<DataType, IndexType> x_rest(x, block_b().columns(), block_a().columns());

        DenseVector<DataType, IndexType> y_first(y, block_a().rows(), 0);
        DenseVector<DataType, IndexType> y_rest(y, block_d().rows(), block_a().rows());

        block_a().apply(r_first, x_first, y_first, alpha);
        block_b().apply(r_first, x_rest, r_first, alpha);
        block_d().apply(r_rest, x_first, y_rest, alpha);
      }

      /**
      * \brief Applies this matrix onto a vector.
      *
      * This function performs
      *  \f[r \leftarrow y + \alpha\cdot this^\top \cdot x \f]
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
      void apply_transposed(VectorTypeR& r, const VectorTypeL& x, const VectorTypeR& y, DataType alpha = DataType(1)) const
      {
        block_a().apply_transposed(r.template at<0>(), x.template at<0>(), y.template at<0>(), alpha);
        block_d().apply_transposed(r.template at<0>(), x.template at<1>(), r.template at<0>(), alpha);
        block_b().apply_transposed(r.template at<1>(), x.template at<0>(), y.template at<1>(), alpha);
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
        const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->columns(), "Vector size of y does not match!");

        DenseVector<DataType, IndexType> r_first(r, block_a().columns(), 0);
        DenseVector<DataType, IndexType> r_rest(r, block_b().columns(), block_a().columns());

        DenseVector<DataType, IndexType> x_first(x, block_a().rows(), 0);
        DenseVector<DataType, IndexType> x_rest(x, block_d().rows(), block_a().rows());

        DenseVector<DataType, IndexType> y_first(y, block_a().columns(), 0);
        DenseVector<DataType, IndexType> y_rest(y, block_b().columns(), block_a().columns());

        block_a().apply_transposed(r_first, x_first, y_first, alpha);
        block_d().apply_transposed(r_first, x_rest, r_first, alpha);
        block_b().apply_transposed(r_rest, x_first, y_rest, alpha);
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

      void set_line_reverse(const Index row, const typename MatrixA_::DataType * const pval_set, const Index stride = 1)
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

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return (3 * sizeof(std::uint64_t)) + this->block_a().get_checkpoint_size(config) + this->block_b().get_checkpoint_size(config) + this->block_d().get_checkpoint_size(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        std::uint64_t isize = *(std::uint64_t*) data.data(); //get size of checkpointed block a
        std::vector<char>::iterator start = std::begin(data) + sizeof(std::uint64_t); //get iterator at the beginning of block a
        std::vector<char>::iterator end = std::begin(data) + sizeof(std::uint64_t) + (int) isize; //get iterator at the beginning of block a
        std::vector<char> buffer_a(start, end); //copy the data of block a to fresh vector
        this->block_a().restore_from_checkpoint_data(buffer_a);

        data.erase(std::begin(data), end);
        isize = *(std::uint64_t*) data.data();
        start = std::begin(data) + sizeof(std::uint64_t);
        end = std::begin(data) + sizeof(std::uint64_t) + (int) isize;
        std::vector<char> buffer_b(start, end);
        this->block_b().restore_from_checkpoint_data(buffer_b);

        data.erase(std::begin(data), end);
        this->block_d().restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        std::size_t old_size = data.size();
        data.insert(std::end(data), sizeof(std::uint64_t),0); //add placeholder
        std::uint64_t ireal_size = this->block_a().set_checkpoint_data(data, config); //add data of block a to the overall checkpoint
        std::uint64_t ret_size = ireal_size;
        char* csize = reinterpret_cast<char*>(&ireal_size);
        for(std::size_t i(0) ; i < sizeof(std::uint64_t) ; ++i)  //overwrite the guessed datalength
        {
          data[old_size + i] = csize[i];
        }

        old_size = data.size();
        data.insert(std::end(data), sizeof(std::uint64_t), 0); //add placeholder
        ireal_size = this->block_b().set_checkpoint_data(data, config); //add data of block b to the overall checkpoint
        ret_size += ireal_size;
        csize = reinterpret_cast<char*>(&ireal_size);
        for(std::size_t i(0) ; i < sizeof(std::uint64_t) ; ++i)  //overwrite the guessed datalength
        {
          data[old_size + i] = csize[i];
        }

        return 2*sizeof(std::uint64_t) + ret_size + this->block_d().set_checkpoint_data(data, config); //add data of block d to the overall checkpoint
      }

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
      friend bool operator== (const SaddlePointMatrix & a, const SaddlePointMatrix & b)
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

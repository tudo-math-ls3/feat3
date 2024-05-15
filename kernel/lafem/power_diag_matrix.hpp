// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP
#define KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP 1

// includes, FEAT
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/meta_element.hpp>
#include <kernel/lafem/container.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Power-Diag-Matrix meta class template
     *
     * This class template implements a diagonal composition of \e n sub-matrices of the same class.
     * This can be interpreted as a diagonal m-by-m matrix of other matrices.
     *
     * \tparam SubType_
     * The type of the sub-matrix.
     *
     * \tparam blocks_
     * The number of sub-matrix blocks.
     *
     * \author Peter Zajac
     */
    template<
      typename SubType_,
      int blocks_>
    class PowerDiagMatrix
    {
      // Note: the case = 1 is specialized below
      static_assert(blocks_ > 1, "invalid block size");

      // declare this class template as a friend for recursive inheritance
      template<typename, int>
      friend class PowerDiagMatrix;

      /// rest-class typedef
      typedef PowerDiagMatrix<SubType_, blocks_-1> RestClass;

    public:
      /// sub-matrix type
      typedef SubType_ SubMatrixType;
      /// sub-matrix data type
      typedef typename SubMatrixType::DataType DataType;
      /// sub-matrix index type
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutId layout_id = SubMatrixType::layout_id;
      /// Compatible L-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeL, blocks_> VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, blocks_> VectorTypeR;
      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = PowerDiagMatrix<typename SubType_::template ContainerType<DT2_, IT2_>, blocks_>;

      /// this typedef lets you create a matrix container with new Datatape and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

      /// number of row blocks (vertical size)
      static constexpr int num_row_blocks = blocks_;
      /// number of column blocks (horizontal size)
      static constexpr int num_col_blocks = blocks_;

    protected:
      /// the first sub-matrix
      SubMatrixType _first;
      /// the remaining part
      RestClass _rest;

      /// base-class constructor; this one is protected for a reason
      explicit PowerDiagMatrix(SubMatrixType&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

    public:
      /// default ctor
      PowerDiagMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerDiagMatrix(const SparseLayout<IndexType, layout_id>& layout) :
        _first(layout),
        _rest(layout)
      {
      }

      /// move ctor
      PowerDiagMatrix(PowerDiagMatrix&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a power-diag-point-matrix based on the source file.
       */
      explicit PowerDiagMatrix(FileMode mode, String filename)
      {
        read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * \note This constructor is used internally when reading a file
       *
       * Creates a power-diag-matrix based on the source filestream.
       */
      explicit PowerDiagMatrix(FileMode mode, std::istream& file, String directory = "")
      {
        String line;
        do {
          if (file.eof())
            XABORTM("Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        SubMatrixType tmp_first(mode, directory + line);
        _first = std::move(tmp_first);

        RestClass tmp_rest(mode, file, directory);
        _rest = std::move(tmp_rest);
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
          XABORTM("Unable to open Matrix file " + filename);

        String line;
        std::getline(file, line);
        if (line.find("%%MatrixMarket powerdiagmatrix coordinate real general") == String::npos)
          XABORTM("Input-file is not a complatible file");

        PowerDiagMatrix other(mode, file, directory);

        _first = std::move(other._first);
        _rest = std::move(other._rest);

        file.close();
      }

      /// move-assign operator
      PowerDiagMatrix& operator=(PowerDiagMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerDiagMatrix(const PowerDiagMatrix&) = delete;
      /// deleted copy-assign operator
      PowerDiagMatrix& operator=(const PowerDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerDiagMatrix()
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

        file << "%%MatrixMarket powerdiagmatrix coordinate real general" << "\n";
        for (Index i(1); i <= blocks_; ++i)
        {
          file << filename << "_pd" << i << suffix << "\n";
        }

        file.close();

        this->write_out_submatrices(mode, directory, filename, suffix);
      }

      /**
       * \brief Write out submatrices to file.
       *
       * \param[in] mode The used file format.
       * \param[in] directory The directory of the matrix-files.
       * \param[in] prefix The prefix of the matrix-files.
       * \param[in] suffix The suffix of the matrix-files.
       */
      void write_out_submatrices(FileMode mode, String directory, String prefix, String suffix, Index length = blocks_) const
      {
        _first.write_out(mode, directory + prefix + "_pd" + stringify(length + 1 - blocks_) + suffix);
        _rest.write_out_submatrices(mode, directory, prefix, suffix, length);
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      PowerDiagMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return PowerDiagMatrix(_first.clone(mode), _rest.clone(mode));
      }

      /** \brief Clone operation
       *
       * Become a copy of a given container.
       *
       * \param[in] other The source container.
       * \param[in] clone_mode The actual cloning procedure
       *
       */
      void clone(const PowerDiagMatrix & other, CloneMode clone_mode = CloneMode::Weak)
      {
        (*this)=other.clone(clone_mode);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
      }

      /**
       * \brief Creates and returns the transpose of this matrix.
       */
      PowerDiagMatrix transpose() const
      {
        return PowerDiagMatrix(_first.transpose(), _rest.transpose());
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
      SubMatrixType& at()
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        static_assert((0 <= i_) && (i_ < blocks_), "invalid sub-matrix index");
        return PowerElement<i_, SubMatrixType>::get(*this);
      }

      /** \copydoc at() */
      template<int i_, int j_>
      const SubMatrixType& at() const
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        static_assert((0 <= i_) && (i_ < blocks_), "invalid sub-matrix index");
        return PowerElement<i_, SubMatrixType>::get(*this);
      }

      /**
       * \brief Returns a sub-matrix block.
       *
       * \param[in] i, j
       * The indices of the sub-matrix block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-matrix at position (i,j).
       */
      SubMatrixType& get(int i, int j)
      {
        XASSERTM((i == j) && (0 <= i) && (i < blocks_), "invalid sub-matrix index");
        return (i == 0) ? _first : _rest.get(i-1, j-1);
      }

      /** \copydoc get() */
      const SubMatrixType& get(int i, int j) const
      {
        XASSERTM((i == j) && (0 <= i) && (i < blocks_), "invalid sub-matrix index");
        return (i == 0) ? _first : _rest.get(i-1, j-1);
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return sizeof(std::uint64_t) + _first.get_checkpoint_size(config) + _rest.get_checkpoint_size(config); //sizeof(std::uint64_t) bits needed to store lenght of checkpointed _first
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        std::uint64_t isize = *(std::uint64_t*) data.data(); //get size of checkpointed _first
        std::vector<char>::iterator start = std::begin(data) + sizeof(std::uint64_t);
        std::vector<char>::iterator last_of_first = std::begin(data) + sizeof(std::uint64_t) + (int) isize;
        std::vector<char> buffer_first(start, last_of_first);
        _first.restore_from_checkpoint_data(buffer_first);

        data.erase(std::begin(data), last_of_first);
        _rest.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        std::size_t old_size = data.size();
        data.insert(std::end(data), sizeof(std::uint64_t), 0); //add placeholder
        std::uint64_t ireal_size = _first.set_checkpoint_data(data, config); //add data of _first to the overall checkpoint and save its size
        char* csize = reinterpret_cast<char*>(&ireal_size);
        for(std::size_t i(0) ; i < sizeof(std::uint64_t) ; ++i)  //overwrite the guessed datalength
        {
          data[old_size + i] = csize[i];
        }

        return sizeof(std::uint64_t) + ireal_size + _rest.set_checkpoint_data(data, config); //generate and add checkpoint data for the _rest
      }

      /// \cond internal
      SubMatrixType& first()
      {
        return _first;
      }

      const SubMatrixType& first() const
      {
        return _first;
      }

      RestClass& rest()
      {
        return _rest;
      }

      const RestClass& rest() const
      {
        return _rest;
      }

      int row_blocks() const
      {
        return num_row_blocks;
      }

      int col_blocks() const
      {
        return num_col_blocks;
      }
      /// \endcond

      /**
       * \brief Returns the total number of rows in this matrix.
       *
       * \returns Matrix row count if perspective_ = false.
       * \returns Raw matrix row count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        return first().template rows<perspective_>() + rest().template rows<perspective_>();
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
        return first().template columns<perspective_>() + rest().template columns<perspective_>();
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
        return first().template used_elements<perspective_>() + rest().template used_elements<perspective_>();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("PowerDiagMatrix<") + SubMatrixType::name() + "," + stringify(blocks_) + ">";
      }

      Index size() const
      {
        return rows() * columns();
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] value
       * The value to which the matrix' entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
      }

      /**
       * \brief Free all allocated arrays
       */
      void clear()
      {
        first().clear();
        rest().clear();
      }

      /// extract main diagonal vector from matrix
      void extract_diag(VectorTypeL& diag) const
      {
        first().extract_diag(diag.first());
        rest().extract_diag(diag.rest());
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
        first().apply(r.first(), x.first());
        rest().apply(r.rest(), x.rest());
      }

      void apply(DenseVector<DataType , IndexType>& r, const DenseVector<DataType , IndexType>& x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        DenseVector<DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        first().apply(r_first, x_first);
        rest().apply(r_rest, x_rest);
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
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r.first(), x.first(), y.first(), alpha);
        rest().apply(r.rest(), x.rest(), y.rest(), alpha);
      }

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
                 const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

        DenseVector<DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        DenseVector<DataType, IndexType> y_first(y, first().rows(), 0);
        DenseVector<DataType, IndexType> y_rest(y, rest().rows(), first().rows());

        first().apply(r_first, x_first, y_first, alpha);
        rest().apply(r_rest, x_rest, y_rest, alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(first().create_vector_l(), rest().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(first().create_vector_r(), rest().create_vector_r());
      }

      void scale_rows(const PowerDiagMatrix& a, const VectorTypeL& w)
      {
        first().scale_rows(a.first(), w.first());
        rest().scale_rows(a.rest(), w.rest());
      }

      void scale_cols(const PowerDiagMatrix& a, const VectorTypeR& w)
      {
        first().scale_cols(a.first(), w.first());
        rest().scale_cols(a.rest(), w.rest());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        const Index brows(this->first().template rows<Perspective::pod>());

        if (row < brows)
        {
          return this->first().get_length_of_line(row);
        }
        else
        {
          return this->rest().get_length_of_line(row - brows);
        }
      }

      /// \cond internal
      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const Index brows(this->first().template rows<Perspective::pod>());
        const Index bcolumns(this->first().template columns<Perspective::pod>());

        if (row < brows)
        {
          this->first().set_line(row, pval_set, pcol_set, col_start, stride);
        }
        else
        {
          this->rest().set_line(row - brows, pval_set, pcol_set, col_start + bcolumns, stride);
        }
      }

      void set_line_reverse(const Index row, const DataType * const pval_set, const Index stride = 1)
      {
        const Index brows(this->first().template rows<Perspective::pod>());

        if (row < brows)
        {
          this->first().set_line_reverse(row, pval_set, stride);
        }
        else
        {
          this->rest().set_line_reverse(row - brows, pval_set, stride);
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
      template <typename SubType2_>
      void convert(const PowerDiagMatrix<SubType2_, blocks_> & other)
      {
        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }

      template <typename SubType2_>
      void convert_reverse(PowerDiagMatrix<SubType2_, blocks_> & other) const
      {
        this->first().convert_reverse(other.first());
        this->rest().convert_reverse(other.rest());
      }

      /**
       * \brief PowerDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      friend bool operator== (const PowerDiagMatrix & a, const PowerDiagMatrix & b)
      {
        return (a.name() == b.name()) && (a.first() == b.first()) && (a.rest() == b.rest());
      }
    };

    /// \cond internal
    template<typename SubType_>
    class PowerDiagMatrix<SubType_, 1>
    {
      template<typename, int>
      friend class PowerDiagMatrix;

    public:
      typedef SubType_ SubMatrixType;
      typedef typename SubMatrixType::DataType DataType;
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutId layout_id = SubMatrixType::layout_id;
      /// Compatible L-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeL, 1> VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, 1> VectorTypeR;
      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = PowerDiagMatrix<typename SubType_::template ContainerType<DT2_, IT2_>, 1>;

      /// this typedef lets you create a matrix container with new Datatape and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

      static constexpr int num_row_blocks = 1;
      static constexpr int num_col_blocks = 1;

    protected:
      SubMatrixType _first;

      /// base-class constructor; this one is protected for a reason
      explicit PowerDiagMatrix(SubMatrixType&& the_first) :
        _first(std::move(the_first))
      {
      }

    public:
      /// default ctor
      PowerDiagMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerDiagMatrix(const SparseLayout<IndexType, layout_id>& layout) :
        _first(layout)
      {
      }

      /// move ctor
      PowerDiagMatrix(PowerDiagMatrix&& other) :
        _first(std::move(other._first))
      {
      }

      /// move-assign operator
      PowerDiagMatrix& operator=(PowerDiagMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      /// file-input ctor
      explicit PowerDiagMatrix(FileMode mode, const String& filename)
      {
        read_from(mode, filename);
      }

      /// filestream-input ctor
      explicit PowerDiagMatrix(FileMode mode, std::istream& file, const String& directory = "")
      {
        String line;
        do {
          if (file.eof())
            XABORTM("Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        SubMatrixType tmp_first(mode, directory + line);
        _first = std::move(tmp_first);
      }

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

        PowerDiagMatrix other(mode, file, directory);

        _first = std::move(other._first);

        file.close();
      }

      /// deleted copy-ctor
      PowerDiagMatrix(const PowerDiagMatrix&) = delete;
      /// deleted copy-assign operator
      PowerDiagMatrix& operator=(const PowerDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerDiagMatrix()
      {
      }

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

        file << "%%MatrixMarket powerdiagmatrix coordinate real general" << "\n";
        file << filename << "_pd" << 1 << suffix << "\n";

        file.close();

        this->write_out_submatrices(mode, directory, filename, suffix);
      }

      void write_out_submatrices(FileMode mode, const String& directory, const String& prefix, const String& suffix, Index length = 1) const
      {
        _first.write_out(mode, directory + prefix + "_pd" + stringify(length) + suffix);
      }

      PowerDiagMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return PowerDiagMatrix(_first.clone(mode));
      }

      void clone(const PowerDiagMatrix & other, CloneMode clone_mode = CloneMode::Weak)
      {
        (*this)=other.clone(clone_mode);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes();
      }

      PowerDiagMatrix transpose() const
      {
        return PowerDiagMatrix(_first.transpose());
      }

      template<int i, int j>
      SubMatrixType& at()
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _first;
      }

      template<int i, int j>
      const SubMatrixType& at() const
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _first;
      }

      SubMatrixType& get(int i, int j)
      {
        XASSERTM((i == 0) && (j == 0), "invalid sub-matrix index");
        return _first;
      }

      const SubMatrixType& get(int i, int j) const
      {
        XASSERTM((i == 0) && (j == 0), "invalid sub-matrix index");
        return _first;
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return _first.get_checkpoint_size(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        _first.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        return _first.set_checkpoint_data(data, config);
      }

      SubMatrixType& first()
      {
        return _first;
      }

      const SubMatrixType& first() const
      {
        return _first;
      }

      int row_blocks() const
      {
        return 1;
      }

      int col_blocks() const
      {
        return 1;
      }

      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        return first().template rows<perspective_>();
      }

      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        return first().template columns<perspective_>();
      }

      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        return first().template used_elements<perspective_>();
      }

      static String name()
      {
        return String("PowerDiagMatrix<") + SubMatrixType::name() + "," + stringify(1) + ">";
      }

      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return rows<perspective_>() * columns<perspective_>();
      }

      void format(DataType value = DataType(0))
      {
        first().format(value);
      }

      void clear()
      {
        first().clear();
      }

      void extract_diag(VectorTypeL& diag) const
      {
        first().extract_diag(diag.first());
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().apply(r.first(), x.first());
      }

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        first().apply(r, x);
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r.first(), x.first(), y.first(), alpha);
      }

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
                 const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        first().apply(r, x, y, alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(first().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(first().create_vector_r());
      }

      void scale_rows(const PowerDiagMatrix& a, const VectorTypeL& w)
      {
        first().scale_rows(a.first(), w.first());
      }

      void scale_cols(const PowerDiagMatrix& a, const VectorTypeR& w)
      {
        first().scale_cols(a.first(), w.first());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        return this->first().get_length_of_line(row);
      }

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        this->first().set_line(row, pval_set, pcol_set, col_start, stride);
      }

      void set_line_reverse(const Index row, const DataType * const pval_set, const Index stride = 1)
      {
        this->first().set_line_reverse(row, pval_set, stride);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename SubType2_>
      void convert(const PowerDiagMatrix<SubType2_, 1> & other)
      {
        this->first().convert(other.first());
      }

      template <typename SubType2_>
      void convert_reverse(PowerDiagMatrix<SubType2_, 1> & other) const
      {
        this->first().convert_reverse(other.first());
      }

      /**
       * \brief PowerDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      friend bool operator== (const PowerDiagMatrix & a, const PowerDiagMatrix & b)
      {
        return (a.name() == b.name()) && (a.first() == b.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP

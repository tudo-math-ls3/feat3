// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/meta_element.hpp>
#include <kernel/lafem/container.hpp>


#include <fstream>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Power-Row-Matrix meta class template
     *
     * This class template implements a horizontal composition of \e n sub-matrices of the same class.
     * This can be interpreted as a dense 1-by-m matrix of other matrices.
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
    class PowerRowMatrix
    {
      // Note: the case = 1 is specialized below
      static_assert(blocks_ > 1, "invalid block size");

      // declare this class template as a friend for recursive inheritance
      template<typename, int>
      friend class PowerRowMatrix;

      /// rest-class typedef
      typedef PowerRowMatrix<SubType_, blocks_-1> RestClass;

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
      typedef typename SubMatrixType::VectorTypeL VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, blocks_> VectorTypeR;
      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = PowerRowMatrix<typename SubType_::template ContainerType<DT2_, IT2_>, blocks_>;

      /// number of row blocks (vertical size)
      static constexpr int num_row_blocks = 1;
      /// number of column blocks (horizontal size)
      static constexpr int num_col_blocks = blocks_;

    protected:
      /// the first sub-matrix
      SubMatrixType _first;
      /// the remaining part
      RestClass _rest;

      /// base-class constructor; this one is protected for a reason
      explicit PowerRowMatrix(SubMatrixType&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

    public:
      /// default ctor
      PowerRowMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerRowMatrix(const SparseLayout<IndexType, layout_id>& layout) :
        _first(layout),
        _rest(layout)
      {
      }

      /// move ctor
      PowerRowMatrix(PowerRowMatrix&& other) :
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
       * Creates a power-row-point-matrix based on the source file.
       */
      explicit PowerRowMatrix(FileMode mode, String filename)
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
       * Creates a power-row-matrix based on the source filestream.
       */
      explicit PowerRowMatrix(FileMode mode, std::istream& file, String directory = "")
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
        if (line.find("%%MatrixMarket powerrowmatrix coordinate real general") == String::npos)
          XABORTM("Input-file is not a complatible file");

        PowerRowMatrix other(mode, file, directory);

        _first = std::move(other._first);
        _rest = std::move(other._rest);

        file.close();
      }

      /// move-assign operator
      PowerRowMatrix& operator=(PowerRowMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerRowMatrix(const PowerRowMatrix&) = delete;
      /// deleted copy-assign operator
      PowerRowMatrix& operator=(const PowerRowMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerRowMatrix()
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

        file << "%%MatrixMarket powerrowmatrix coordinate real general\n";
        for (Index i(1); i <= blocks_; ++i)
        {
          file << filename << "_pr" << i << suffix << "\n";
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
        _first.write_out(mode, directory + prefix + "_pr" + stringify(length + 1 - blocks_) + suffix);
        _rest.write_out_submatrices(mode, directory, prefix, suffix, length);
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      PowerRowMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return PowerRowMatrix(_first.clone(mode), _rest.clone(mode));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
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
        static_assert(i_ == 0, "invalid sub-matrix index");
        static_assert((0 <= j_) && (j_ < blocks_), "invalid sub-matrix index");
        return PowerElement<j_, SubMatrixType>::get(*this);
      }

      /** \copydoc at() */
      template<int i_, int j_>
      const SubMatrixType& at() const
      {
        static_assert(i_ == 0, "invalid sub-matrix index");
        static_assert((0 <= j_) && (j_ < blocks_), "invalid sub-matrix index");
        return PowerElement<j_, SubMatrixType>::get(*this);
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
        XASSERTM((i == 0) && (0 <= j) && (j < blocks_), "invalid sub-matrix index");
        return (j == 0) ? _first : _rest.get(i, j-1);
      }

      /** \copydoc get() */
      const SubMatrixType& get(int i, int j) const
      {
        XASSERTM((i == 0) && (0 <= j) && (j < blocks_), "invalid sub-matrix index");
        return (j == 0) ? _first : _rest.get(i, j-1);
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
       * \returns Matrix row count if raw = false.
       * \returns Raw matrix row count if raw = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        return first().template rows<perspective_>();
      }

      /**
       * \brief Returns the total number of columns in this matrix.
       *
       * \returns Matrix column count if raw = false.
       * \returns Raw matrix column count if raw = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        return first().template columns<perspective_>() + rest().template columns<perspective_>();
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \returns Matrix non zero element count if raw = false.
       * \returns Raw matrix non zero element count if raw = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        return first().template used_elements<perspective_>() + rest().template used_elements<perspective_>();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("PowerRowMatrix<") + SubMatrixType::name() + "," + stringify(blocks_) + ">";
      }

      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return rows(perspective_) * columns(perspective_);
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
        first().apply(r, x.first());
        rest().apply(r, x.rest(), r, DataType(1));
      }

      void apply(DenseVector<DataType , IndexType>& r, const DenseVector<DataType , IndexType>& x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        DenseVector<DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        first().apply(r, x_first);
        rest().apply(r, x_rest, r, DataType(1));
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
        first().apply_transposed(r.first(), x);
        rest().apply_transposed(r.rest(), x);
      }

      void apply_transposed(DenseVector<DataType , IndexType>& r, const DenseVector<DataType , IndexType>& x) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");

        DenseVector<DataType, IndexType> r_first(r, first().columns(), 0);
        DenseVector<DataType, IndexType> r_rest(r, rest().columns(), first().columns());

        first().apply_transposed(r_first, x);
        rest().apply_transposed(r_rest, x);
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
        first().apply(r, x.first(), y, alpha);
        rest().apply(r, x.rest(), r, alpha);
      }

      void apply(DenseVector<DataType , IndexType>& r, const DenseVector<DataType , IndexType>& x,
                 const DenseVector<DataType , IndexType>& y, DataType alpha = DataType(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

        DenseVector<DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        first().apply(r, x_first, y, alpha);
        rest().apply(r, x_rest, r, alpha);
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
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(VectorTypeR& r, const VectorTypeL& x, const VectorTypeR& y, DataType alpha = DataType(1)) const
      {
        first().apply_transposed(r.first(), x, y.first(), alpha);
        rest().apply_transposed(r.rest(), x, y.rest(), alpha);
      }

      void apply_transposed(DenseVector<DataType , IndexType>& r, const DenseVector<DataType , IndexType>& x,
        const DenseVector<DataType , IndexType>& y, DataType alpha = DataType(1)) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->columns(), "Vector size of y does not match!");

        DenseVector<DataType, IndexType> r_first(r, first().columns(), 0);
        DenseVector<DataType, IndexType> r_rest(r, rest().columns(), first().columns());

        DenseVector<DataType, IndexType> y_first(y, first().columns(), 0);
        DenseVector<DataType, IndexType> y_rest(y, rest().columns(), first().columns());

        first().apply(r_first, x, y_first, alpha);
        rest().apply(r_rest, x, y_rest, alpha);
      }

      /**
      * \brief Retrieve the maximum relative difference of this matrix and another one
      * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
      *
      * \return The largest relative difference.
      */
      DataType max_rel_diff(const PowerRowMatrix& x) const
      {
        DataType max_rel_diff = Math::max(this->first().max_rel_diff(x.first()), this->rest().max_rel_diff(x.rest()));
        return max_rel_diff;
      }

      /**
       * \brief Checks if the structural layout of this matrix matches that of another matrix.
       * This excludes comparison of the actual data values.
       *
       * \param[in] x The matrix to compare this matrix to
       *
       * \returns true if the layouts match, false otherwise.
       */
      bool same_layout(const PowerRowMatrix& x) const
      {
        return (this->name() == x.name()) && (this->first().same_layout(x.first())) && (this->rest().same_layout(x.rest()));
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return first().create_vector_l();
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(first().create_vector_r(), rest().create_vector_r());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        return this->first().get_length_of_line(row) + this->rest().get_length_of_line(row);
      }

      /// \cond internal
      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const Index length_of_base(this->first().get_length_of_line(row));

        this->first().set_line(row, pval_set, pcol_set, col_start, stride);
        this->rest().set_line(row, pval_set + stride * length_of_base, pcol_set + stride * length_of_base, col_start + this->first().template columns<Perspective::pod>(), stride);
      }

      void set_line_reverse(const Index row, const DataType * const pval_set, const Index stride = 1)
      {
        const Index length_of_base(this->first().get_length_of_line(row));

        this->first().set_line_reverse(row, pval_set, stride);
        this->rest().set_line_reverse(row, pval_set + stride * length_of_base, stride);
      }
      /// \endcond

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
          data[old_size +i] = csize[i];
        }

        return sizeof(std::uint64_t) + ireal_size + _rest.set_checkpoint_data(data, config); //generate and add checkpoint data for the _rest
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename SubType2_>
      void convert(const PowerRowMatrix<SubType2_, blocks_> & other)
      {
        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }

      template <typename SubType2_>
      void convert_reverse(PowerRowMatrix<SubType2_, blocks_> & other) const
      {
        this->first().convert_reverse(other.first());
        this->rest().convert_reverse(other.rest());
      }
    };

    /// \cond internal
    template<typename SubType_>
    class PowerRowMatrix<SubType_, 1>
    {
      template<typename, int>
      friend class PowerRowMatrix;

    public:
      typedef SubType_ SubMatrixType;
      typedef typename SubMatrixType::DataType DataType;
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutId layout_id = SubMatrixType::layout_id;
      /// Compatible L-vector type
      typedef typename SubMatrixType::VectorTypeL VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, 1> VectorTypeR;
      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = PowerRowMatrix<typename SubType_::template ContainerType<DT2_, IT2_>, 1>;

      static constexpr int num_row_blocks = 1;
      static constexpr int num_col_blocks = 1;

    protected:
      SubMatrixType _first;

      /// base-class constructor; this one is protected for a reason
      explicit PowerRowMatrix(SubMatrixType&& the_first) :
        _first(std::move(the_first))
      {
      }

    public:
      /// default ctor
      PowerRowMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerRowMatrix(const SparseLayout<IndexType, layout_id>& layout) :
        _first(layout)
      {
      }

      /// move ctor
      PowerRowMatrix(PowerRowMatrix&& other) :
        _first(std::move(other._first))
      {
      }

      /// file-input ctor
      explicit PowerRowMatrix(FileMode mode, const String& filename)
      {
        read_from(mode, filename);
      }

      /// filestream-input ctor
      explicit PowerRowMatrix(FileMode mode, std::istream& file, const String& directory = "")
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

        PowerRowMatrix other(mode, file, directory);

        _first = std::move(other._first);

        file.close();
      }

      /// move-assign operator
      PowerRowMatrix& operator=(PowerRowMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerRowMatrix(const PowerRowMatrix&) = delete;
      /// deleted copy-assign operator
      PowerRowMatrix& operator=(const PowerRowMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerRowMatrix()
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

        file << "%%MatrixMarket powerrowmatrix coordinate real general\n";
        file << filename << "_pr" << 1 << suffix << "\n";

        file.close();

        this->write_out_submatrices(mode, directory, filename, suffix);
      }

      void write_out_submatrices(FileMode mode, const String& directory, const String& prefix, const String& suffix, Index length = 1) const
      {
        _first.write_out(mode, directory + prefix + "_pr" + stringify(length) + suffix);
      }

      PowerRowMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return PowerRowMatrix(_first.clone(mode));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes();
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
        return String("PowerRowMatrix<") + SubMatrixType::name() + "," + stringify(1) + ">";
      }

      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return rows(perspective_) * columns(perspective_);
      }

      void format(DataType value = DataType(0))
      {
        first().format(value);
      }

      void clear()
      {
        first().clear();
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().apply(r, x.first());
      }

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        first().apply(r, x);
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r, x.first(), y, alpha);
      }

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
                 const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        first().apply(r, x, y, alpha);
      }

      void apply_transposed(VectorTypeR& r, const VectorTypeL& x) const
      {
        first().apply_transposed(r.first(), x);
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        first().apply_transposed(r, x);
      }

      void apply_transposed(VectorTypeR& r, const VectorTypeL& x, const VectorTypeR& y, DataType alpha = DataType(1)) const
      {
        first().apply_transposed(r.first(), x, y.first(), alpha);
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
        const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        first().apply_transposed(r, x, y, alpha);
      }

      /**
      * \brief Retrieve the maximum relative difference of this matrix and another one
      * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
      *
      * \return The largest relative difference.
      */
      DataType max_rel_diff(const PowerRowMatrix& x) const
      {
        return this->first().max_rel_diff(x.first());
      }

      /**
       * \brief Checks if the structural layout of this matrix matches that of another matrix.
       * This excludes comparison of the actual data values.
       *
       * \param[in] x The matrix to compare this matrix to
       *
       * \returns true if the layouts match, false otherwise.
       */
      bool same_layout(const PowerRowMatrix& x) const
      {
        return (this->name() == x.name()) && (this->first().same_layout(x.first()));
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

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename SubType2_>
      void convert(const PowerRowMatrix<SubType2_, 1> & other)
      {
        this->first().convert(other.first());
      }

      template <typename SubType2_>
      void convert_reverse(PowerRowMatrix<SubType2_, 1> & other) const
      {
        this->first().convert_reverse(other.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/meta_element.hpp>


#include <fstream>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Tuple-Diag-Matrix meta class template
     *
     * This class template implements a diagonal composition of \e n sub-matrices of arbitrary classes.
     * This can be interpreted as a diagonal m-by-m matrix of other matrices.
     *
     * \tparam First_, ...Rest_
     * A sequence of (meta) matrix classes which are to be composed.
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleDiagMatrix
    {
      // declare this class template as a friend for recursive inheritance
      template<typename, typename...>
      friend class TupleDiagMatrix;

      /// rest-class typedef
      typedef TupleDiagMatrix<Rest_...> RestClass;

    public:
      /// sub-matrix data type
      typedef typename First_::DataType DataType;
      /// sub-matrix index type
      typedef typename First_::IndexType IndexType;

      /// Compatible L-vector type
      typedef TupleVector<typename First_::VectorTypeL, typename Rest_::VectorTypeL...> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR, typename Rest_::VectorTypeR...> VectorTypeR;

      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleDiagMatrix<
        typename First_::template ContainerType<DT2_, IT2_>,
        typename Rest_::template ContainerType<DT2_, IT2_>...>;

      /// this typedef lets you create a matrix container with different Data and Index types
      template <typename DT2_, typename IT2_>
      using ContainerTypeByDI = TupleDiagMatrix<
        typename First_::template ContainerType<DT2_, IT2_>,
        typename Rest_::template ContainerType<DT2_, IT2_>...>;

      /// number of row blocks (vertical size)
      static constexpr int num_row_blocks = RestClass::num_row_blocks + 1;
      /// number of column blocks (horizontal size)
      static constexpr int num_col_blocks = RestClass::num_col_blocks + 1;

    protected:
      /// the first sub-matrix
      First_ _first;
      /// the remaining part
      RestClass _rest;

      /// base-class constructor; this one is protected for a reason
      explicit TupleDiagMatrix(First_&& the_first, RestClass&& the_rest) :
        _first(std::forward<First_>(the_first)),
        _rest(std::forward<RestClass>(the_rest))
      {
      }

      /// Returns a list of all sub-matrix type names
      static String sub_name_list()
      {
        return First_::name() + "," + RestClass::sub_name_list();
      }

    public:
      /// default ctor
      TupleDiagMatrix()
      {
      }

      /// Sub-Matrix emplacement constructor
      explicit TupleDiagMatrix(First_&& the_first, Rest_&&... the_rest) :
        _first(std::forward<First_>(the_first)),
        _rest(std::forward<Rest_>(the_rest)...)
      {
      }

      /// move ctor
      TupleDiagMatrix(TupleDiagMatrix&& other) :
        _first(std::forward<First_>(other._first)),
        _rest(std::forward<RestClass>(other._rest))
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a tuple-diag-point-matrix based on the source file.
       */
      explicit TupleDiagMatrix(FileMode mode, const String& filename)
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
       * Creates a tuple-diag-matrix based on the source filestream.
       */
      explicit TupleDiagMatrix(FileMode mode, std::istream& file, const String& directory = "")
      {
        String line;
        do {
          if (file.eof())
            XABORTM("Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        _first = First_(mode, directory + line);
        _rest = RestClass(mode, file, directory);
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
        if (line.find("%%MatrixMarket tuplediagmatrix coordinate real general") == String::npos)
          XABORTM("Input-file is not a complatible file");

        TupleDiagMatrix other(mode, file, directory);

        _first = std::move(other._first);
        _rest = std::move(other._rest);

        file.close();
      }

      /// move-assign operator
      TupleDiagMatrix& operator=(TupleDiagMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::forward<First_>(other._first);
          _rest = std::forward<RestClass>(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleDiagMatrix(const TupleDiagMatrix&) = delete;
      /// deleted copy-assign operator
      TupleDiagMatrix& operator=(const TupleDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~TupleDiagMatrix()
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

        file << "%%MatrixMarket tuplediagmatrix coordinate real general" << "\n";
        for (Index i(1); i <= num_row_blocks; ++i)
        {
          file << filename << "_td" << i << suffix << "\n";
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
      void write_out_submatrices(FileMode mode, const String& directory, const String& prefix, const String& suffix, Index length = num_row_blocks) const
      {
        _first.write_out(mode, directory + prefix + "_td" + stringify(length + 1 - num_row_blocks) + suffix);
        _rest.write_out_submatrices(mode, directory, prefix, suffix, length);
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      TupleDiagMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleDiagMatrix(_first.clone(mode), _rest.clone(mode));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
      }

      /**
       * \brief Creates and returns the transpose of this matrix.
       */
      TupleDiagMatrix transpose() const
      {
        return TupleDiagMatrix(_first.transpose(), _rest.transpose());
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
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<int i_, int j_>
      const typename TupleElement<i_, First_, Rest_...>::Type& at() const
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /// \cond internal
      First_& first()
      {
        return _first;
      }

      const First_& first() const
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

      /// Returns the total number of rows in this matrix.
      Index rows() const
      {
        return first().rows() + rest().rows();
      }

      /// Returns the total number of columns in this matrix.
      Index columns() const
      {
        return first().columns() + rest().columns();
      }

      /// Returns the total number of non-zeros in this matrix.
      Index used_elements() const
      {
        return first().used_elements() + rest().used_elements();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleDiagMatrix<") + sub_name_list() + ">";
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

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
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
        first().apply_transposed(r.first(), x.first());
        rest().apply_transposed(r.rest(), x.rest());
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");

        DenseVector<DataType, IndexType> r_first(r, first().columns(), 0);
        DenseVector<DataType, IndexType> r_rest(r, rest().columns(), first().columns());

        DenseVector<DataType, IndexType> x_first(x, first().rows(), 0);
        DenseVector<DataType, IndexType> x_rest(x, rest().rows(), first().rows());

        first().apply_transposed(r_first, x_first);
        rest().apply_transposed(r_rest, x_rest);
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
        first().apply_transposed(r.first(), x.first(), y.first(), alpha);
        rest().apply_transposed(r.rest(), x.rest(), y.rest(), alpha);
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
        const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->columns(), "Vector size of y does not match!");

        DenseVector<DataType, IndexType> r_first(r, first().columns(), 0);
        DenseVector<DataType, IndexType> r_rest(r, rest().columns(), first().columns());

        DenseVector<DataType, IndexType> x_first(x, first().rows(), 0);
        DenseVector<DataType, IndexType> x_rest(x, rest().rows(), first().rows());

        DenseVector<DataType, IndexType> y_first(y, first().columns(), 0);
        DenseVector<DataType, IndexType> y_rest(y, rest().columns(), first().columns());

        first().apply_transposed(r_first, x_first, y_first, alpha);
        rest().apply_transposed(r_rest, x_rest, y_rest, alpha);
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

      void scale_rows(const TupleDiagMatrix& a, const VectorTypeL& w)
      {
        first().scale_rows(a.first(), w.first());
        rest().scale_rows(a.rest(), w.rest());
      }

      void scale_cols(const TupleDiagMatrix& a, const VectorTypeR& w)
      {
        first().scale_cols(a.first(), w.first());
        rest().scale_cols(a.rest(), w.rest());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        const Index brows(this->first().rows());

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
        const Index brows(this->first().rows());
        const Index bcolumns(this->first().columns());

        if (row < brows)
        {
          this->first().set_line(row, pval_set, pcol_set, col_start, stride);
        }
        else
        {
          this->rest().set_line(row - brows, pval_set, pcol_set, col_start + bcolumns, stride);
        }
      }

      void set_line_reverse(const Index row, DataType * const pval_set, const Index stride = 1)
      {
        const Index brows(this->first().rows());

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
      template <typename First2_, typename... Rest2_>
      void convert(const TupleDiagMatrix<First2_, Rest2_...>& other)
      {
        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }

      template <typename First2_, typename... Rest2_>
      void convert_reverse(TupleDiagMatrix<First2_, Rest2_...>& other) const
      {
        this->first().convert_reverse(other.first());
        this->rest().convert_reverse(other.rest());
      }

      /**
       * \brief TupleDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       *
       */
      friend bool operator== (const TupleDiagMatrix & a, const TupleDiagMatrix & b)
      {
        return (a.name() == b.name()) && (a.first() == b.first()) && (a.rest() == b.rest());
      }
    };

    /// \cond internal
    template<typename First_>
    class TupleDiagMatrix<First_>
    {
      template<typename,typename...>
      friend class TupleDiagMatrix;

    public:
      typedef typename First_::DataType DataType;
      typedef typename First_::IndexType IndexType;
      /// Compatible L-vector type
      typedef TupleVector<typename First_::VectorTypeL> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR> VectorTypeR;

      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleDiagMatrix<typename First_::template ContainerType<DT2_, IT2_> >;

      static constexpr int num_row_blocks = 1;
      static constexpr int num_col_blocks = 1;

    protected:
      First_ _first;

      static String sub_name_list()
      {
        return First_::name();
      }

    public:
      /// default ctor
      TupleDiagMatrix()
      {
      }

      explicit TupleDiagMatrix(First_&& the_first) :
        _first(std::forward<First_>(the_first))
      {
      }

      /// move ctor
      TupleDiagMatrix(TupleDiagMatrix&& other) :
        _first(std::forward<First_>(other._first))
      {
      }

      /// file-input ctor
      explicit TupleDiagMatrix(FileMode mode, const String& filename)
      {
        read_from(mode, filename);
      }

      /// filestream-input ctor
      explicit TupleDiagMatrix(FileMode mode, std::istream& file, const String& directory = "")
      {
        String line;
        do {
          if (file.eof())
            XABORTM("Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        _first = First_(mode, directory + line);
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

        TupleDiagMatrix other(mode, file, directory);

        _first = std::move(other._first);

        file.close();
      }

      /// move-assign operator
      TupleDiagMatrix& operator=(TupleDiagMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::forward<First_>(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleDiagMatrix(const TupleDiagMatrix&) = delete;
      /// deleted copy-assign operator
      TupleDiagMatrix& operator=(const TupleDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~TupleDiagMatrix()
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

        file << "%%MatrixMarket tuplediagmatrix coordinate real general" << "\n";
        file << filename << "_td" << 1 << suffix << "\n";

        file.close();

        this->write_out_submatrices(mode, directory, filename, suffix);
      }

      void write_out_submatrices(FileMode mode, const String& directory, const String& prefix, const String& suffix, Index length = 1) const
      {
        _first.write_out(mode, directory + prefix + "_td" + stringify(length) + suffix);
      }

      TupleDiagMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleDiagMatrix(_first.clone(mode));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes();
      }

      TupleDiagMatrix transpose() const
      {
        return TupleDiagMatrix(_first.transpose());
      }

      template<int i, int j>
      typename TupleElement<i, First_>::Type& at()
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _first;
      }

      template<int i, int j>
      const typename TupleElement<i, First_>::Type& at() const
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _first;
      }

      First_& first()
      {
        return _first;
      }

      const First_& first() const
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

      Index rows() const
      {
        return first().rows();
      }

      Index columns() const
      {
        return first().columns();
      }

      Index used_elements() const
      {
        return first().used_elements();
      }

      void format(DataType value = DataType(0))
      {
        first().format(value);
      }

      void clear()
      {
        first().clear();
      }

      Index size() const
      {
        return rows() * columns();
      }

      void extract_diag(VectorTypeL& diag) const
      {
        first().extract_diag(diag.first());
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().apply(r.first(), x.first());
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r.first(), x.first(), y.first(), alpha);
      }

      void apply_transposed(VectorTypeR& r, const VectorTypeL& x) const
      {
        first().apply_transposed(r.first(), x.first());
      }

      void apply_transposed(VectorTypeR& r, const VectorTypeL& x, const VectorTypeR& y, DataType alpha = DataType(1)) const
      {
        first().apply_transposed(r.first(), x.first(), y.first(), alpha);
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

      void scale_rows(const TupleDiagMatrix& a, const VectorTypeL& w)
      {
        first().scale_rows(a.first(), w.first());
      }

      void scale_cols(const TupleDiagMatrix& a, const VectorTypeR& w)
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

      void set_line_reverse(const Index row, DataType * const pval_set, const Index stride = 1) const
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
      template <typename First2_>
      void convert(const TupleDiagMatrix<First2_>& other)
      {
        this->first().convert(other.first());
      }

      template <typename First2_>
      void convert_reverse(TupleDiagMatrix<First2_>& other) const
      {
        this->first().convert_reverse(other.first());
      }

      /**
       * \brief TupleDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      friend bool operator== (const TupleDiagMatrix & a, const TupleDiagMatrix & b)
      {
        return (a.name() == b.name()) && (a.first() == b.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

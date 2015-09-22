#pragma once
#ifndef KERNEL_LAFEM_TUPLE_DIAG_MATRIX_HPP
#define KERNEL_LAFEM_TUPLE_DIAG_MATRIX_HPP 1

// includes, FEAST
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/meta_element.hpp>

#include <fstream>

namespace FEAST
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
      /// sub-matrix memory type
      typedef typename First_::MemType MemType;
      /// sub-matrix data type
      typedef typename First_::DataType DataType;
      /// sub-matrix index type
      typedef typename First_::IndexType IndexType;

      /// Compatible L-vector type
      typedef TupleVector<typename First_::VectorTypeL, typename Rest_::VectorTypeL...> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR, typename Rest_::VectorTypeR...> VectorTypeR;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleDiagMatrix<
        typename First_::template ContainerType<Mem2_, DT2_, IT2_>,
        typename Rest_::template ContainerType<Mem2_, DT2_, IT2_>...>;

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
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
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
        _first(std::move(the_first)),
        _rest(std::move(the_rest...))
      {
      }

      /// move ctor
      TupleDiagMatrix(TupleDiagMatrix&& other) :
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
       * Creates a tuple-diag-point-matrix based on the source file.
       */
      explicit TupleDiagMatrix(FileMode mode, String filename)
      {
        CONTEXT("When creating TupleDiagMatrix");

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
      explicit TupleDiagMatrix(FileMode mode, std::istream& file, String directory = "")
      {
        CONTEXT("When creating TupleDiagMatrix");

        String line;
        do {
          if (file.eof())
            throw InternalError(__func__, __FILE__, __LINE__, "Wrong Input-file");
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
      void read_from(FileMode mode, String filename)
      {
        CONTEXT("When reading in TupleDiagMatrix");

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
        if (line.find("%%MatrixMarket tuplediagmatrix coordinate real general") == String::npos)
          throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a complatible file");

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
          _first = std::move(other._first);
          _rest = std::move(other._rest);
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
        CONTEXT("When writing out TupleDiagMatrix");

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

        file << "%%MatrixMarket tuplediagmatrix coordinate real general" << std::endl;
        for (Index i(1); i <= num_row_blocks; ++i)
        {
          file << filename << "_td" << i << suffix << std::endl;
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
      void write_out_submatrices(FileMode mode, String directory, String prefix, String suffix, Index length = num_row_blocks) const
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
       * \brief Clears this matrix.
       *
       * \param[in] value
       * The value to which the matrix' entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
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

      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, rest().columns(), first().columns());

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

      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x,
                 const DenseVector<MemType, DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        DenseVector<MemType, DataType, IndexType> y_first(y, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> y_rest(y, rest().rows(), first().rows());

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

      void scale_rows(const TupleDiagMatrix& a, const VectorTypeL& w)
      {
        first().scale_rows(a.first(), w.first());
        rest().scale_rows(a.rest(), w.rest());
      }

      void scale_cols(const TupleDiagMatrix& a, const VectorTypeL& w)
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
      /// \endcond

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename First2_, typename... Rest2_>
      void convert(const TupleDiagMatrix<First2_, Rest2_...>& other)
      {
        CONTEXT("When converting TupleDiagMatrix");

        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }

      /**
       * \brief TupleDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       *
       * \compilerhack Intel C++ 14 variadic template bug
       */
#if defined(FEAST_COMPILER_INTEL) && (FEAST_COMPILER_INTEL < 1500)
      template<typename F2_, typename... R2_>
      friend bool operator== (const TupleDiagMatrix & a, const TupleDiagMatrix<F2_, R2_...>& b)
#else
      template <typename Mem2_>
      friend bool operator== (const TupleDiagMatrix & a, const ContainerType<Mem2_> & b)
#endif
      {
        CONTEXT("When comparing TupleDiagMatrices");

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
      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;
      typedef typename First_::IndexType IndexType;
      /// Compatible L-vector type
      typedef TupleVector<typename First_::VectorTypeL> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR> VectorTypeR;

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class TupleDiagMatrix<typename First_::template ContainerType<Mem2_, DT2_, IT2_> >;

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
        _first(std::move(the_first))
      {
      }

      /// move ctor
      TupleDiagMatrix(TupleDiagMatrix&& other) :
        _first(std::move(other._first))
      {
      }

      /// file-input ctor
      explicit TupleDiagMatrix(FileMode mode, String filename)
      {
        CONTEXT("When creating TupleDiagMatrix");

        read_from(mode, filename);
      }

      /// filestream-input ctor
      explicit TupleDiagMatrix(FileMode mode, std::istream& file, String directory = "")
      {
        CONTEXT("When creating TupleDiagMatrix");

        String line;
        do {
          if (file.eof())
            throw InternalError(__func__, __FILE__, __LINE__, "Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        _first = First_(mode, directory + line);
      }

      void read_from(FileMode mode, String filename)
      {
        CONTEXT("When reading in TupleDiagMatrix");

        String directory;
        auto found = filename.rfind("/");
        if (found != std::string::npos)
        {
          directory = filename.substr(0, found + 1);
        }

        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);

        TupleDiagMatrix other(mode, file, directory);

        _first = std::move(other._first);

        file.close();
      }

      /// move-assign operator
      TupleDiagMatrix& operator=(TupleDiagMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
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
        CONTEXT("When writing out TupleDiagMatrix");

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

        file << "%%MatrixMarket tuplediagmatrix coordinate real general" << std::endl;
        file << filename << "_td" << 1 << suffix << std::endl;

        file.close();

        this->write_out_submatrices(mode, directory, filename, suffix);
      }

      void write_out_submatrices(FileMode mode, String directory, String prefix, String suffix, Index length = 1) const
      {
        _first.write_out(mode, directory + prefix + "_td" + stringify(length) + suffix);
      }

      TupleDiagMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleDiagMatrix(_first.clone(mode));
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

      void scale_cols(const TupleDiagMatrix& a, const VectorTypeL& w)
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

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename First2_>
      void convert(const TupleDiagMatrix<First2_>& other)
      {
        CONTEXT("When converting TupleDiagMatrix");

        this->first().convert(other.first());
      }

      /**
       * \brief TupleDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_>
      friend bool operator== (const TupleDiagMatrix & a, const ContainerType<Mem2_> & b)
      {
        CONTEXT("When comparing TupleDiagMatrices");

        return (a.name() == b.name()) && (a.first() == b.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_DIAG_MATRIX_HPP

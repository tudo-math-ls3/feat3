#pragma once
#ifndef KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP
#define KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP 1

// includes, FEAST
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/meta_element.hpp>

namespace FEAST
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
      Index blocks_>
    class PowerDiagMatrix
    {
      // declare this class template as a friend for recursive inheritance
      template<typename, Index>
      friend class PowerDiagMatrix;

      /// rest-class typedef
      typedef PowerDiagMatrix<SubType_, blocks_-1> RestClass;

    public:
      /// sub-matrix type
      typedef SubType_ SubMatrixType;
      /// sub-matrix memory type
      typedef typename SubMatrixType::MemType MemType;
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
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class PowerDiagMatrix<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, blocks_>;

        /// number of row blocks (vertical size)
      static constexpr Index num_row_blocks = blocks_;
        /// number of column blocks (horizontal size)
      static constexpr Index num_col_blocks = blocks_;

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
      explicit PowerDiagMatrix(const SparseLayout<MemType, IndexType, layout_id>& layout) :
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
        if (line.find("%%MatrixMarket powerdiagmatrix coordinate real general") == String::npos)
          throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a complatible file");

        PowerDiagMatrix other(mode, file, directory);

        _first = std::move(other._first);
        _rest = std::move(other._rest);

        file.close();
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a power-diag-matrix based on the source filestream.
       */
      explicit PowerDiagMatrix(FileMode mode, std::istream& file, String directory = "")
      {
        CONTEXT("When creating PowerDiagMatrix");

        String line;
        do {
          if (file.eof())
            throw InternalError(__func__, __FILE__, __LINE__, "Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        SubMatrixType tmp_first(mode, directory + line);
        _first = std::move(tmp_first);

        RestClass tmp_rest(mode, file, directory);
        _rest = std::move(tmp_rest);
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
        CONTEXT("When writing out PowerDiagMatrix");

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

        file << "%%MatrixMarket powerdiagmatrix coordinate real general" << std::endl;
        for (Index i(1); i <= blocks_; ++i)
        {
          file << filename << "_pd" << i << suffix << std::endl;
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
      PowerDiagMatrix clone() const
      {
        return PowerDiagMatrix(_first.clone(), _rest.clone());
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
      template<Index i_, Index j_>
      SubMatrixType& at()
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        static_assert(i_ < blocks_, "invalid sub-matrix index");
        return PowerElement<i_, SubMatrixType>::get(*this);
      }

      /** \copydoc at() */
      template<Index i_, Index j_>
      const SubMatrixType& at() const
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        static_assert(i_ < blocks_, "invalid sub-matrix index");
        return PowerElement<i_, SubMatrixType>::get(*this);
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

      Index row_blocks() const
      {
        return Index(num_row_blocks);
      }

      Index col_blocks() const
      {
        return Index(num_col_blocks);
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
        return String("PowerDiagMatrix<") + SubMatrixType::name() + "," + stringify(blocks_) + ">";
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
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().template apply<Algo_>(r.first(), x.first());
        rest().template apply<Algo_>(r.rest(), x.rest());
      }

      template<typename Algo_>
      void apply(DenseVector<MemType, DataType , IndexType>& r, const DenseVector<MemType, DataType , IndexType>& x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        first().template apply<Algo_>(r_first, x_first);
        rest().template apply<Algo_>(r_rest, x_rest);
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
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().template apply<Algo_>(r.first(), x.first(), y.first(), alpha);
        rest().template apply<Algo_>(r.rest(), x.rest(), y.rest(), alpha);
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

        DenseVector<MemType, DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        DenseVector<MemType, DataType, IndexType> y_first(y, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> y_rest(y, rest().rows(), first().rows());

        first().template apply<Algo_>(r_first, x_first, y_first, alpha);
        rest().template apply<Algo_>(r_rest, x_rest, y_rest, alpha);
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

      template<typename Algo_>
      void scale_rows(const PowerDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_rows<Algo_>(a.first(), w.first());
        rest().template scale_rows<Algo_>(a.rest(), w.rest());
      }

      template<typename Algo_>
      void scale_cols(const PowerDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_cols<Algo_>(a.first(), w.first());
        rest().template scale_cols<Algo_>(a.rest(), w.rest());
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
#ifdef FEAST_COMPILER_MICROSOFT
      template <typename SubType_>
      void convert(const PowerDiagMatrix<SubType_, blocks_>& other)
#else
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const ContainerType<Mem2_, DT2_, IT2_> & other)
#endif
      {
        CONTEXT("When converting PowerDiagMatrix");

        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }

      /**
       * \brief PowerDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
#ifdef FEAST_COMPILER_MICROSOFT
      template <typename SubType_>
      friend bool operator== (const PowerDiagMatrix & a, const PowerDiagMatrix<SubType_, blocks_> & b)
#else
      template <typename Mem2_>
      friend bool operator== (const PowerDiagMatrix & a, const ContainerType<Mem2_> & b)
#endif
      {
        CONTEXT("When comparing PowerDiagMatrices");

        return (a.name() == b.name()) && (a.first() == b.first()) && (a.rest() == b.rest());
      }
    };

    /// \cond internal
    template<typename SubType_>
    class PowerDiagMatrix<SubType_, 1>
    {
      template<typename, Index>
      friend class PowerDiagMatrix;

    public:
      typedef SubType_ SubMatrixType;
      typedef typename SubMatrixType::MemType MemType;
      typedef typename SubMatrixType::DataType DataType;
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutId layout_id = SubMatrixType::layout_id;
      /// Compatible L-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeL, 1> VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, 1> VectorTypeR;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class PowerDiagMatrix<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, 1>;

      static constexpr Index num_row_blocks = 1;
      static constexpr Index num_col_blocks = 1;

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
      explicit PowerDiagMatrix(const SparseLayout<MemType, IndexType, layout_id>& layout) :
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
      explicit PowerDiagMatrix(FileMode mode, String filename)
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

        PowerDiagMatrix other(mode, file, directory);

        _first = std::move(other._first);

        file.close();
      }

      /// filestream-input ctor
      explicit PowerDiagMatrix(FileMode mode, std::istream& file, String directory = "")
      {
        CONTEXT("When creating PowerDiagMatrix");

        String line;
        do {
          if (file.eof())
            throw InternalError(__func__, __FILE__, __LINE__, "Wrong Input-file");
          std::getline(file, line);
          line.trim_me();
        } while (line.find("%%") == 0 || line == "");

        SubMatrixType tmp_first(mode, directory + line);
        _first = std::move(tmp_first);
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
        CONTEXT("When writing out PowerDiagMatrix");

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

        file << "%%MatrixMarket powerdiagmatrix coordinate real general" << std::endl;
        file << filename << "_pd" << 1 << suffix << std::endl;

        file.close();

        this->write_out_submatrices(mode, directory, filename, suffix);
      }

      void write_out_submatrices(FileMode mode, String directory, String prefix, String suffix, Index length = 1) const
      {
        _first.write_out(mode, directory + prefix + "_pd" + stringify(length) + suffix);
      }

      PowerDiagMatrix clone() const
      {
        return PowerDiagMatrix(_first.clone());
      }

      PowerDiagMatrix transpose() const
      {
        return PowerDiagMatrix(_first.transpose());
      }

      template<Index i, Index j>
      SubMatrixType& at()
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _first;
      }

      template<Index i, Index j>
      const SubMatrixType& at() const
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
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

      Index row_blocks() const
      {
        return Index(1);
      }

      Index col_blocks() const
      {
        return Index(1);
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

      static String name()
      {
        return String("PowerDiagMatrix<") + SubMatrixType::name() + "," + stringify(1) + ">";
      }

      void format(DataType value = DataType(0))
      {
        first().format(value);
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().template apply<Algo_>(r.first(), x.first());
      }

      template<typename Algo_>
      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x) const
      {
        first().template apply<Algo_>(r, x);
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().template apply<Algo_>(r.first(), x.first(), y.first(), alpha);
      }

      template<typename Algo_>
      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x,
                 const DenseVector<MemType, DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        first().template apply<Algo_>(r, x, y, alpha);
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

      template<typename Algo_>
      void scale_rows(const PowerDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_rows<Algo_>(a.first(), w.first());
      }

      template<typename Algo_>
      void scale_cols(const PowerDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_cols<Algo_>(a.first(), w.first());
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
#ifdef FEAST_COMPILER_MICROSOFT
      template <typename SubType_>
      void convert(const PowerDiagMatrix<SubType_, Index(1)>& other)
#else
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const ContainerType<Mem2_, DT2_, IT2_> & other)
#endif
      {
        CONTEXT("When converting PowerDiagMatrix");

        this->first().convert(other.first());
      }

      /**
       * \brief PowerDiagMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
#ifdef FEAST_COMPILER_MICROSOFT
      template <typename SubType_>
      friend bool operator== (const PowerDiagMatrix & a, const PowerDiagMatrix<SubType_, Index(1)> & b)
#else
      template <typename Mem2_>
      friend bool operator== (const PowerDiagMatrix & a, const ContainerType<Mem2_> & b)
#endif
      {
        CONTEXT("When comparing PowerDiagMatrices");

        return (a.name() == b.name()) && (a.first() == b.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP

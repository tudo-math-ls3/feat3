// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/container.hpp>


#include <fstream>

namespace FEAT
{
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      // helper class for extract_diag function of PowerFullMatrix
      template<int n_>
      struct ExtDiagPFM;
    }
    /// \endcond

    /**
     * \brief Power-Full-Matrix meta class template
     *
     * This class template implements a composition of \e m x\e n sub-matrices of the same class.
     * This can be interpreted as an m-by-n matrix of other matrices.
     *
     * \tparam SubType_
     * The type of the sub-matrix.
     *
     * \tparam width_
     * The number of sub-matrix blocks per row.
     *
     * \tparam height_
     * The number of sub-matrix blocks per column.
     *
     * \author Peter Zajac
     */
    template<
      typename SubType_,
      int width_,
      int height_>
    class PowerFullMatrix
    {
      static_assert((width_ > 0) && (height_ > 0), "invalid matrix dimensions");

    public:
      /// container-class typedef
      typedef PowerColMatrix<PowerRowMatrix<SubType_, width_>, height_> ContClass;
      /// sub-matrix type
      typedef SubType_ SubMatrixType;
      /// sub-matrix data type
      typedef typename SubMatrixType::DataType DataType;
      /// sub-matrix index type
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutId layout_id = SubMatrixType::layout_id;
      /// Compatible L-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeL, height_> VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, width_> VectorTypeR;
      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = PowerFullMatrix<
        typename SubType_::template ContainerType<DT2_, IT2_>, width_, height_>;

      /// number of row blocks (vertical size)
      static constexpr int num_row_blocks = height_;
      /// number of column blocks (horizontal size)
      static constexpr int num_col_blocks = width_;

    protected:
      // the container
      ContClass _container;

    protected:
      /// base-class emplacement constructor
      explicit PowerFullMatrix(ContClass&& cont) :
        _container(std::move(cont))
      {
      }

    public:
      /// default ctor
      PowerFullMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerFullMatrix(const SparseLayout<IndexType, layout_id>& layout) :
        _container(layout)
      {
      }

      /// move ctor
      PowerFullMatrix(PowerFullMatrix&& other) :
        _container(std::move(other._container))
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a power-full-point-matrix based on the source file.
       */
      explicit PowerFullMatrix(FileMode mode, const String& filename)
      {
        ContClass other(mode, filename);
        _container = std::move(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a power-full-matrix based on the source filestream.
       */
      explicit PowerFullMatrix(FileMode mode, std::istream& file, const String& directory = "")
      {
        ContClass other(mode, file, directory);
        _container = std::move(other);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, const String& filename)
      {
        ContClass other(mode, filename);
        _container = std::move(other);
      }

      /// move-assign operator
      PowerFullMatrix& operator=(PowerFullMatrix&& other)
      {
        if(this != &other)
        {
          _container = std::move(other._container);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerFullMatrix(const PowerFullMatrix&) = delete;
      /// deleted copy-assign operator
      PowerFullMatrix& operator=(const PowerFullMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerFullMatrix()
      {
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, const String& filename) const
      {
        _container.write_out(mode, filename);
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      PowerFullMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return PowerFullMatrix(_container.clone(mode));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _container.bytes();
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
        static_assert((0 <= i_) && (i_ < height_), "invalid sub-matrix index");
        static_assert((0 <= j_) && (j_ < width_), "invalid sub-matrix index");
        return _container.template at<i_, 0>().template at<0, j_>();
      }

      /** \copydoc at() */
      template<int i_, int j_>
      const SubMatrixType& at() const
      {
        static_assert((0 <= i_) && (i_ < height_), "invalid sub-matrix index");
        static_assert((0 <= j_) && (j_ < width_), "invalid sub-matrix index");
        return _container.template at<i_, 0>().template at<0, j_>();
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
        XASSERTM((0 <= i) && (i < height_), "invalid sub-matrix row index");
        XASSERTM((0 <= j) && (j < width_), "invalid sub-matrix column index");
        return _container.get(i, 0).get(0, j);
      }

      /** \copydoc get() */
      const SubMatrixType& get(int i, int j) const
      {
        XASSERTM((0 <= i) && (i < height_), "invalid sub-matrix row index");
        XASSERTM((0 <= j) && (j < width_), "invalid sub-matrix column index");
        return _container.get(i, 0).get(0, j);
      }

      /// \cond internal
      int row_blocks() const
      {
        return num_row_blocks;
      }

      int col_blocks() const
      {
        return num_col_blocks;
      }

      Index get_length_of_line(const Index row) const
      {
        return _container.get_length_of_line(row);
      }

      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        _container.set_line(row, pval_set, pcol_set, col_start, stride);
      }

      void set_line_reverse(const Index row, DataType * const pval_set, const Index stride = 1)
      {
        _container.set_line_reverse(row, pval_set, stride);
      }

      ContClass& get_container()
      {
        return _container;
      }

      const ContClass& get_container() const
      {
        return _container;
      }
      /// \endcond

      VectorTypeL create_vector_l() const
      {
        return _container.create_vector_l();
      }

      VectorTypeR create_vector_r() const
      {
        return _container.create_vector_r();
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
        return _container.template rows<perspective_>();
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
        return _container.template columns<perspective_>();
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
        return _container.template used_elements<perspective_>();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("PowerFullMatrix<") + SubMatrixType::name() + "," + stringify(width_) + "," + stringify(height_) + ">";
      }

      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return rows<perspective_>() * columns<perspective_>();
      }

      void format(DataType value = DataType(0))
      {
        _container.format(value);
      }

      /// extract main diagonal vector from matrix
      void extract_diag(VectorTypeL& diag) const
      {
        static_assert(width_ == height_, "cannot extract diagonal from rectangular matrix");
        Intern::ExtDiagPFM<width_>::extract_diag(*this, diag);
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        _container.apply(r, x);
      }

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        _container.apply(r, x);
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        _container.apply(r, x, y, alpha);
      }

      void apply(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
                 const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        _container.apply(r, x, y, alpha);
      }

      void apply_transposed(VectorTypeR& r, const VectorTypeL& x) const
      {
        _container.apply_transposed(r, x);
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x) const
      {
        _container.apply_transposed(r, x);
      }

      void apply_transposed(VectorTypeR& r, const VectorTypeL& x, const VectorTypeR& y, DataType alpha = DataType(1)) const
      {
        _container.apply_transposed(r, x, y, alpha);
      }

      void apply_transposed(DenseVector<DataType, IndexType>& r, const DenseVector<DataType, IndexType>& x,
        const DenseVector<DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        _container.apply_transposed(r, x, y, alpha);
      }

      /**
      * \brief Retrieve the maximum relative difference of this matrix and another one
      * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
      *
      * \return The largest relative difference.
      */
      DataType max_rel_diff(const PowerFullMatrix& x) const
      {
        return this->_container.max_rel_diff(x._container);
      }

      /**
       * \brief Checks if the structural layout of this matrix matches that of another matrix.
       * This excludes comparison of the actual data values.
       *
       * \param[in] x The matrix to compare this matrix to
       *
       * \returns true if the layouts match, false otherwise.
       */
      bool same_layout(const PowerFullMatrix& x) const
      {
        return this->_container.same_layout(x._container);
      }

      template <typename SubType2_>
      void convert(const PowerFullMatrix<SubType2_, width_, height_> & other)
      {
        _container.convert(other._container);
      }

      template <typename SubType2_>
      void convert_reverse(PowerFullMatrix<SubType2_, width_, height_> & other) const
      {
        _container.convert_reverse(other._container);
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return _container.get_checkpoint_size(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        _container.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        return _container.set_checkpoint_data(data, config);
      }
    }; // class PowerFullMatrix

    /// \cond internal
    namespace Intern
    {
      // JacobiHelper specialization for PowerFullMatrix
      template<int n_>
      struct ExtDiagPFM
      {
        static_assert(n_ > 1, "invalid block size");
        template<typename SMT_, int m_>
        static void extract_diag(
          const PowerFullMatrix<SMT_, m_, m_>& matrix,
          PowerVector<typename SMT_::VectorTypeL, m_>& diag)
        {
          ExtDiagPFM<n_-1>::extract_diag(matrix, diag);
          matrix.template at<n_-1,n_-1>().extract_diag(diag.template at<n_-1>());
        }
      };

      template<>
      struct ExtDiagPFM<1>
      {
        template<typename SMT_, int m_>
        static void extract_diag(
          const PowerFullMatrix<SMT_, m_, m_>& matrix,
          PowerVector<typename SMT_::VectorTypeL, m_>& diag)
        {
          matrix.template at<0,0>().extract_diag(diag.template at<0>());
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

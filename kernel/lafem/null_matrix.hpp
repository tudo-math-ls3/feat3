// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/lafem/base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      template<typename DT_, typename IT_, int size_>
      struct BlockedVectorHelper
      {
        static_assert(size_ > 1, "invalid block size");
        typedef DenseVectorBlocked<DT_, IT_, size_> VectorType;
      };

      template<typename DT_, typename IT_>
      struct BlockedVectorHelper<DT_, IT_, 1>
      {
        typedef DenseVector<DT_, IT_> VectorType;
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Null Matrix implementation
     *
     * This class implements a null matrix in the common LAFEM matrix interface.
     * The main purpose of this class is to act as a "null block" inside a bigger
     * meta-matrix, which is typically a LAFEM::TupleMatrix.
     *
     * This matrix is a "hybrid-blocked-scalar" implementation, i.e. it can act
     * as both a "blocked" null-matrix class as well as a "scalar" null-matrix class
     * by setting \p block_height_ = \p block_width_ = 1. The compatible L/R vector
     * types are chosen to be \p DenseVector or \p DenseVectorBlocked depending on
     * whether the corresponding block dimension is 1 or greater than 1 (in analogy
     * to the SparseMatrixBCSR class).
     *
     * Also note that the \p DT_ and \p IT_ template arguments are only
     * required for the conformance with the LAFEM matrix interface.
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_, int block_height_ = 1, int block_width_ = 1>
    class NullMatrix
    {
      static_assert(block_height_ > 0, "invalid block size");
      static_assert(block_width_ > 0, "invalid block size");

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our block height
      static constexpr int block_height = block_height_;
      /// Our block width
      static constexpr int block_width = block_width_;
      /// Value type, meaning the type of each block
      typedef Tiny::Matrix<DataType, block_height, block_width> ValueType;

      /// ImageIterator typedef for Adjactor interface implementation
      typedef const IT_* ImageIterator;

      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = NullMatrix<DT2_, IT2_, block_height_, block_width_>;

      /// this typedef lets you create a matrix container with new Datatype and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

      /// Compatible L-vector type
      typedef typename Intern::BlockedVectorHelper<DT_, IT_, block_height_>::VectorType VectorTypeL;
      /// Compatible R-vector type
      typedef typename Intern::BlockedVectorHelper<DT_, IT_, block_width_>::VectorType VectorTypeR;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

    private:
      /// matrix dimensions
      Index _num_rows, _num_cols;

    public:

      /**
       * \brief Constructor
       *
       * Creates an empty matrix.
       */
      explicit NullMatrix() :
        _num_rows(0u), _num_cols(0u)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       *
       * \note This matrix does not allocate any memory
       */
      explicit NullMatrix(Index rows_in, Index columns_in) :
        _num_rows(rows_in), _num_cols(columns_in)
      {
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      NullMatrix(NullMatrix && other) :
        _num_rows(other._num_rows), _num_cols(other._num_cols)
      {
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      NullMatrix & operator= (NullMatrix && other)
      {
        this->_num_rows = other._num_rows;
        this->_num_cols = other._num_cols;

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a clone of this container.
       *
       * \param[in] clone_mode The actual cloning procedure.
       * \returns The created clone.
       *
       */
      NullMatrix clone(CloneMode DOXY(clone_mode) = CloneMode::Weak) const
      {
        return NullMatrix(_num_rows, _num_cols);
      }

      /** \brief Clone operation
       *
       * Create a clone of another container.
       *
       * \param[in] other The source container to create the clone from.
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      template<typename DT2_, typename IT2_>
      void clone(
        const NullMatrix<DT2_, IT2_, block_height_, block_width_> & other,
        CloneMode DOXY(clone_mode) = CloneMode::Weak)
      {
        this->_num_rows = other.num_rows();
        this->_num_cols = other.num_cols();
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename DT2_, typename IT2_>
      void convert(const NullMatrix<DT2_, IT2_, block_height_, block_width_> & other)
      {
        this->_num_rows = other.num_rows();
        this->_num_cols = other.num_cols();
      }

      /// Retrieve matrix row count
      Index num_rows() const
      {
        return this->_num_rows;
      }

      /// Retrieve matrix column count
      Index num_cols() const
      {
        return this->_num_cols;
      }

      /// Retrieve matrix non-zero element count
      Index num_nzes() const
      {
        return Index(0);
      }

      /// Retrieve matrix row count
      Index num_rows_raw() const
      {
        return this->_num_rows * Index(block_height_);
      }

      /// Retrieve matrix column count
      Index num_cols_raw() const
      {
        return this->_num_cols * Index(block_width_);
      }

      /// Retrieve matrix non-zero element count
      Index num_nzes_raw() const
      {
        return Index(0);
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "NullMatrix";
      }

      /**
       * \brief Resizes the matrix to different dimensions.
       *
       * \param[in] rows_in, columns_in
       * The new dimensions of the matrix.
       */
      void resize(Index rows_in, Index columns_in)
      {
        this->_num_rows = rows_in;
        this->_num_cols = columns_in;
      }

      /**
       * \brief Returns the total amount of bytes allocated.
       *
       * \returns 0
       */
      std::size_t bytes() const
      {
        return std::size_t(0);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const NullMatrix & DOXY(x), bool DOXY(full) = false)
      {
        // nothing to do here
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] value The value to be set (defaults to 0)
       */
      void format(const DT_ DOXY(value) = DT_(0))
      {
        // nothing to do here
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow y + \alpha~ x\f$
       *
       * \param[in] x The first summand matrix to be scaled.
       * \param[in] y The second summand matrix
       * \param[in] alpha A scalar to multiply x with.
       *
       * \warning All three matrices must have the same non zero layout.
       * This operation assumes this silently and does not check this on its own!
       */
      void axpy(
                const NullMatrix & x,
                const NullMatrix & y,
                const DT_ DOXY(alpha) = DT_(1))
      {
        XASSERTM(x.num_rows() == y.num_rows(), "Matrix rows do not match!");
        XASSERTM(x.num_rows() == this->num_rows(), "Matrix rows do not match!");
        XASSERTM(x.num_cols() == y.num_cols(), "Matrix columns do not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Matrix columns do not match!");
        XASSERTM(x.num_nzes() == y.num_nzes(), "Matrix used_elements do not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Matrix used_elements do not match!");

        // nothing to do here
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const NullMatrix & x, const DT_ DOXY(alpha))
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Row count does not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Column count does not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");

        // nothing to do here
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        return DT_(0);
      }

      /**
       * \brief Computes the 2-norm for every row
       *
       * \param[in] row_norms
       * For every row, this left-vector will contain its 2-norm
       */
      void row_norm2(VectorTypeL& row_norms) const
      {
        XASSERTM(row_norms.size() == this->num_rows(), "Matrix/Vector dimension mismatch");

        row_norms.format();
      }

      /**
       * \brief Computes the square of the 2-norm for every row
       *
       * \param[out] row_norms
       * For every row, this left-vector will contain the square of its 2-norm
       */
      void row_norm2sqr(VectorTypeL& row_norms) const
      {
        XASSERTM(row_norms.size() == this->num_rows(), "Matrix/Vector dimension mismatch");

        row_norms.format();
      }

      /**
       * \brief Calculate \f$this^\top \f$
       *
       * \return The transposed matrix
       *
       * \note The resulting matrix has transposed block dimensions, too.
       */
      NullMatrix<DT_, IT_, block_width_, block_height_> transpose() const
      {
        return NullMatrix<DT_, IT_, block_width_, block_height_>(_num_cols, _num_rows);
      }

      /**
       * \brief Calculate \f$this \leftarrow x^\top \f$
       *
       * \param[in] x The matrix to be transposed.
       */
      void transpose(const NullMatrix<DT_, IT_, block_width_, block_height_> & x)
      {
        x = this->transpose();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<DT_, IT_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->num_rows_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols_raw(), "Vector size of x does not match!");

        r.format();
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      */
      void apply_transposed(DenseVector<DT_, IT_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->num_cols_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows_raw(), "Vector size of x does not match!");

        r.format();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<DT_, IT_, block_height_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols_raw(), "Vector size of x does not match!");

        r.format();
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      */
      void apply_transposed(DenseVectorBlocked<DT_, IT_, block_width_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows_raw(), "Vector size of x does not match!");

        r.format();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<DT_, IT_> & r, const DenseVectorBlocked<DT_, IT_, block_width_> & x) const
      {
        XASSERTM(r.size() == this->num_rows_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");

        r.format();
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      */
      void apply_transposed(DenseVector<DT_, IT_> & r, const DenseVectorBlocked<DT_, IT_, block_height_> & x) const
      {
        XASSERTM(r.size() == this->num_cols_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");

        r.format();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<DT_, IT_, block_height_> & r, const DenseVectorBlocked<DT_, IT_, block_width_> & x) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");

        r.format();
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      */
      void apply_transposed(DenseVectorBlocked<DT_, IT_, block_width_> & r, const DenseVectorBlocked<DT_, IT_, block_height_> & x) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");

        r.format();
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<DT_, IT_> & r,
                 const DenseVector<DT_, IT_> & x,
                 const DenseVector<DT_, IT_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols_raw(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows_raw(), "Vector size of y does not match!");

        r.copy(y);
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(
        DenseVector<DT_, IT_> & r,
        const DenseVector<DT_, IT_> & x,
        const DenseVector<DT_, IT_> & y,
        const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_cols_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows_raw(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols_raw(), "Vector size of y does not match!");

        r.copy(y);
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<DT_, IT_, block_height_> & r,
                 const DenseVector<DT_, IT_> & x,
                 const DenseVectorBlocked<DT_, IT_, block_height_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols_raw(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows(), "Vector size of y does not match!");

        r.copy(y);
      }
      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(
        DenseVectorBlocked<DT_, IT_, block_width_> & r,
        const DenseVector<DT_, IT_> & x,
        const DenseVectorBlocked<DT_, IT_, block_width_> & y,
        const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows_raw(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols(), "Vector size of y does not match!");

        r.copy(y);
      }


      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<DT_, IT_> & r,
                 const DenseVectorBlocked<DT_, IT_, block_width_> & x,
                 const DenseVector<DT_, IT_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows_raw(), "Vector size of y does not match!");

        r.copy(y);
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(
        DenseVector<DT_, IT_> & r,
        const DenseVectorBlocked<DT_, IT_, block_height_> & x,
        const DenseVector<DT_, IT_> & y,
        const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_cols_raw(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols_raw(), "Vector size of y does not match!");

        r.copy(y);
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<DT_, IT_, block_height_> & r,
                 const DenseVectorBlocked<DT_, IT_, block_width_> & x,
                 const DenseVectorBlocked<DT_, IT_, block_height_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows(), "Vector size of y does not match!");

        r.copy(y);
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(
        DenseVectorBlocked<DT_, IT_, block_width_> & r,
        const DenseVectorBlocked<DT_, IT_, block_height_> & x,
        const DenseVectorBlocked<DT_, IT_, block_width_> & y,
        const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols(), "Vector size of y does not match!");

        r.copy(y);
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<DT_, IT_, block_height_> & r,
                 const DenseVectorBlocked<DT_, IT_, block_width_> & x,
                 const DenseVector<DT_, IT_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows_raw(), "Vector size of y does not match!");

        r.copy(y);
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(
        DenseVectorBlocked<DT_, IT_, block_width_> & r,
        const DenseVectorBlocked<DT_, IT_, block_height_> & x,
        const DenseVector<DT_, IT_> & y,
        const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols_raw(), "Vector size of y does not match!");

        r.copy(y);
      }

      ///@}

      /// \copydoc lump_rows()
      void lump_rows(VectorTypeL& lump) const
      {
        XASSERTM(lump.size() == num_rows(), "lump vector size does not match matrix row count!");
        lump.format();
      }

      /**
       * \brief Returns the lumped rows vector
       *
       * Each entry in the returned lumped rows vector contains the
       * the sum of all matrix elements in the corresponding row.
       *
       * \returns
       * The lumped vector.
       */
      VectorTypeL lump_rows() const
      {
        VectorTypeL lump = create_vector_l();
        lump_rows(lump);
        return lump;
      }


      /// \copydoc extract_diag()
      void extract_diag(VectorTypeL & diag) const
      {
        XASSERTM(diag.size() == num_rows(), "diag size does not match matrix row count!");
        XASSERTM(num_rows() == num_cols(), "matrix is not square!");
        diag.format();
      }

      /// extract main diagonal vector from matrix
      VectorTypeL extract_diag() const
      {
        VectorTypeL diag = create_vector_l();
        extract_diag(diag);
        return diag;
      }

      /// \cond internal
      // Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(this->num_rows());
      }

      // Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->num_cols());
      }

      Index row_degree(const Index) const
      {
        return Index(0);
      }

      template<typename IT2_>
      Index get_row_col_indices(const Index, IT2_ * const, const IT2_) const
      {
        return Index(0);
      }

      template<typename DT2_>
      Index get_row_values(const Index, DT2_ * const) const
      {
        return Index(0);
      }

      template<typename DT2_>
      Index set_row_values(const Index, const DT2_ * const)
      {
        return Index(0);
      }

      /// \endcond

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      uint64_t get_checkpoint_size()
      {
        return size_t(0);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & /*data*/)
      {
        // nothing to do here
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      void set_checkpoint_data(std::vector<char>& /*data*/)
      {
        // nothing to do here
      }

      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:
      /** \copydoc Adjactor::get_num_nodes_domain() */
      inline Index get_num_nodes_domain() const
      {
        return num_rows();
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      inline Index get_num_nodes_image() const
      {
        return num_cols();
      }

      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        XASSERTM(domain_node < num_rows(), "Domain node index out of range");
        return nullptr;
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        XASSERTM(domain_node < num_rows(), "Domain node index out of range");
        return nullptr;
      }
    }; // class NullMatrix
  } // namespace LAFEM
} // namespace FEAT

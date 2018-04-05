#pragma once
#ifndef KERNEL_LAFEM_NULL_MATRIX_HPP
#define KERNEL_LAFEM_NULL_MATRIX_HPP 1

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
      template<typename Mem_, typename DT_, typename IT_, int size_>
      struct BlockedVectorHelper
      {
        static_assert(size_ > 1, "invalid block size");
        typedef DenseVectorBlocked<Mem_, DT_, IT_, size_> VectorType;
      };

      template<typename Mem_, typename DT_, typename IT_>
      struct BlockedVectorHelper<Mem_, DT_, IT_, 1>
      {
        typedef DenseVector<Mem_, DT_, IT_> VectorType;
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
     * by setting \p BlockHeight_ = \p BlockWidth_ = 1. The compatible L/R vector
     * types are chosen to be \p DenseVector or \p DenseVectorBlocked depending on
     * whether the corresponding block dimension is 1 or greater than 1 (in analogy
     * to the SparseMatrixBCSR class).
     *
     * Also note that the \p Mem_, \p DT_ and \p IT_ template arguments are only
     * required for the conformance with the LAFEM matrix interface.
     *
     * \author Peter Zajac
     */
    template<typename Mem_, typename DT_, typename IT_, int BlockHeight_ = 1, int BlockWidth_ = 1>
    class NullMatrix
    {
      static_assert(BlockHeight_ > 0, "invalid block size");
      static_assert(BlockWidth_ > 0, "invalid block size");

    public:
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our block height
      static constexpr int BlockHeight = BlockHeight_;
      /// Our block width
      static constexpr int BlockWidth = BlockWidth_;
      /// Value type, meaning the type of each block
      typedef Tiny::Matrix<DataType, BlockHeight, BlockWidth> ValueType;

      /// ImageIterator typedef for Adjactor interface implementation
      typedef const IT_* ImageIterator;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = class NullMatrix<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_>;

      /// this typedef lets you create a matrix container with new Memory, Datatype and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /// Compatible L-vector type
      typedef typename Intern::BlockedVectorHelper<Mem_, DT_, IT_, BlockHeight_>::VectorType VectorTypeL;
      /// Compatible R-vector type
      typedef typename Intern::BlockedVectorHelper<Mem_, DT_, IT_, BlockWidth_>::VectorType VectorTypeR;

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
      template<typename Mem2_, typename DT2_, typename IT2_>
      void clone(
        const NullMatrix<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_> & other,
        CloneMode DOXY(clone_mode) = CloneMode::Weak)
      {
        this->_num_rows = other._num_rows;
        this->_num_cols = other._num_cols;
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const NullMatrix<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_> & other)
      {
        this->_num_rows = other._num_rows;
        this->_num_cols = other._num_cols;
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count if perspective_ = false, e.g. count every block as one row.
       * \returns Raw matrix row count if perspective_ = true, e.g. row_count * BlockHeight_.
       */
      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        if (perspective_ == Perspective::pod)
          return this->_num_rows * Index(BlockHeight_);
        else
          return this->_num_rows;
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count if perspective_ = false, e.g. count every block as one column.
       * \returns Raw matrix column count if perspective_ = true, e.g. column_count * BlockWidth_.
       */
      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        if (perspective_ == Perspective::pod)
          return this->_num_cols * Index(BlockWidth_);
        else
          return this->_num_cols;
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count if perspective_ = false, e.g. count every block as one entry.
       * \returns Raw non zero element count if perspective_ = true, e.g. used_elements * BlockHeight_ * BlockWidth_.
       */
      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
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
      void format(const DT_ DOXY(alpha) = DT_(1))
      {
        // nothing to do here
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      template <typename Mem2_>
      void copy(const NullMatrix<Mem2_, DT_, IT_, BlockHeight_, BlockWidth_> & DOXY(x), bool DOXY(full) = false)
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
       * \warning All three matrices must have the same non zero layout. This operation assumes this silently and does not check this on its own!
       */
      void axpy(
                const NullMatrix & x,
                const NullMatrix & y,
                const DT_ DOXY(alpha) = DT_(1))
      {
        XASSERTM(x.rows() == y.rows(), "Matrix rows do not match!");
        XASSERTM(x.rows() == this->rows(), "Matrix rows do not match!");
        XASSERTM(x.columns() == y.columns(), "Matrix columns do not match!");
        XASSERTM(x.columns() == this->columns(), "Matrix columns do not match!");
        XASSERTM(x.used_elements() == y.used_elements(), "Matrix used_elements do not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Matrix used_elements do not match!");

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
        XASSERTM(x.rows() == this->rows(), "Row count does not match!");
        XASSERTM(x.columns() == this->columns(), "Column count does not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Nonzero count does not match!");

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
        XASSERTM(row_norms.size() == this->rows(), "Matrix/Vector dimension mismatch");

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
        XASSERTM(row_norms.size() == this->rows(), "Matrix/Vector dimension mismatch");

        row_norms.format();
      }

      /**
       * \brief Computes the square of the 2-norm for every row, where every row is scaled by a vector
       *
       * \param[out] row_norms
       * For every (scaled) row, this left-vector will contain the square of its 2-norm
       *
       * \param[in] scal
       * The scaling vector
       *
       * This computes
       * \f[
       *    row\_norms_i = \sum_{j=0}^{n-1} scal_j (this_{ij})^2
       * \f]
       * and is used to compute
       * \f[
       *   \mathrm{tr}(B^T \mathrm{diag}(A) B)
       * \f]
       *
       */
      void row_norm2sqr(VectorTypeL& row_norms, const VectorTypeR& scal) const
      {
        XASSERTM(row_norms.size() == this->rows(), "Matrix/Vector dimension mismatch");
        XASSERTM(scal.size() == this->columns(), "Matrix/scalings dimension mismatch");

        row_norms.format();
      }

      /**
       * \brief Calculate \f$this^\top \f$
       *
       * \return The transposed matrix
       *
       * \note The resulting matrix has transposed block dimensions, too.
       */
      NullMatrix<Mem_, DT_, IT_, BlockWidth_, BlockHeight_> transpose() const
      {
        return NullMatrix<Mem_, DT_, IT_, BlockWidth_, BlockHeight_>(_num_cols, _num_rows);
      }

      /**
       * \brief Calculate \f$this \leftarrow x^\top \f$
       *
       * \param[in] x The matrix to be transposed.
       */
      void transpose(const NullMatrix<Mem_, DT_, IT_, BlockWidth_, BlockHeight_> & x)
      {
        x = this->transpose();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");

        r.format();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");

        r.format();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        r.format();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r, const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

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
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows<Perspective::pod>(), "Vector size of y does not match!");

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
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockHeight_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

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
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows<Perspective::pod>(), "Vector size of y does not match!");

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
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockHeight_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

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
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ DOXY(alpha) = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows<Perspective::pod>(), "Vector size of y does not match!");

        r.copy(y);
      }

      ///@}

      /// \copydoc lump_rows()
      void lump_rows(VectorTypeL& lump) const
      {
        XASSERTM(lump.size() == rows(), "lump vector size does not match matrix row count!");
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
        XASSERTM(diag.size() == rows(), "diag size does not match matrix row count!");
        XASSERTM(rows() == columns(), "matrix is not square!");
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
        return VectorTypeL(this->rows());
      }

      // Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->columns());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index) const
      {
        return Index(0);
      }

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index, DT_ * const, IT_ * const, const Index, const Index = 1) const
      {
        // nothing to do here
      }

      void set_line_reverse(const Index, DT_ * const, const Index = 1)
      {
        // nothing to do here
      }
      /// \endcond

      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:
      /** \copydoc Adjactor::get_num_nodes_domain() */
      inline Index get_num_nodes_domain() const
      {
        return rows();
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      inline Index get_num_nodes_image() const
      {
        return columns();
      }

      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        XASSERTM(domain_node < rows(), "Domain node index out of range");
        return nullptr;
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        XASSERTM(domain_node < rows(), "Domain node index out of range");
        return nullptr;
      }
    }; // class NullMatrix
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_NULL_MATRIX_HPP

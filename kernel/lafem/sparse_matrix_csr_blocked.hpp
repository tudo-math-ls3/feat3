#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_CSR_BLOCKED_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_CSR_BLOCKED_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <fstream>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief CSR based blocked sparse matrix.
     *
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse row format.
     * Every non zero element is a small dense matrix block on itself, represented as Tiny::Matrix objects.
     * To be consistent with the external interface, the layout information are stored in block-scoped coordinates.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: column index per non zero element \n
     * _indices[1]: row start index (including matrix end index)\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) (multiple of blocksize)\n
     * _scalar_dt[0]: zero element
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
    class SparseMatrixCSRBlocked : public Container<Mem_, DT_, IT_>
    {
    private:
      Index & _size()
      {
        return this->_scalar_index.at(0);
      }

      Index & _rows()
      {
        return this->_scalar_index.at(1);
      }

      Index & _columns()
      {
        return this->_scalar_index.at(2);
      }

      Index & _used_elements()
      {
        return this->_scalar_index.at(3);
      }

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Our block height
      static constexpr int BlockHeight = BlockHeight_;
      /// Our block width
      static constexpr int BlockWidth = BlockWidth_;
      /// Our used layout type
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_csr;
      /// Value type, meaning the type of each block
      typedef Tiny::Matrix<DataType, BlockHeight, BlockWidth> ValueType;
      /// Compatible L-vector type
      typedef DenseVectorBlocked<Mem_, DT_, IT_, BlockHeight> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth> VectorTypeR;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixCSRBlocked() :
        Container<Mem_, DT_, IT_> (0)
      {
        CONTEXT("When creating SparseMatrixCSRBlocked");
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixCSRBlocked(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
        Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
      {
        CONTEXT("When creating SparseMatrixCSRBlocked");
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          Util::MemoryPool<Mem_>::instance()->increase_memory(i);

        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(raw_used_elements()));
        this->_elements_size.push_back(raw_used_elements());
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The.graph to create the matrix from
       *
       * Creates a CSR blocked matrix based on a given adjacency graph representing the sparsity pattern.
       */
      explicit SparseMatrixCSRBlocked(const Adjacency::Graph & graph) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixCSR");

        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nnze = graph.get_num_indices();

        // Create temporary vectors. Row and column pointer are block wise
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vrow_ptr(num_rows+1);
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vcol_idx(num_nnze);
        // The data array has to account for the block size
        LAFEM::DenseVector<Mem::Main, DT_, IT_> vdata(num_nnze*Index(BlockHeight_*BlockWidth_), DT_(0));

        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());
        IT_ * prow_ptr(vrow_ptr.elements());
        IT_ * pcol_idx(vcol_idx.elements());

        // build row-end
        prow_ptr[0] = IT_(dom_ptr[0]);
        for(Index i(0); i < num_rows; ++i)
          prow_ptr[i+1] = IT_(dom_ptr[i+1]);

        // build col-idx
        for(Index i(0); i < num_nnze; ++i)
          pcol_idx[i] = IT_(img_idx[i]);

        // build the matrix
        this->assign(SparseMatrixCSRBlocked<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_>(num_rows, num_cols, vcol_idx, vdata, vrow_ptr));
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] col_ind_in Vector with column indices.
       * \param[in] val_in Vector with non zero elements.
       * \param[in] row_ptr_in Vector with start indices of all rows into the val/col_ind arrays.
       * Note that this vector must also contain the end index of the last row and thus has a size of row_count + 1.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixCSRBlocked(const Index rows_in, const Index columns_in,
                                      DenseVector<Mem_, IT_, IT_> & col_ind_in, DenseVector<Mem_, DT_, IT_> & val_in, DenseVector<Mem_, IT_, IT_> & row_ptr_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        CONTEXT("When creating SparseMatrixCSRBlocked");
        ASSERT(val_in.size() % (BlockHeight_ * BlockWidth_) == 0, "Error: " + stringify(val_in.size()) + " not multiple of container blocksize!");
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(val_in.size() / Index(BlockHeight_ * BlockWidth_));
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(col_ind_in.elements());
        this->_indices_size.push_back(col_ind_in.size());
        this->_indices.push_back(row_ptr_in.elements());
        this->_indices_size.push_back(row_ptr_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixCSRBlocked(SparseMatrixCSRBlocked && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixCSRBlocked>(other))
      {
        CONTEXT("When moving SparseMatrixCSRBlocked");
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixCSRBlocked & operator= (SparseMatrixCSRBlocked && other)
      {
        CONTEXT("When moving SparseMatrixCSRBlocked");

        this->move(std::forward<SparseMatrixCSRBlocked>(other));

        return *this;
      }

      InsertWeakClone( SparseMatrixCSRBlocked );

      /** \brief Shallow copy operation
       *
       * Create a shallow copy of itself.
       *
       */
      SparseMatrixCSRBlocked shared() const
      {
        CONTEXT("When sharing SparseMatrixCSRBlocked");
        SparseMatrixCSRBlocked r;
        r.assign(*this);
        return r;
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixCSRBlocked<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_> & other)
      {
        CONTEXT("When converting SparseMatrixCSRBlocked");
        this->assign(other);
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixCSRBlocked & operator= (const SparseLayout<Mem_, IT_, layout_id> & layout_in)
      {
        CONTEXT("When assigning SparseMatrixCSRBlocked");

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();
        this->_scalar_dt.clear();

        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          Util::MemoryPool<Mem_>::instance()->increase_memory(i);

        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(raw_used_elements()));
        this->_elements_size.push_back(raw_used_elements());

        return *this;
      }

      /**
       * \brief Retrieve specific matrix element.
       *
       * \param[in] row The row of the matrix element.
       * \param[in] col The column of the matrix element.
       *
       * \returns Specific matrix element.
       */
      Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> operator()(Index row, Index col) const
      {
        CONTEXT("When retrieving SparseMatrixCSRBlocked element");

        ASSERT(row < rows(), "Error: " + stringify(row) + " exceeds sparse matrix csr row size " + stringify(rows()) + " !");
        ASSERT(col < columns(), "Error: " + stringify(col) + " exceeds sparse matrix csr column size " + stringify(columns()) + " !");

        for (Index i(Util::MemoryPool<Mem_>::get_element(this->_indices.at(1), row)) ; i < Util::MemoryPool<Mem_>::get_element(this->_indices.at(1), row + 1) ; ++i)
        {
          if (Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i) == col)
          {
            Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> t;
            Util::MemoryPool<Mem_>::download((DT_*)t.v, this->_elements.at(0) + i * Index(BlockHeight_*BlockWidth_), Index(BlockHeight_*BlockWidth_));
            return t;
          }
          if (Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i) > col)
            break; //return zero element
        }

        Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> ze((DT_(zero_element())));
        return ze;
      }

      /**
       * \brief Retrieve convenient sparse matrix layout object.
       *
       * \return An object containing the sparse matrix layout.
       */
      SparseLayout<Mem_, IT_, layout_id> layout() const
      {
        return SparseLayout<Mem_, IT_, layout_id>(this->_indices, this->_indices_size, this->_scalar_index);
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      const Index & rows() const
      {
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      const Index & columns() const
      {
        return this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve raw row count.
       *
       * \returns Raw matrix row count.
       */
      Index raw_rows() const
      {
        return this->_scalar_index.at(1) * Index(BlockHeight_);
      }

      /**
       * \brief Retrieve raw column count.
       *
       * \returns Raw matrix column count.
       */
      Index raw_columns() const
      {
        return this->_scalar_index.at(2) * BlockWidth_;
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      Index used_elements() const override
      {
        return this->_scalar_index.at(3);
      }

      /// The raw number of non zero elements of type DT_
      Index raw_used_elements() const
      {
        return used_elements() * Index(BlockHeight_ * BlockWidth_);
      }

      /**
       * \brief Retrieve column indices array.
       *
       * \returns Column indices array.
       */
      IT_ * col_ind()
      {
        return this->_indices.at(0);
      }

      /// \copydoc col_ind()
      /// const version.
      IT_ const * col_ind() const
      {
        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve non zero element array.
       *
       * \returns Non zero element array.
       */
      Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> * val()
      {
        return (Tiny::Matrix<DT_, BlockHeight_, BlockWidth_>*)this->_elements.at(0);
      }

      /// \copydoc val()
      /// const version.
      Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> const * val() const
      {
        return (Tiny::Matrix<DT_, BlockHeight_, BlockWidth_>*)this->_elements.at(0);
      }

      /**
       * \brief Get a pointer to the raw non zero element array.
       *
       * \returns Pointer to the raw non zero element array.
       */
      DT_ * raw_val()
      {
        return this->_elements.at(0);
      }

      /// \copydoc raw_val()
      /// const version.
      DT_ const * raw_val() const
      {
        return this->_elements.at(0);
      }

      /**
       * \brief Retrieve row start index array.
       *
       * \returns Row start index array.
       */
      IT_ * row_ptr()
      {
        return this->_indices.at(1);
      }

      /// \copydoc row_ptr()
      /// const version.
      IT_ const * row_ptr() const
      {
        return this->_indices.at(1);
      }

      /**
       * \brief Retrieve non zero element.
       *
       * \returns Non zero element.
       */
      const DT_ zero_element() const
      {
        return this->_scalar_dt.at(0);
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseMatrixCSRBlocked";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      void copy(const SparseMatrixCSRBlocked & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      template <typename Mem2_>
      void copy(const SparseMatrixCSRBlocked<Mem2_, DT_, IT_, BlockHeight_, BlockWidth_> & x)
      {
        this->_copy_content(x);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow y + \alpha x\f$
       *
       * \param[in] x The first summand matrix to be scaled.
       * \param[in] y The second summand matrix
       * \param[in] alpha A scalar to multiply x with.
       */
      void axpy(
                const SparseMatrixCSRBlocked & x,
                const SparseMatrixCSRBlocked & y,
                const DT_ alpha = DT_(1))
      {
        if (x.rows() != y.rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.used_elements() != y.used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
          Arch::Sum<Mem_>::value(this->raw_val(), x.raw_val(), y.raw_val(), this->raw_used_elements());
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          Arch::Difference<Mem_>::value(this->raw_val(), y.raw_val(), x.raw_val(), this->raw_used_elements());
        // r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
          Arch::Axpy<Mem_>::dv(this->raw_val(), alpha, x.raw_val(), y.raw_val(), this->raw_used_elements());
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const SparseMatrixCSRBlocked & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");

        Arch::Scale<Mem_>::value(this->val(), x.val(), alpha, this->raw_used_elements());
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        return Arch::Norm2<Mem_>::value(this->raw_val(), this->raw_used_elements());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        if (r.size() != this->raw_rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->raw_columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.elements(), this->raw_val(), this->col_ind(), this->row_ptr(),
                                                                                             x.elements(), this->rows(), columns(), used_elements());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->raw_columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.raw_elements(), this->raw_val(), this->col_ind(), this->row_ptr(),
                                                                                             x.elements(), this->rows(), columns(), used_elements());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x) const
      {
        if (r.size() != this->raw_rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.elements(), this->raw_val(), this->col_ind(), this->row_ptr(),
                                                                                             x.raw_elements(), this->rows(), columns(), used_elements());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r, const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.raw_elements(), this->raw_val(), this->col_ind(), this->row_ptr(),
                                                                                             x.raw_elements(), this->rows(), columns(), used_elements());
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->raw_rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->raw_columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->raw_rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.elements(), y.elements(), this->raw_val(), this->col_ind(),
                                                                                        this->row_ptr(), x.elements(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.elements(), alpha, x.elements(), y.elements(),
                                                                                      this->raw_val(), this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockHeight_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->raw_columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.raw_elements(), y.raw_elements(), this->raw_val(), this->col_ind(),
                                                                                        this->row_ptr(), x.elements(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          //r.copy(y);
          r.convert(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.raw_elements(), alpha, x.elements(), y.raw_elements(),
                                                                                      this->raw_val(), this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->raw_rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->raw_rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.elements(), y.elements(), this->raw_val(), this->col_ind(),
                                                                                        this->row_ptr(), x.raw_elements(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          //r.copy(y);
          r.convert(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.elements(), alpha, x.raw_elements(), y.elements(),
                                                                                      this->raw_val(), this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockHeight_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.raw_elements(), y.raw_elements(), this->raw_val(), this->col_ind(),
                                                                                        this->row_ptr(), x.raw_elements(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(r.raw_elements(), alpha, x.raw_elements(), y.raw_elements(),
                                                                                      this->raw_val(), this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }
      }
      ///@}

      /**
       * \brief SparseMatrixCSRBlocked comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_>
      friend bool operator== (const SparseMatrixCSRBlocked & a, const SparseMatrixCSRBlocked<Mem2_, DT_, IT_, BlockHeight_, BlockWidth_> & b)
      {
        CONTEXT("When comparing SparseMatrixCSRBlockeds");

        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;
        if (a.zero_element() != b.zero_element())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && a.get_indices().size() == 0 && b.get_elements().size() == 0 && b.get_indices().size() == 0)
          return true;

        IT_ * col_ind_a;
        IT_ * col_ind_b;
        DT_ * val_a;
        DT_ * val_b;
        IT_ * row_ptr_a;
        IT_ * row_ptr_b;

        bool ret(true);

        if(std::is_same<Mem::Main, Mem_>::value)
        {
          col_ind_a = (IT_*)a.col_ind();
          val_a = (DT_*)a.val();
          row_ptr_a = (IT_*)a.row_ptr();
        }
        else
        {
          col_ind_a = new IT_[a.used_elements()];
          Util::MemoryPool<Mem_>::instance()->template download<IT_>(col_ind_a, a.col_ind(), a.used_elements());
          val_a = new DT_[a.raw_used_elements()];
          Util::MemoryPool<Mem_>::instance()->template download<DT_>(val_a, a.raw_val(), a.raw_used_elements());
          row_ptr_a = new IT_[a.rows() + 1];
          Util::MemoryPool<Mem_>::instance()->template download<IT_>(row_ptr_a, a.row_ptr(), a.rows() + 1);
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
        {
          col_ind_b = (IT_*)b.col_ind();
          val_b = (DT_*)b.val();
          row_ptr_b = (IT_*)b.row_ptr();
        }
        else
        {
          col_ind_b = new IT_[b.used_elements()];
          Util::MemoryPool<Mem2_>::instance()->template download<IT_>(col_ind_b, b.col_ind(), b.used_elements());
          val_b = new DT_[b.raw_used_elements()];
          Util::MemoryPool<Mem2_>::instance()->template download<DT_>(val_b, b.raw_val(), b.raw_used_elements());
          row_ptr_b = new IT_[b.rows() + 1];
          Util::MemoryPool<Mem2_>::instance()->template download<IT_>(row_ptr_b, b.row_ptr(), b.rows() + 1);
        }

        for (Index i(0) ; i < a.used_elements() ; ++i)
        {
          if (col_ind_a[i] != col_ind_b[i])
          {
            ret = false;
            break;
          }
        }
        if (ret)
        {
          for (Index i(0) ; i < a.raw_used_elements() ; ++i)
          {
            if (val_a[i] != val_b[i])
            {
              ret = false;
              break;
            }
          }
        }
        if (ret)
        {
          for (Index i(0) ; i < a.rows() + 1; ++i)
          {
            if (row_ptr_a[i] != row_ptr_b[i])
            {
              ret = false;
              break;
            }
          }
        }

        if(! std::is_same<Mem::Main, Mem_>::value)
        {
          delete[] col_ind_a;
          delete[] val_a;
          delete[] row_ptr_a;
        }
        if(! std::is_same<Mem::Main, Mem2_>::value)
        {
          delete[] col_ind_b;
          delete[] val_b;
          delete[] row_ptr_b;
        }

        return ret;
      }

      /**
       * \brief SparseMatrixCSRBlocked streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCSRBlocked & b)
      {
        lhs << "[" << std::endl;
        for (Index i(0) ; i < b.rows() ; ++i)
        {
          for (int k(0) ; k < BlockHeight_ ; ++k)
          {
            lhs << "[";
            for (Index j(0) ; j < b.columns() ; ++j)
            {
              for (int l(0) ; l < BlockWidth_ ; ++l)
                lhs << "  " << b(i, j).v[k][l];
            }
            lhs << "]" << std::endl;
          }
        }
        lhs << "]" << std::endl;

        return lhs;
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
      /// \endcond
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_MATRIX_CSR_BLOCKED_HPP

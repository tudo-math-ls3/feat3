#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_BAND_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_BAND_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/adjacency/graph.hpp>

#include <fstream>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * Supported File modes.
     */
    enum class FiniteElementType
    {
      fe_q2 = 0,
      fe_q1
    };

    template <typename Mem_, typename DT_, FiniteElementType BT_, typename IT_>
    class SparseMatrixBand : public Container<Mem_, DT_, IT_>, public MatrixBase
    {
    };

    /**
     * \brief sparse band matrix for Q1-FE.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its diagonal entries matching to a Q1-discretization on lexicographical numbering.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: number of bands \n
     * _scalar_index[5]: number of nodes per row \n
     * _scalar_index[6]: number of nodes per column \n
     * _scalar_dt[0]: zero element
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_, typename IT_>
    class SparseMatrixBand<Mem_, DT_, FiniteElementType::fe_q1, IT_> : public Container<Mem_, DT_, IT_>, public MatrixBase
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

      Index & _num_of_bands()
      {
        return this->_scalar_index.at(4);
      }

      Index & _nodes_per_row()
      {
        return this->_scalar_index.at(5);
      }

      Index & _nodes_per_column()
      {
        return this->_scalar_index.at(6);
      }

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Compatible L-vector type
      typedef DenseVector<MemType, DataType> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<MemType, DataType> VectorTypeR;
      /// ImageIterator typedef for Adjactor interface implementation
      typedef const IT_* ImageIterator;


      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixBand() :
        Container<Mem_, DT_, IT_> (0)
      {
        CONTEXT("When creating SparseMatrixBand");
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] nodes_per_row number of nodes per row \n
       * \param[in] nodes_per_column number of nodes per column \n
       * \param[in] val_in Vector with non zero elements.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBand(const Index nodes_per_row_in, const Index nodes_per_column_in,
                                DenseVector<Mem_, DT_, IT_> & val_in) :
        Container<Mem_, DT_, IT_>(nodes_per_row_in * nodes_per_row_in * nodes_per_column_in * nodes_per_column_in)
      {
        if (val_in.size() != 9 * nodes_per_row_in * nodes_per_column_in)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Size of vector does not match to number of nodes!");
        }

        CONTEXT("When creating SparseMatrixBand");
        this->_scalar_index.push_back(nodes_per_row_in * nodes_per_column_in);
        this->_scalar_index.push_back(nodes_per_row_in * nodes_per_column_in);
        this->_scalar_index.push_back(9 * nodes_per_row_in * nodes_per_column_in - 6 * nodes_per_row_in - 2);
        this->_scalar_index.push_back(9);
        this->_scalar_index.push_back(nodes_per_row_in);
        this->_scalar_index.push_back(nodes_per_column_in);
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in number of nodes per row/column.
       * \param[in] val_in Vector with non zero elements.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBand(const Index size_in,
                                DenseVector<Mem_, DT_, IT_> & val_in) :
        Container<Mem_, DT_, IT_>(size_in * size_in * size_in * size_in)
      {
        const Index size_in_sqr(Math::sqr(size_in));

        if (val_in.size() != 9 * size_in_sqr)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Size of vector does not match to number of nodes!");
        }
        CONTEXT("When creating SparseMatrixBand");
        this->_scalar_index.push_back(size_in_sqr);
        this->_scalar_index.push_back(size_in_sqr);
        this->_scalar_index.push_back(9 * size_in_sqr - 6 * size_in - 2);
        this->_scalar_index.push_back(9);
        this->_scalar_index.push_back(size_in);
        this->_scalar_index.push_back(size_in);
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixBand(SparseMatrixBand && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixBand>(other))
      {
        CONTEXT("When moving SparseMatrixBand");
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixBand & operator= (SparseMatrixBand && other)
      {
        CONTEXT("When moving SparseMatrixBand");

        this->move(std::forward<SparseMatrixBand>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a deep copy of itself.
       *
       */
      SparseMatrixBand clone()
      {
        SparseMatrixBand t;
        t.clone(*this);
        return t;
      }

      using Container<Mem_, DT_, IT_>::clone;

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixBand<Mem2_, DT2_, FiniteElementType::fe_q1, IT2_> & other)
      {
        CONTEXT("When converting SparseMatrixBand");
        this->assign(other);
      }

      /**
       * \brief Retrieve specific matrix element.
       *
       * \param[in] row The row of the matrix element.
       * \param[in] col The column of the matrix element.
       *
       * \returns Specific matrix element.
       */
      DT_ operator()(Index row, Index col) const
      {
        CONTEXT("When retrieving SparseMatrixBand element");

        ASSERT(row < rows(), "Error: " + stringify(row) + " exceeds sparse matrix band row size " + stringify(rows()) + " !");
        ASSERT(col < columns(), "Error: " + stringify(col) + " exceeds sparse matrix band column size " + stringify(columns()) + " !");

        const Index tnpr(this->_scalar_index.at(5));
        const Index trows(this->_scalar_index.at(1));

        if (row == col + tnpr + 1)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), row);
        }
        else if (row == col + tnpr)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), trows + row);
        }
        else if (row == col + tnpr - 1)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), 2 * trows + row);
        }
        else if (row == col + 1)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), 3 * trows + row);
        }
        else if (row == col)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), 4 * trows + row);
        }
        else if (row + 1 == col)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), 5 * trows + row);
        }
        else if (row + tnpr - 1 == col)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), 6 * trows + row);
        }
        else if (row + tnpr == col)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), 7 * trows + row);
        }
        else if (row + tnpr + 1 == col)
        {
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), 8 * trows + row);
        }

        return this->_scalar_dt.at(0);
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
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      const Index & used_elements() const override
      {
        return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve number of bands.
       *
       * \returns Number of bands.
       */
      const Index & num_of_bands() const
      {
        return this->_scalar_index.at(4);
      }

      /**
       * \brief Retrieve number of nodes per row.
       *
       * \returns Number of nodes per row.
       */
      const Index & nodes_per_row() const
      {
        return this->_scalar_index.at(5);
      }

      /**
       * \brief Retrieve number of nodes per column.
       *
       * \returns Number of nodes per column.
       */
      const Index & nodes_per_column() const
      {
        return this->_scalar_index.at(6);
      }

      /**
       * \brief Retrieve element array.
       *
       * \returns Non zero element array.
       */
      DT_ * val()
      {
        return this->_elements.at(0);
      }

      DT_ const * val() const
      {
        return this->_elements.at(0);
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
        return "SparseMatrixBand";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      void copy(const SparseMatrixBand & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      template <typename Mem2_>
      void copy(const SparseMatrixBand<Mem2_, DT_, FiniteElementType::fe_q1, IT_> & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Calculate \f$this \leftarrow y + \alpha x\f$
       *
       * \param[in] x The first summand matrix to be scaled.
       * \param[in] y The second summand matrix
       * \param[in] alpha A scalar to multiply x with.
       */
      template <typename Algo_>
      void axpy(
                const SparseMatrixBand & x,
                const SparseMatrixBand & y,
                const DT_ alpha = DT_(1))
      {
        if (x.nodes_per_row() != y.nodes_per_row())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix number of nodes per row do not match!");
        if (x.nodes_per_row() != this->nodes_per_row())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix number of nodes per row do not match!");
        if (x.nodes_per_column() != y.nodes_per_column())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix number of nodes per column do not match!");
        if (x.nodes_per_column() != this->nodes_per_column())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix number of nodes per column do not match!");

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
          Arch::Sum<Mem_, Algo_>::value(this->val(), x.val(), y.val(), this->used_elements());
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          Arch::Difference<Mem_, Algo_>::value(this->val(), y.val(), x.val(), this->used_elements());
        // r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
          Arch::Axpy<Mem_, Algo_>::dv(this->val(), alpha, x.val(), y.val(), _rows() * _num_of_bands());
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      template <typename Algo_>
      void scale(const SparseMatrixBand & x, const DT_ alpha)
      {
        if (x.nodes_per_row() != this->nodes_per_row())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix number of nodes per row do not match!");
        if (x.nodes_per_column() != this->nodes_per_column())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix number of nodes per column do not match!");

        Arch::Scale<Mem_, Algo_>::value(this->val(), x.val(), alpha, _rows() * _num_of_bands());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      template<typename Algo_>
      void apply(DenseVector<Mem_,DT_, IT_>& r, const DenseVector<Mem_, DT_, IT_>& x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        Arch::ProductMatVec<Mem_, Algo_>::band_q2_d2(r.elements(),
                                                     this->val(),
                                                     x.elements(),
                                                     this->_scalar_index.at(5),
                                                     this->_scalar_index.at(6));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      template<typename Algo_>
      void apply(DenseVector<Mem_,DT_, IT_>& r,
                 const DenseVector<Mem_, DT_, IT_>& x,
                 const DenseVector<Mem_, DT_, IT_>& y,
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
          Arch::Defect<Mem_, Algo_>::band_q2_d2(r.elements(),
                                                y.elements(),
                                                this->val(),
                                                x.elements(),
                                                this->_scalar_index.at(5),
                                                this->_scalar_index.at(6));
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_, Algo_>::band_q2_d2(r.elements(),
                                              alpha,
                                              x.elements(),
                                              y.elements(),
                                              this->val(),
                                              this->_scalar_index.at(5),
                                              this->_scalar_index.at(6));
        }
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(this->rows());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->columns());
      }
    };

    /**
     * \brief SparseMatrixBand streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_, typename IT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixBand<Mem_, DT_, FiniteElementType::fe_q1, IT_> & b)
    {
      lhs << "[" << std::endl;
      for (Index i(0) ; i < b.rows() ; ++i)
      {
        lhs << "[";
        for (Index j(0) ; j < b.columns() ; ++j)
        {
          lhs << "  " << b(i, j);
        }
        lhs << "]" << std::endl;
      }
      lhs << "]" << std::endl;

      return lhs;
    }

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_MATRIX_BAND_HPP

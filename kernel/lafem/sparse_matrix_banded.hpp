#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP 1

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
     * \brief sparse banded matrix
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its diagonal entries
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: vector of offsets (main diagonal has offset rows-1)\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: number of bands \n
     * _scalar_dt[0]: zero element
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_, typename IT_>
    class SparseMatrixBanded : public Container<Mem_, DT_, IT_>, public MatrixBase
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

      Index & _num_of_offsets()
      {
        return this->_scalar_index.at(4);
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
      explicit SparseMatrixBanded() :
        Container<Mem_, DT_, IT_> (0)
      {
        CONTEXT("When creating SparseMatrixBanded");
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] val_in The vector with non zero elements.
       * \param[in] offsets_in The vector of offsets.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBanded(const Index rows_in, const Index columns_in,
                                  DenseVector<Mem_, DT_, IT_> & val_in,
                                  DenseVector<Mem_, IT_> & offsets_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        if (val_in.size() != rows_in * offsets_in.size())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Size of values does not match to number of offsets and row count!");
        }

        const Index offset0(offsets_in(0));
        const Index offset1(offsets_in(offsets_in.size() - 1));

        const Index row0(             rows_in - offset0 - 1);
        const Index row1(columns_in + rows_in - offset1 - 1);

        if (offset0 < rows_in - 1 && offset1 >= rows_in && row0 >= row1)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Offsets can't be used for this matrix-dimensions!");
        }

        CONTEXT("When creating SparseMatrixBanded");
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(0); // TODO: Muss noch berechnet werden
        this->_scalar_index.push_back(offsets_in.size());
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(offsets_in.elements());
        this->_indices_size.push_back(offsets_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixBanded(SparseMatrixBanded && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixBanded>(other))
      {
        CONTEXT("When moving SparseMatrixBanded");
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixBanded & operator= (SparseMatrixBanded && other)
      {
        CONTEXT("When moving SparseMatrixBanded");

        this->move(std::forward<SparseMatrixBanded>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a deep copy of itself.
       *
       */
      SparseMatrixBanded clone()
      {
        SparseMatrixBanded t;
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
      void convert(const SparseMatrixBanded<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting SparseMatrixBanded");
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
        CONTEXT("When retrieving SparseMatrixBanded element");

        ASSERT(row < rows(), "Error: " + stringify(row) + " exceeds sparse matrix banded row size " + stringify(rows()) + " !");
        ASSERT(col < columns(), "Error: " + stringify(col) + " exceeds sparse matrix banded column size " + stringify(columns()) + " !");

        const Index trows(this->_scalar_index.at(1));

        for (Index i(0); i < this->_scalar_index.at(4); ++i)
        {
          const Index toffset(MemoryPool<Mem_>::get_element(this->_indices.at(0), i));
          if (row + toffset + 1 == col + trows)
          {
            return MemoryPool<Mem_>::get_element(this->_elements.at(0), i * trows + row);
          }
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
       * \brief Retrieve number of offsets.
       *
       * \returns Number of offsets.
       */
      const Index & num_of_offsets() const
      {
        return this->_scalar_index.at(4);
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
       * \brief Retrieve offsets array.
       *
       * \returns offsets array.
       */
      DT_ * offsets()
      {
        return this->_indices.at(0);
      }

      DT_ const * offsets() const
      {
        return this->_indices.at(0);
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
        return "SparseMatrixBanded";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      void copy(const SparseMatrixBanded & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      template <typename Mem2_>
      void copy(const SparseMatrixBanded<Mem2_, DT_, IT_> & x)
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
                const SparseMatrixBanded & x,
                const SparseMatrixBanded & y,
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
        if (x.num_of_offsets() != y.num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.num_of_offsets() != this->num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.used_elements() != y.used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
          Arch::Sum<Mem_, Algo_>::value(this->val(), x.val(), y.val(), this->rows() * this->num_of_offsets());
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          Arch::Difference<Mem_, Algo_>::value(this->val(), y.val(), x.val(), this->rows() * this->num_of_offsets());
        // r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
          Arch::Axpy<Mem_, Algo_>::dv(this->val(), alpha, x.val(), y.val(), this->rows() * this->num_of_offsets());
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      template <typename Algo_>
      void scale(const SparseMatrixBanded & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.num_of_offsets() != this->num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        Arch::Scale<Mem_, Algo_>::value(this->val(), x.val(), alpha, this->rows() * this->num_of_offsets());
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

        Arch::ProductMatVec<Mem_, Algo_>::banded(r.elements(),
                                                 this->val(),
                                                 x.elements(),
                                                 this->_scalar_index.at(5),
                                                 this->_scalar_index.at(6)); // TODO
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
          Arch::Defect<Mem_, Algo_>::banded(r.elements(),
                                            y.elements(),
                                            this->val(),
                                            x.elements(),
                                            this->_scalar_index.at(5),
                                            this->_scalar_index.at(6)); // TODO
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_, Algo_>::banded(r.elements(),
                                          alpha,
                                          x.elements(),
                                          y.elements(),
                                          this->val(),
                                          this->_scalar_index.at(5),
                                          this->_scalar_index.at(6)); // TODO
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
     * \brief SparseMatrixBanded streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_, typename IT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixBanded<Mem_, DT_, IT_> & b)
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

#endif // KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP

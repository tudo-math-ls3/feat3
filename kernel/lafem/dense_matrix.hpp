#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/matrix_base.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Dense data matrix class template.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a matrix of continuous data in memory. \n\n
     * Data survey: \n
     * _elements[0]: raw number values \n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class DenseMatrix : public Container<Mem_, DT_>, public MatrixBase
    {
      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Mem_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit DenseMatrix() :
          Container<Mem_, DT_> (0)
        {
          CONTEXT("When creating DenseMatrix");

          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         *
         * Creates a matrix with given dimensions.
         */
        explicit DenseMatrix(Index rows, Index columns) :
          Container<Mem_, DT_>(rows * columns)
        {
          CONTEXT("When creating DenseMatrix");

          this->_scalar_index.at(0) = rows * columns;
          this->_scalar_index.push_back(rows);
          this->_scalar_index.push_back(columns);

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(0) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         * \param[in] value The value, each element will be set to.
         *
         * Creates a matrix with given dimensions and value.
         */
        explicit DenseMatrix(Index rows, Index columns, DT_ value) :
          Container<Mem_, DT_>(rows * columns)
        {
          CONTEXT("When creating DenseMatrix");

          this->_scalar_index.at(0) = rows * columns;
          this->_scalar_index.push_back(rows);
          this->_scalar_index.push_back(columns);
          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(0) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(0));
          MemoryPool<Mem_>::instance()->set_memory(this->_elements.at(0), value, this->_scalar_index.at(0));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a shallow copy of a given matrix.
         */
        DenseMatrix(const DenseMatrix<Mem_, DT_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying DenseMatrix");
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a copy of a given matrix from another memory architecture.
         */
        template <typename Arch2_, typename DT2_>
        DenseMatrix(const DenseMatrix<Arch2_, DT2_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying DenseMatrix");
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        DenseMatrix<Mem_, DT_> clone()
        {
          CONTEXT("When cloning DenseMatrix");

          DenseMatrix<Mem_, DT_> t;
          ((Container<Mem_, DT_>&)t).clone(*this);
          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns another matrix to the target matrix.
         */
        DenseMatrix<Mem_, DT_> & operator= (const DenseMatrix<Mem_, DT_> & other)
        {
          CONTEXT("When assigning DenseMatrix");

          assign(other);

          return *this;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns a matrix from another memory architecture to the target matrix.
         */
        template <typename Arch2_, typename DT2_>
        DenseMatrix<Mem_, DT_> & operator= (const DenseMatrix<Arch2_, DT2_> & other)
        {
          CONTEXT("When assigning DenseMatrix");

          assign(other);

          return *this;
        }

        /**
         * \brief Get a pointer to the data array.
         *
         * \returns Pointer to the data array.
         */
        DT_ * elements()
        {
          return this->_elements.at(0);
        }

        const DT_ * elements() const
        {
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         *
         * \returns Specific matrix element.
         */
        const DT_ operator()(Index row, Index col) const
        {
          CONTEXT("When retrieving DenseMatrix element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->columns()) + " !");
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), row * this->columns() + col);
        }

        /**
         * \brief Set specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         * \param[in] value The value to be set.
         */
        void operator()(Index row, Index col, DT_ value)
        {
          CONTEXT("When setting DenseMatrix element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->columns()) + " !");
          MemoryPool<Mem_>::set_memory(this->_elements.at(0) + row * this->columns() + col, value);
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
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String type_name()
        {
          return "DenseMatrix";
        }
    };

    /**
     * \brief DenseMatrix comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Mem_, typename Arch2_, typename DT_> bool operator== (const DenseMatrix<Mem_, DT_> & a, const DenseMatrix<Arch2_, DT_> & b)
    {
      CONTEXT("When comparing DenseMatrices");

      if (a.size() != b.size())
        return false;
      if (a.get_elements().size() != b.get_elements().size())
        return false;
      if (a.get_indices().size() != b.get_indices().size())
        return false;
      if (a.rows() != b.rows())
        return false;
      if (a.columns() != b.columns())
        return false;

      for (Index i(0) ; i < a.rows() ; ++i)
        for (Index j(0) ; j < a.columns() ; ++j)
        if (a(i, j) != b(i, j))
          return false;

      return true;
    }

    /**
     * \brief DenseMatrix streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_>
      std::ostream &
      operator<< (std::ostream & lhs, const DenseMatrix<Mem_, DT_> & b)
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

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP

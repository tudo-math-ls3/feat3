#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Dense data matrix class template.
     *
     * \tparam Arch_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a matrix of continuous data in memory.
     *
     * \author Dirk Ribbrock
     */
    template <typename Arch_, typename DT_>
    class DenseMatrix : public Container<Arch_, DT_>
    {
      private:
        /// Pointer to our elements.
        DT_ * _pelements;
        /// Our row count.
        Index _rows;
        /// Our column count.
        Index _columns;

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Arch_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit DenseMatrix() :
          Container<Arch_, DT_> (0),
          _rows(0),
          _columns(0)
        {
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
          Container<Arch_, DT_>(rows * columns)
        {
          this->_size = rows * columns;
          this->_rows = rows;
          this->_columns = columns;

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(this->_size * sizeof(DT_)));
          this->_elements_size.push_back(this->_size);
          _pelements = this->_elements.at(0);
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
          Container<Arch_, DT_>(rows * columns)
        {
          this->_size = rows * columns;
          this->_rows = rows;
          this->_columns = columns;
          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(this->_size * sizeof(DT_)));
          this->_elements_size.push_back(this->_size);
          _pelements = this->_elements.at(0);

          MemoryPool<Arch_>::instance()->set_memory(_pelements, value, this->_size);
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a shallow copy of a given matrix.
         */
        DenseMatrix(const DenseMatrix<Arch_, DT_> & other) :
          Container<Arch_, DT_>(other),
          _rows(other._rows),
          _columns(other._columns)
        {
          _pelements = this->_elements.at(0);
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
          Container<Arch_, DT_>(other),
          _rows(other.rows()),
          _columns(other.columns())
        {
          _pelements = this->_elements.at(0);
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        DenseMatrix<Arch_, DT_> clone()
        {
          DenseMatrix<Arch_, DT_> t(this->_rows, this->_columns);

          void * pdest(t.elements());
          const void * psrc(this->elements());
          MemoryPool<Arch_>::copy(pdest, psrc, this->_size * sizeof(DT_));

          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns another matrix to the target matrix.
         */
        DenseMatrix<Arch_, DT_> & operator= (const DenseMatrix<Arch_, DT_> & other)
        {
          if (this == &other)
            return *this;

          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();

          std::vector<DT_ *> new_elements = other.get_elements();
          std::vector<Index *> new_indices = other.get_indices();

          this->_elements.assign(new_elements.begin(), new_elements.end());
          this->_indices.assign(new_indices.begin(), new_indices.end());
          this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
          this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());

          _pelements = this->_elements.at(0);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

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
        DenseMatrix<Arch_, DT_> & operator= (const DenseMatrix<Arch2_, DT2_> & other)
        {
          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();


          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(other.size() * sizeof(DT_)));
          this->_elements_size.push_back(this->_size);
          this->_pelements = this->_elements.at(0);

          Index src_size(other.get_elements_size().at(0) * sizeof(DT2_));
          Index dest_size(other.get_elements_size().at(0) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_elements().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(0), temp, dest_size);
          ::free(temp);

          return *this;
        }

        /**
         * \brief Get a pointer to the data array.
         *
         * \returns Pointer to the data array.
         */
        DT_ * elements()
        {
          return _pelements;
        }

        const DT_ * elements() const
        {
          return _pelements;
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
          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->_columns) + " !");
          return MemoryPool<Arch_>::get_element(_pelements, row * this->_columns + col);
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
          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->_columns) + " !");
          MemoryPool<Arch_>::modify_element(_pelements, row * this->_columns + col, value);
        }

        /**
         * \brief Retrieve matrix row count.
         *
         * \returns Matrix row count.
         */
        const Index & rows() const
        {
          return this->_rows;
        }

        /**
         * \brief Retrieve matrix column count.
         *
         * \returns Matrix column count.
         */
        const Index & columns() const
        {
          return this->_columns;
        }
    };

    /**
     * \brief DenseMatrix comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Arch_, typename Arch2_, typename DT_> bool operator== (const DenseMatrix<Arch_, DT_> & a, const DenseMatrix<Arch2_, DT_> & b)
    {
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
    template <typename Arch_, typename DT_>
      std::ostream &
      operator<< (std::ostream & lhs, const DenseMatrix<Arch_, DT_> & b)
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

#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/util/graph.hpp>

#include <vector>
#include <map>
#include <algorithm>


namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename DT_>
    class SparseMatrixCOO : public Container<Arch_, DT_>
    {
    };

    /**
     * \brief Coordinate based sparse matrix.
     *
     * \tparam Arch_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements alongside their coordinates in a single heap.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_>
    class SparseMatrixCOO<Mem::Main, DT_> : public Container<Mem::Main, DT_>
    {
      private:
        /// Row count.
        Index _rows;
        /// Column count.
        Index _columns;
        /// One dimensional mapping of matrix coordinates to non zero elements.
        std::map<Index, DT_> _elements;
        /// Our non zero element.
        DT_ _zero_element;
        /// Our non zero element count.
        Index _used_elements;

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Mem::Main MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit SparseMatrixCOO() :
          Container<Mem::Main, DT_> (0),
          _rows(0),
          _columns(0),
          _zero_element(DT_(0)),
          _used_elements(0)
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
        explicit SparseMatrixCOO(Index rows, Index columns) :
          Container<Mem::Main, DT_>(rows * columns),
          _zero_element(DT_(0)),
          _used_elements(0)
        {
          this->_rows = rows;
          this->_columns = columns;
        }

        /**
         * \brief Constructor
         *
         * \param[in] graph The Graph, the matrix will be created from.
         *
         * Creates a matrix from a given graph.
         */
        explicit SparseMatrixCOO(const Graph & graph) :
          Container<Mem::Main, DT_>(graph.get_num_nodes_domain() * graph.get_num_nodes_image()),
          _zero_element(DT_(0)),
          _used_elements(graph.get_num_indices())
        {
          this->_rows = graph.get_num_nodes_domain();
          this->_columns = graph.get_num_nodes_image();

          for (Index i(0) ; i < this->_rows ; ++i)
          {
            typename Graph::ImageIterator it(graph.image_begin(i));
            typename Graph::ImageIterator jt(graph.image_end(i));
            for( ; it != jt ; ++it)
            {
              (*this)(i, *it, DT_(0));
            }
          }
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a shallow copy of a given matrix.
         */
        SparseMatrixCOO(const SparseMatrixCOO<Mem::Main, DT_> & other) :
          Container<Mem::Main, DT_>(other),
          _rows(other._rows),
          _columns(other._columns),
          _zero_element(other._zero_element),
          _used_elements(other.used_elements())
        {
            _elements.insert(other._elements.begin(), other._elements.end());
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a copy of a given matrix from another memory architecture.
         */
        template <typename Arch2_, typename DT2_>
        SparseMatrixCOO(const SparseMatrixCOO<Arch2_, DT2_> & other) :
          Container<Mem::Main, DT_>(other),
          _rows(other._rows),
          _columns(other._columns),
          _zero_element(other._zero_element),
          _used_elements(other.used_elements())
        {
            _elements.insert(other._elements.begin(), other._elements.end());
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixCOO<Mem::Main, DT_> clone()
        {
          SparseMatrixCOO<Mem::Main, DT_> t(this->_rows, this->_columns);
          t._elements.insert(this->_elements.begin(), this->_elements.end());
          return t;
        }

        /**
         * \brief Set specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         * \param[in] value The value to be set.
         */
        void operator()(Index row, Index col, DT_ val)
        {
          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds sparse matrix coo row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds sparse matrix coo column size " + stringify(this->_columns) + " !");

          typename std::map<Index, DT_>::const_iterator it(_elements.find(row * _columns + col));
          if (it == _elements.end())
            ++_used_elements;

          _elements[row * _columns + col] = val;
        }

        /**
         * \brief Set specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         * \param[in] value The value to be set.
         */
        void insert(Index row, Index col, DT_ val)
        {
          (*this)(row, col, val);
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
          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds sparse matrix coo row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds sparse matrix coo column size " + stringify(this->_columns) + " !");

          typename std::map<Index, DT_>::const_iterator it(_elements.find(row * _columns + col));
          if (it == _elements.end())
            return _zero_element;
          else
          {
            return it->second;
          }
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

        /**
         * \brief Retrieve non zero element count.
         *
         * \returns Non zero element count.
         */
        const Index & used_elements() const
        {
          return this->_used_elements;
        }

        /**
         * \brief Retrieve non zero element map.
         *
         * \returns Non zero element map.
         */
        const std::map<Index, DT_> & elements() const
        {
          return _elements;
        }

        /**
         * \brief Retrieve non zero element.
         *
         * \returns Non zero element.
         */
        const DT_ zero_element() const
        {
          return _zero_element;
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String name()
        {
          return "SparseMatrixCOO";
        }
    };

    /**
     * \brief SparseMatrixCOO comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Arch_, typename DT_> bool operator== (const SparseMatrixCOO<Arch_, DT_> & a, const SparseMatrixCOO<Arch_, DT_> & b)
    {
      if (a.rows() != b.rows())
        return false;
      if (a.columns() != b.columns())
        return false;

      for (Index i(0) ; i < a.rows() ; ++i)
      {
        for (Index j(0) ; j < a.columns() ; ++j)
        {
          if (a(i, j) != b(i, j))
            return false;
        }
      }

      return true;
    }

    /**
     * \brief SparseMatrixCOO streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Arch_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixCOO<Arch_, DT_> & b)
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

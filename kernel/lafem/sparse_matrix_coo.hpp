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
#include <fstream>
#include <iostream>


namespace FEAST
{
  namespace LAFEM
  {
    //forward declarations
    template <typename Arch_, typename DT_>
    class SparseMatrixELL;

    template <typename Arch_, typename DT_>
    class SparseMatrixCSR;

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
          CONTEXT("When creating SparseMatrixCOO");
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
          _rows(rows),
          _columns(columns),
          _zero_element(DT_(0)),
          _used_elements(0)
        {
          CONTEXT("When creating SparseMatrixCOO");
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source ell matrix.
         *
         * Creates a matrix from the given source matrix.
         */
        explicit SparseMatrixCOO(const SparseMatrixELL<Mem::Main, DT_> & other) :
          Container<Mem::Main, DT_>(other.rows() * other.columns()),
          _rows(other.rows()),
          _columns(other.columns()),
          _zero_element(DT_(0)),
          _used_elements(other.used_elements())
        {
          CONTEXT("When creating SparseMatrixCOO from SparseMatrixELL");

          for (Index i(0) ; i < other.num_cols_per_row() * other.stride() ; ++i)
          {
            if (other.Ax()[i] != DT_(0))
            {
              (*this)(i%other.stride(), other.Aj()[i], other.Ax()[i]);
            }
          }
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source csr matrix.
         *
         * Creates a matrix from the given source matrix.
         */
        explicit SparseMatrixCOO(const SparseMatrixCSR<Mem::Main, DT_> & other) :
          Container<Mem::Main, DT_>(other.rows() * other.columns()),
          _rows(other.rows()),
          _columns(other.columns()),
          _zero_element(DT_(0)),
          _used_elements(other.used_elements())
        {
          for (Index row(0) ; row < _rows ; ++row)
          {
            const Index end(other.row_ptr_end()[row]);
            for (Index i(other.row_ptr()[row]) ; i < end ; ++i)
            {
              if (other.val()[i] != DT_(0))
              {
                (*this)(row, other.col_ind()[i], other.val()[i]);
              }
            }
          }
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
          CONTEXT("When creating SparseMatrixCOO");

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
         * \brief Constructor
         *
         * \param[in] filename The source file to be read in.
         *
         * Creates a matrix based on the source file.
         */
        explicit SparseMatrixCOO(FileMode mode, String filename) :
          Container<Mem::Main, DT_>(0),
          _rows(0),
          _columns(0),
          _zero_element(DT_(0))
        {
          CONTEXT("When creating SparseMatrixELL");

          if (mode != fm_m)
                throw InternalError("Filemode not supported!");

          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);

          std::vector<Index> rowsv;
          std::vector<Index> colsv;
          std::vector<DT_> valsv;

          String line;
          std::getline(file, line);
          while(!file.eof())
          {
            std::getline(file, line);

            if(line.find("]", 0) < line.npos)
              break;

            if(line[line.size()-1] == ';')
              line.resize(line.size()-1);

            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String srow(line, 0, end);
            Index row(atol(srow.c_str()));
            --row;
            _rows = std::max(row+1, _rows);
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            Index col(atol(scol.c_str()));
            --col;
            _columns = std::max(col+1, _columns);
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String sval(line, 0, end);
            DT_ val(atof(sval.c_str()));

            rowsv.push_back(row);
            colsv.push_back(col);
            valsv.push_back(val);

          }
          this->_size = this->_rows * this->_columns;

          for (Index i(0) ; i < rowsv.size() ; ++i)
          {
            (*this)(rowsv.at(i), colsv.at(i), valsv.at(i));
          }



          file.close();
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
          CONTEXT("When copying SparseMatrixCOO");

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
          CONTEXT("When copying SparseMatrixCOO");

          _elements.insert(other._elements.begin(), other._elements.end());
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixCOO<Mem::Main, DT_> clone()
        {
          CONTEXT("When cloning SparseMatrixCOO");

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
          CONTEXT("When setting SparseMatrixCOO element");

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
        /*void insert(Index row, Index col, DT_ val)
        {
          (*this)(row, col, val);
        }*/

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
          CONTEXT("When retrieving SparseMatrixCOO element");

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
        static String type_name()
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
      CONTEXT("When comparing SparseMatrixCOOs");

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

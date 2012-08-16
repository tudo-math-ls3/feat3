#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/absolute.hpp>

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
      private:
        Index _rows;
        Index _columns;
        std::map<Index, DT_> _elements;
        DT_ _zero_element;
        Index _used_elements;

      public:
        typedef DT_ data_type;

        explicit SparseMatrixCOO() :
          Container<Arch_, DT_> (0),
          _rows(0),
          _columns(0),
          _zero_element(DT_(0)),
          _used_elements(0)
        {
        }

        explicit SparseMatrixCOO(Index rows, Index columns) :
          Container<Arch_, DT_>(rows * columns),
          _zero_element(DT_(0)),
          _used_elements(0)
        {
          this->_size = rows * columns;
          this->_rows = rows;
          this->_columns = columns;
        }

        SparseMatrixCOO(const SparseMatrixCOO<Arch_, DT_> & other) :
          Container<Arch_, DT_>(other),
          _rows(other._rows),
          _columns(other._columns),
          _zero_element(other._zero_element),
          _used_elements(other.used_elements())
        {
            _elements.insert(other._elements.begin(), other._elements.end());
        }

        template <typename Arch2_, typename DT2_>
        SparseMatrixCOO(const SparseMatrixCOO<Arch2_, DT2_> & other) :
          Container<Arch_, DT_>(other),
          _rows(other._rows),
          _columns(other._columns),
          _zero_element(other._zero_element),
          _used_elements(other.used_elements())
        {
            _elements.insert(other._elements.begin(), other._elements.end());
        }

        void operator()(Index row, Index col, DT_ val)
        {
          ASSERT(row < this->_rows, "Error: " + stringify(row) + "exceeds sparse matrix coo row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + "exceeds sparse matrix coo column size " + stringify(this->_columns) + " !");

          typename std::map<Index, DT_>::const_iterator it(_elements.find(row * _columns + col));
          if (it == _elements.end())
            ++_used_elements;

          _elements[row * _columns + col] = val;
        }

        const DT_ & operator()(Index row, Index col) const
        {
          ASSERT(row < this->_rows, "Error: " + stringify(row) + "exceeds sparse matrix coo row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + "exceeds sparse matrix coo column size " + stringify(this->_columns) + " !");

          typename std::map<Index, DT_>::const_iterator it(_elements.find(row * _columns + col));
          if (it == _elements.end())
            return _zero_element;
          else
            return it->second;
        }

        const Index & rows() const
        {
          return this->_rows;
        }

        const Index & columns() const
        {
          return this->_columns;
        }

        const Index & used_elements() const
        {
          return this->_used_elements;
        }

        const std::map<Index, DT_> & elements() const
        {
          return _elements;
        }

    };

    template <typename Arch_, typename DT_> bool operator== (const SparseMatrixCOO<Arch_, DT_> & a, const SparseMatrixCOO<Arch_, DT_> & b)
    {
      if (a.rows() != b.rows())
        return false;
      if (a.columns() != b.columns())
        return false;

      if (typeid(DT_) == typeid(Index))
      {
        for (Index i(0) ; i < a.rows() ; ++i)
        {
          for (Index j(0) ; j < a.columns() ; ++j)
          {
            if (a(i, j) != b(i, j))
              return false;
          }
        }
      }
      else
      {
        for (Index i(0) ; i < a.rows() ; ++i)
        {
          for (Index j(0) ; j < a.columns() ; ++j)
          {
            if (Absolute<DT_>::value(a(i, j) - b(i, j)) > std::numeric_limits<DT_>::epsilon())
              return false;
          }
        }
      }

      return true;
    }

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

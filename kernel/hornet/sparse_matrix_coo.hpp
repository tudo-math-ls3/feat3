#pragma once
#ifndef KERNEL_HORNET_DENSE_MATRIX_HPP
#define KERNEL_HORNET_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/hornet/container.hpp>

#include <vector>
#include <map>
#include <algorithm>


namespace FEAST
{
  template <typename Arch_, typename DT_>
  class SparseMatrixCOO : public Container<Arch_, DT_>
  {
    private:
      Index _rows;
      Index _columns;
      std::vector<std::map<unsigned long, DT_> > _elements;
      DT_ _zero_element;

    public:
      typedef DT_ data_type;

      SparseMatrixCOO(Index rows, Index columns) :
        Container<Arch_, DT_>(rows * columns),
        _elements(rows),
        _zero_element(0.)
      {
        this->_size = rows * columns;
        this->_rows = rows;
        this->_columns = columns;
      }

      SparseMatrixCOO(const SparseMatrixCOO<Arch_, DT_> & other) :
        Container<Arch_, DT_>(other),
        _rows(other._rows),
        _columns(other._columns),
        _elements(_rows),
        _zero_element(other._zero_element)
      {
        for (unsigned long i(0) ; i < _rows ; ++i)
          _elements.at(i).insert(other._elements.at(i).begin(), other._elements.at(i).end());
      }

      template <typename Arch2_, typename DT2_>
      SparseMatrixCOO(const SparseMatrixCOO<Arch2_, DT2_> & other) :
        Container<Arch_, DT_>(other),
        _rows(other._rows),
        _columns(other._columns),
        _elements(_rows),
        _zero_element(other._zero_element)
      {
        for (unsigned long i(0) ; i < _rows ; ++i)
          _elements.at(i).insert(other._elements.at(i).begin(), other._elements.at(i).end());
      }

      void operator()(Index row, Index col, DT_ val)
      {
        _elements.at(row)[col] = val;
      }

      const DT_ & operator()(Index row, Index col) const
      {
        typename std::map<unsigned long, DT_>::const_iterator it(_elements.at(row).find(col));
        if (it == _elements.at(row).end())
          return _zero_element;
        else
          return it->second;
      }

      virtual const Index & size() const
      {
        return this->_size;
      }

      virtual const Index & rows() const
      {
        return this->_rows;
      }

      virtual const Index & columns() const
      {
        return this->_columns;
      }
  };
} // namespace FEAST

#endif // KERNEL_HORNET_DENSE_VECTOR_HPP

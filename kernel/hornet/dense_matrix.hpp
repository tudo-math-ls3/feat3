#pragma once
#ifndef KERNEL_HORNET_DENSE_MATRIX_HPP
#define KERNEL_HORNET_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/hornet/container.hpp>


namespace FEAST
{
  template <typename Arch_, typename DT_>
  class DenseMatrix : public Container<Arch_, DT_>
  {
    private:
      DT_ * _pelements;
      Index _rows;
      Index _columns;

    public:
      typedef DT_ data_type;

      DenseMatrix(Index rows, Index columns) :
        Container<Arch_, DT_>(rows * columns)
      {
        this->_size = rows * columns;
        this->_rows = rows;
        this->_columns = columns;

        this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(this->_size * sizeof(DT_)));
        _pelements = this->_elements.at(0);
      }

      DenseMatrix(Index rows, Index columns, DT_ value) :
        Container<Arch_, DT_>(rows * columns)
      {
        this->_size = rows * columns;
        this->_rows = rows;
        this->_columns = columns;
        this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(this->_size * sizeof(DT_)));
        _pelements = this->_elements.at(0);

        //TODO add arch, use memory pool set memory function
        for (Index i(0) ; i < this->_size ; ++i)
        {
          _pelements[i] = value;
        }
      }

      DenseMatrix(const DenseMatrix<Arch_, DT_> & other) :
        Container<Arch_, DT_>(other),
        _rows(other._rows),
        _columns(other._columns)
      {
      }

      template <typename Arch2_, typename DT2_>
      DenseMatrix(const DenseMatrix<Arch2_, DT2_> & other) :
        Container<Arch_, DT_>(other)
      {
      }

      void operator()(Index index, DT_ val)
      {
        _pelements[index] = val;
      }

      void operator()(Index row, Index col, DT_ val)
      {
        _pelements[row * this->rows + col] = val;
      }

      const DT_ & operator[](Index index) const
      {
        return _pelements[index];
      }

      /*const DT_ & operator[](Index row, Index col) const
      {
        return _pelements[row * this->rows + col];
      }*/

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

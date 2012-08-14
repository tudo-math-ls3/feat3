#pragma once
#ifndef KERNEL_HORNET_DENSE_MATRIX_HPP
#define KERNEL_HORNET_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
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
        _pelements = this->_elements.at(0);
      }

      template <typename Arch2_, typename DT2_>
      DenseMatrix(const DenseMatrix<Arch2_, DT2_> & other) :
        Container<Arch_, DT_>(other),
        _rows(other._rows),
        _columns(other._columns)
      {
        _pelements = this->_elements.at(0);
      }

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

        std::vector<DT_ *> new_elements = other.get_elements();
        std::vector<unsigned long *> new_indices = other.get_indices();

        this->_elements.assign(new_elements.begin(), new_elements.end());
        this->_indices.assign(new_indices.begin(), new_indices.end());

        _pelements = this->_elements.at(0);

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

        return *this;
      }

      template <typename Arch2_, typename DT2_>
      DenseMatrix<Arch_, DT_> & operator= (const DenseMatrix<Arch2_, DT2_> & other)
      {
        if (this == &other)
          return *this;

        //TODO copy memory from arch2 to arch
        return *this;
      }

      void operator()(Index row, Index col, DT_ val)
      {
        ASSERT(row < this->_rows, "Error: " + stringify(row) + "exceeds dense matrix row size " + stringify(this->_rows) + " !");
        ASSERT(col < this->_columns, "Error: " + stringify(col) + "exceeds dense matrix column size " + stringify(this->_columns) + " !");
        _pelements[row * this->_rows + col] = val;
      }

      const DT_ & operator()(Index row, Index col) const
      {
        ASSERT(row < this->_rows, "Error: " + stringify(row) + "exceeds dense matrix row size " + stringify(this->_rows) + " !");
        ASSERT(col < this->_columns, "Error: " + stringify(col) + "exceeds dense matrix column size " + stringify(this->_columns) + " !");
        return _pelements[row * this->_rows + col];
      }

      const Index & rows() const
      {
        return this->_rows;
      }

      const Index & columns() const
      {
        return this->_columns;
      }
  };

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
} // namespace FEAST

#endif // KERNEL_HORNET_DENSE_VECTOR_HPP

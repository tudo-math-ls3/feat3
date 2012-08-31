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
    template <typename Arch_, typename DT_>
    class DenseMatrix : public Container<Arch_, DT_>
    {
      private:
        DT_ * _pelements;
        Index _rows;
        Index _columns;

      public:
        typedef DT_ data_type;
        typedef Arch_ arch_type;

        explicit DenseMatrix() :
          Container<Arch_, DT_> (0),
          _rows(0),
          _columns(0)
        {
        }

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
          _rows(other.rows()),
          _columns(other.columns())
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

        void operator()(Index row, Index col, DT_ value)
        {
          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->_columns) + " !");
          MemoryPool<Arch_>::modify_element(_pelements, row * this->_rows + col, value);
        }

        const DT_ operator()(Index row, Index col) const
        {
          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->_columns) + " !");
          return MemoryPool<Arch_>::get_element(_pelements, row * this->_rows + col);
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

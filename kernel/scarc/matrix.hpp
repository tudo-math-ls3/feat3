#pragma once
#ifndef SCARC_GUARD_MATRIX_HH
#define SCARC_GUARD_MATRIX_HH 1

#include<kernel/base_header.hpp>
#include<kernel/lafem/matrix_base.hpp>
#include<kernel/foundation/topology.hpp>

using namespace FEAST::Foundation;
using namespace FEAST::LAFEM;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {

    template<
      template<typename, typename> class StorageType_ = std::vector
      >
    class DynamicAOSMatrix : public MatrixBase
    {
      public:

        ///CTORS
        DynamicAOSMatrix(Index rows, Index columns) :
          _rows(rows),
          _columns(columns),
          _num_blocks(0),
          _block_topology(Topology<Index, StorageType_>()),
          _block_row_start(StorageType_<Index, std::allocator<Index> >()),
          _block_column_start(StorageType_<Index, std::allocator<Index> >()),
          _blocks(StorageType_<MatrixBase*, std::allocator<MatrixBase* > >())
        {
        }

        DynamicAOSMatrix(const DynamicAOSMatrix& other) :
          _rows(other._rows),
          _columns(other._columns),
          _num_blocks(other._num_blocks),
          _block_topology(other._block_topology),
          _block_row_start(other._block_row_start),
          _block_column_start(other._block_column_start),
          _blocks(other._blocks)
        {
        }

        ~DynamicAOSMatrix()
        {
        }

        DynamicAOSMatrix& operator=(const DynamicAOSMatrix& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_rows = rhs._rows;
          this->_columns = rhs._columns;
          this->_num_blocks = rhs._num_blocks;
          this->_block_topology = rhs._block_topology;
          this->_block_row_start = rhs._block_row_start;
          this->_block_column_start = rhs._block_column_start;
          this->_blocks = rhs._blocks;

          return *this;
        }

        const Index & rows() const
        {
          return _rows;
        }

        const Index & columns() const
        {
          return _columns;
        }

        const Index & num_blocks() const
        {
          return _num_blocks;
        }

        void add_block(MatrixBase* block, Index row_offset, Index column_offset)
        {
          _block_row_start.push_back(row_offset);
          _block_column_start.push_back(column_offset);
          _blocks.push_back(block);
          ++_num_blocks;
        }

        void remove_block(Index i)
        {
          _blocks.erase(i);
          _block_row_start.erase(i);
          _block_column_start.erase(i);
          --_num_blocks;
        }

      private:
        ///global dimensions
        Index _rows;
        Index _columns;
        Index _num_blocks;

        ///block topology: opt feature, TODO
        Topology<Index, StorageType_> _block_topology;

        ///block start
        StorageType_<Index, std::allocator<Index> > _block_row_start;
        StorageType_<Index, std::allocator<Index> > _block_column_start;

        ///blocks on process
        StorageType_<MatrixBase*, std::allocator<MatrixBase* > > _blocks;

        const MatrixBase& _get_block(Index i) const
        {
          return _blocks.at(i);
        }

        MatrixBase& _get_block(Index i)
        {
          return _blocks.at(i);
        }
    };


    template<typename T_>
    class OffProcessBlock : public MatrixBase
    {
      public:

        OffProcessBlock(Index rank, Index rows, Index columns) :
          _rows(rows),
          _columns(columns),
          _rank(rank)
        {
        }

        typedef T_ type_;

        const Index & rows() const
        {
          return _rows;
        }

        const Index & columns() const
        {
          return _columns;
        }

        const Index & rank() const
        {
          return _rank;
        }

      private:

        Index _rows;
        Index _columns;

        Index _rank;
    };
  }
}

#endif

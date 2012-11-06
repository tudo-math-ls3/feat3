#pragma once
#ifndef SCARC_GUARD_VECTOR_HH
#define SCARC_GUARD_VECTOR_HH 1

#include<kernel/base_header.hpp>
#include<kernel/lafem/vector_base.hpp>
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
    class DynamicVector : public VectorBase
    {
      public:

        ///CTORS
        DynamicVector(Index size) :
          _rows(size),
          _num_blocks(0),
          _block_topology(Topology<Index, StorageType_>()),
          _block_row_start(StorageType_<Index, std::allocator<Index> >()),
          _blocks(StorageType_<VectorBase*, std::allocator<VectorBase* > >())
        {
        }

        DynamicVector(const DynamicVector& other) :
          _rows(other._rows),
          _num_blocks(other._num_blocks),
          _block_topology(other._block_topology),
          _block_row_start(other._block_row_start),
          _blocks(other._blocks)
        {
        }

        ~DynamicVector()
        {
        }

        DynamicVector& operator=(const DynamicVector& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_rows = rhs._rows;
          this->_num_blocks = rhs._num_blocks;
          this->_block_topology = rhs._block_topology;
          this->_block_row_start = rhs._block_row_start;
          this->_blocks = rhs._blocks;

          return *this;
        }

        const Index & size() const
        {
          return _rows;
        }

        const Index & num_blocks() const
        {
          return _num_blocks;
        }

        ///row_offset = k => sizeof(block_{i-1}) = k
        void add_block(VectorBase& block, Index row_offset)
        {
          _block_row_start.push_back(row_offset);
          _blocks.push_back(&block);
          ++_num_blocks;
        }

        void remove_block(Index i)
        {
          _blocks.erase(i);
          _block_row_start.erase(i);
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

        ///blocks on process
        StorageType_<VectorBase*, std::allocator<VectorBase* > > _blocks;

        const VectorBase& _get_block(Index i) const
        {
          return *(_blocks.at(i));
        }

        VectorBase& _get_block(Index i)
        {
          return *(_blocks.at(i));
        }
    };
  }
}

#endif

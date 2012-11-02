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
      template<typename, typename> class StorageType_
      >
    class Matrix : public MatrixBase
    {
      public:
        Matrix() :
          _rows(0),
          _columns(0),
          _block_topology(Topology<Index, StorageType_>()),
          _blocks(StorageType_<std::shared_ptr<MatrixBase>, std::allocator<std::shared_ptr<MatrixBase> > >())
        {
        }

        const Index & rows() const
        {
          return _rows;
        }

        const Index & columns() const
        {
          return _columns;
        }

      private:
        Index rows;
        Index columns;
        Topology<Index, StorageType_> _block_topology;
        StorageType_<std::shared_ptr<MatrixBase>, std::allocator<std::shared_ptr<MatrixBase> > > _blocks;
    };
  }
}

#endif

#pragma once
#ifndef KERNEL_FOUNDATION_TOPOLOGY_HH
#define KERNEL_FOUNDATION_TOPOLOGY_HH 1

#include <vector>
#include<kernel/base_header.hpp>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief Topology is a dynamic, generic adjacency list
     *
     * long desc missing
     *
     * \tparam IndexType_
     * index type
     *
     * \tparam OuterStorageType__
     * dynamic storage type for storing the index lists
     *
     * \tparam StorageType_
     * type of index lists - for compatibility: storage types must only have size(), push_back() and [] members in common
     *
     *
     * \author Markus Geveler
     */
    template<
      typename IndexType_ = Index,
      template<typename, typename> class OuterStorageType_ = std::vector,
      typename StorageType_ = std::vector<IndexType_> >
    class Topology
    {
      public:
        typedef IndexType_ index_type_;
        typedef StorageType_ storage_type_;
        typedef OuterStorageType_<StorageType_, std::allocator<StorageType_> > compound_storage_type_;

        Topology() :
          _num_polytopes(0),
          _topology(OuterStorageType_<StorageType_, std::allocator<StorageType_> >())
      {
      };

        ~Topology()
        {
        }

        Index size()
        {
          return _num_polytopes;
        }

        void push_back(const StorageType_ s)
        {
          _topology.push_back(s);
          ++_num_polytopes;
        }

        StorageType_ & at(Index i)
        {
          return _topology.at(i);
        }

        StorageType_ & operator[] (Index i)
        {
          return _topology.at(i);
        }

        void push_back()
        {
          StorageType_ s;
          _topology.push_back(s);
          ++_num_polytopes;
        }

      private:
        Index _num_polytopes;
        OuterStorageType_<StorageType_, std::allocator<StorageType_> > _topology;
    };

  }
}
#endif

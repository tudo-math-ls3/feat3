#pragma once
#ifndef KERNEL_FOUNDATION_MAPPED_TOPOLOGY_HH
#define KERNEL_FOUNDATION_MAPPED_TOPOLOGY_HH 1

#include<stdexcept>
#include<kernel/base_header.hpp>
#include<kernel/foundation/topology.hpp>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief MappedTopology is a dynamic, generic Topology of indices with a map containing the associated objects of type ValueType_
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
      class ValueType_,
      typename IndexType_ = Index,
      template<typename, typename> class OuterStorageType_ = std::vector,
      typename StorageType_ = std::vector<IndexType_> >
    class MappedTopology : public Topology<IndexType_, OuterStorageType_, StorageType_>
    {
      public:
        ///CTOR
        MappedTopology() :
          Topology<IndexType_, OuterStorageType_, StorageType_>(),
          _keys(StorageType_()),
          _values(OuterStorageType_<ValueType_, std::allocator<ValueType_> >())
        {
        }

        ///DTOR
        ~MappedTopology()
        {
        }

        ///public access member functions
        void insert(const IndexType_ i, const ValueType_& v)
        {
          //insert with this key to map
          _keys.push_back(i);
          _values.push_back(v);

          //create new adjacencylist
          ((Topology<IndexType_, OuterStorageType_, StorageType_>*)this)->push_back();
        }

        void insert(const IndexType_ i, const ValueType_& v, const StorageType_& adjacencies)
        {
          //insert with this key to map
          _keys.push_back(i);
          _values.push_back(v);

          //take over adjacencylist
          ((Topology<IndexType_, OuterStorageType_, StorageType_>*)this)->push_back(adjacencies);
        }

        ValueType_ & find(IndexType_ key)
        {
          for(Index i(0) ; i < _keys.size() ; ++i)
          {
            if(_keys.at(i) == key)
              return _values.at(i);
          }
          throw std::out_of_range("Foundation::MappedTopology: Key not valid!");
        }

        ValueType_ & at(IndexType_ i)
        {
          return _values.at(i);
        }

        ValueType_ & operator[](IndexType_ i)
        {
          return _values.at(i);
        }

      private:
        StorageType_ _keys;
        OuterStorageType_<ValueType_, std::allocator<ValueType_> > _values;
    };

  }
}
#endif

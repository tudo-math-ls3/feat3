#pragma once
#ifndef KERNEL_FOUNDATION_TOPOLOGY_HH
#define KERNEL_FOUNDATION_TOPOLOGY_HH 1

#include<vector>
#include<kernel/base_header.hpp>
#include<kernel/foundation/functor.hpp>
#include<kernel/foundation/smart_pointer.hpp>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief Topology is a dynamic, generic adjacency list
     *
     * Topology is used in many ways:
     * (1) Contains polytope adjacencies for a set of polytopes without caring for the polytope level relation in top-level meshes.
     * (2) Describes Network topologies.
     * (3) Describes patch topologies.
     * (4) Describes top-level mesh topologies in parallel environment.
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
        ///type exports
        typedef IndexType_ index_type_;
        typedef StorageType_ storage_type_;
        typedef OuterStorageType_<StorageType_, std::allocator<StorageType_> > compound_storage_type_;

        ///CTOR
        Topology() :
          _num_polytopes(0),
          _topology(),
          _history()
        {
        };

        ///Copy CTOR
        Topology(const Topology& other) :
          _num_polytopes(other._num_polytopes),
          _topology(other._topology),
          _history()
        {
        };

        ///DTOR
        ~Topology()
        {
        }

        /**
         * \brief member function retrieves number of polytopes / size of the topology
         */
        Index size()
        {
          return _num_polytopes;
        }

        /**
         * \brief member function inserts a given polytope list to end of topology
         *
         * \param[in] s
         * The polytope list to be inserted
         */
        void push_back(const StorageType_& s)
        {
          _topology.push_back(s);
          _history.push_back(SmartPointer<FunctorBase>(new PushBackFunctor<compound_storage_type_, index_type_, StorageType_>(_topology, _num_polytopes, s)));
          ++_num_polytopes;
        }

        void erase(IndexType_ i)
        {
          _history.push_back(SmartPointer<FunctorBase>(new EraseFunctor<compound_storage_type_, index_type_, StorageType_>(_topology, i, _topology.at(i))));
          _topology.erase(_topology.begin() + i);
          --_num_polytopes;
        }

        /**
         * \brief member function retrieves polytope list for given polytope
         *
         * \param[in] i
         * The index of the polytope whose adjacency list is to be returned
         */
        StorageType_ & at(Index i)
        {
          return _topology.at(i);
        }

        /**
         * \brief operator overload to [] retrieves polytope list for given polytope
         *
         * \param[in] i
         * The index of the polytope whose adjacency list is to be returned
         */
        StorageType_ & operator[] (Index i)
        {
          return _topology.at(i);
        }

        /**
         * \brief member function inserts an empty list to end of topology
         *
         */
        void push_back()
        {
          StorageType_ s;
          _topology.push_back(s);
          _history.push_back(SmartPointer<FunctorBase>(new EmptyPushBackFunctor<compound_storage_type_, index_type_, storage_type_>(_topology, _num_polytopes, s)));
          ++_num_polytopes;
        }

        void erase()
        {
          _history.push_back(SmartPointer<FunctorBase>(new EmptyEraseFunctor<compound_storage_type_, index_type_, StorageType_>(_topology, _num_polytopes - 1, _topology.at(_num_polytopes - 1))));
          _topology.erase(_topology.end() - 1);
          --_num_polytopes;
        }

        Topology& operator=(const Topology& rhs)
        {
          if(this == &rhs)
            return *this;

          this-> _num_polytopes = rhs._num_polytopes;
          this-> _topology = rhs._topology;

          return *this;
        }

        OuterStorageType_<SmartPointer<FunctorBase>, std::allocator<SmartPointer<FunctorBase> > >& get_history()
        {
          return _history;
        }

        Index end()
        {
          return _num_polytopes;
        }

        Index begin()
        {
          return 0;
        }

      private:
        ///current size of the topology
        Index _num_polytopes;
        ///data
        OuterStorageType_<StorageType_, std::allocator<StorageType_> > _topology;
        ///history
        OuterStorageType_<SmartPointer<FunctorBase>, std::allocator<SmartPointer<FunctorBase> > > _history;
    };
  }
}
#endif

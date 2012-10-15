#pragma once
#ifndef KERNEL_FOUNDATION_TOPOLOGY_HH
#define KERNEL_FOUNDATION_TOPOLOGY_HH 1

#include<vector>
#include<kernel/base_header.hpp>
#include<kernel/foundation/functor.hpp>
#include<kernel/foundation/buffer.hpp>
#include<kernel/archs.hpp>

using namespace FEAST::Archs;

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
    class Topology :
      public Bufferable<BufferedData<OuterStorageType_> >,
      public Communicateable<BufferedData<OuterStorageType_> >
    {
      public:
        ///type exports
        typedef IndexType_ index_type_;
        typedef StorageType_ storage_type_;
        typedef OuterStorageType_<StorageType_, std::allocator<StorageType_> > compound_storage_type_;
        typedef Topology<IndexType_, OuterStorageType_, StorageType_> exact_type_;
        typedef BufferedData<OuterStorageType_> buffer_type_;

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

        const Index size() const
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
          _history.add_functor(new PushBackFunctor<compound_storage_type_, IndexType_, StorageType_>(_topology, _num_polytopes, s));
          _history.get_functors().at(_history.size() - 1).get()->execute();
          ++_num_polytopes;
        }

        void erase(IndexType_ i)
        {
          _history.add_functor(new EraseFunctor<compound_storage_type_, IndexType_, StorageType_>(_topology, i, _topology.at(i)));
          _history.get_functors().at(_history.size() - 1).get()->execute();
          --_num_polytopes;
        }

        ///insert value into list i (not an insert in the STL-sense)
        void insert(IndexType_ i, IndexType_ value)
        {
          _history.add_functor(new PushBackFunctor<storage_type_, IndexType_, IndexType_>(_topology.at(i), IndexType_(_topology.at(i).size()), value));
          _history.get_functors().at(_history.size() - 1).get()->execute();
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

        const StorageType_ & at(Index i) const
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

        const StorageType_ & operator[] (Index i) const
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
          _history.add_functor(new PushBackFunctor<compound_storage_type_, IndexType_, StorageType_>(_topology, _num_polytopes, s));
          _history.get_functors().at(_history.size() - 1).get()->execute();
          ++_num_polytopes;
        }

        void erase()
        {
          _history.add_functor(new EraseFunctor<compound_storage_type_, IndexType_, StorageType_>(_topology, _num_polytopes - 1, _topology.at(_num_polytopes - 1)));
          _history.get_functors().at(_history.size() - 1).get()->execute();
          --_num_polytopes;
        }

        Topology& operator=(const Topology& rhs)
        {
          if(this == &rhs)
            return *this;

          this-> _num_polytopes = rhs._num_polytopes;
          this-> _topology = rhs._topology;
          this-> _history = CompoundFunctor<OuterStorageType_>();

          return *this;
        }

        CompoundFunctor<OuterStorageType_>& get_history()
        {
          return _history;
        }

        const CompoundFunctor<OuterStorageType_>& get_history() const
        {
          return _history;
        }

        OuterStorageType_<StorageType_, std::allocator<StorageType_> >& get_topology()
        {
          return _topology;
        }

        const OuterStorageType_<StorageType_, std::allocator<StorageType_> >& get_topology() const
        {
          return _topology;
        }

        ///implementation of Bufferable interface
        virtual BufferedData<OuterStorageType_> buffer(Index estimated_size_increase = 0)
        {
          BufferedData<OuterStorageType_> result;

          result.get().push_back(BufferedSharedArray<IndexType_>::create(3)); //sizes, 'row pointers', data
          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[0] = 3;

          //how many row-pointers?
          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[1] = _topology.size() + estimated_size_increase;
          result.get().push_back(BufferedSharedArray<IndexType_>::create(_topology.size() + estimated_size_increase));

          IndexType_ final_size(0);
          for(IndexType_ i(0) ; i < _topology.size() ; ++i)
          {
            (*(BufferedSharedArray<IndexType_>*)((result.get().at(1).get())))[i] = _topology.at(i).size() + estimated_size_increase;
            final_size += _topology.at(i).size() + estimated_size_increase;
          }
          final_size *= _topology.size() + estimated_size_increase;

          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[2] = final_size;;

          result.get().push_back(BufferedSharedArray<IndexType_>::create(final_size));

          return result;
        }

        virtual void to_buffer(BufferedData<OuterStorageType_>& buffer)
        {
          IndexType_ head(0);

          for(IndexType_ i(0) ; i < _topology.size() ; ++i)
          {
            for(IndexType_ j(0) ; j < _topology.at(i).size() ; ++j)
            {
              (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(2).get())))[head] = _topology.at(i).at(j);
              ++head;
            }
          }
        }

        virtual void from_buffer(const BufferedData<OuterStorageType_>& buffer)
        {
          _topology.clear();
          _history.get_functors().clear();

          IndexType_ left(0);
          IndexType_ head(0);
          for(Index i(0) ; i < (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(0).get())))[1] ; ++i) //times the real #polytopes
          {
            push_back();
            for(Index j(0) ; j < (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(1).get())))[i]; ++j)
            {
              _topology.at(i).push_back((*(BufferedSharedArray<IndexType_>*)((buffer.get().at(2).get())))[head]);
              ++head;
            }
          }
        }

        ///implementation of Bufferable interface
        void send_recv(BufferedData<OuterStorageType_>& sendbuffers,
                       int destrank,
                       BufferedData<OuterStorageType_>& recvbuffers,
                       int sourcerank)
        {
#ifndef FEAST_SERIAL_MODE

          //get sizes
          Comm<Parallel>::send_recv(((BufferedSharedArray<Index>*)sendbuffers.get().at(0).get())->get(),
                                    (*(BufferedSharedArray<Index>*)((sendbuffers.get().at(0).get())))[0],
                                    destrank,
                                    ((BufferedSharedArray<Index>*)recvbuffers.get().at(0).get())->get(),
                                    (*(BufferedSharedArray<Index>*)((recvbuffers.get().at(0).get())))[0],
                                    sourcerank);

          //get row_ptrs
          IndexType_ recv_num_polytopes((*(BufferedSharedArray<Index>*)((recvbuffers.get().at(0).get())))[1]);

          Comm<Parallel>::send_recv(((BufferedSharedArray<Index>*)sendbuffers.get().at(1).get())->get(),
                                    (*(BufferedSharedArray<Index>*)((sendbuffers.get().at(0).get())))[1],
                                    destrank,
                                    ((BufferedSharedArray<Index>*)recvbuffers.get().at(1).get())->get(),
                                    recv_num_polytopes,
                                    sourcerank);
          //get data
          IndexType_ recv_datasize((*(BufferedSharedArray<Index>*)((recvbuffers.get().at(0).get())))[2]);

          Comm<Parallel>::send_recv(((BufferedSharedArray<Index>*)sendbuffers.get().at(2).get())->get(),
                                    (*(BufferedSharedArray<Index>*)((sendbuffers.get().at(0).get())))[2],
                                    destrank,
                                    ((BufferedSharedArray<Index>*)recvbuffers.get().at(2).get())->get(),
                                    recv_datasize,
                                    sourcerank);
#else
          ///TODO
#endif
        }

      private:
        ///current size of the topology
        Index _num_polytopes;
        ///data
        OuterStorageType_<StorageType_, std::allocator<StorageType_> > _topology;
        ///history
        CompoundFunctor<OuterStorageType_> _history;
    };
  }
}
#endif

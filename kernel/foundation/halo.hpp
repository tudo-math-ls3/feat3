#pragma once
#ifndef KERNEL_FOUNDATION_HALO_HH
#define KERNEL_FOUNDATION_HALO_HH 1

#include<vector>
#include<kernel/foundation/mesh.hpp>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/buffer.hpp>
#include<kernel/archs.hpp>

using namespace FEAST::Archs;

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief Halo class template represents the interface between two meshes or patches of any kind
     *
     * A Halo is used for communication between patches of any kind. Currently, its policy is inspired by mesh patches.
     * The overlap parameter is used in communication routines to destinct between minimal overlap and element layer overlap.
     *
     * \tparam delta_
     * overlap parameter: delta = 0 => minimal overlap, delta = k > 0 => element layer overlap with k layers.
     *
     * \tparam MeshType_
     * type of the mesh or patch
     *
     * \tparam StorageType_
     * type for inner storage of indices
     *
     * \tparam IndexType_
     * type of the indices
     *
     * \author Markus Geveler
     */
    template<unsigned delta_,
             PolytopeLevels level_,
             typename MeshType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IndexType_ = Index>
    class Halo :
      public Communicateable<BufferedData<StorageType_> >,
      public Bufferable<BufferedData<StorageType_> >
    {
      public:
        ///type exports:
        typedef IndexType_ index_type_;
        typedef MeshType_ mesh_type_;
        typedef BufferedData<StorageType_> buffer_type_;

        ///CTOR
        Halo(MeshType_ & mesh, IndexType_ other, PolytopeLevels level = pl_edge) : //TODO move to template
          _halo_elements(),
          _halo_element_counterparts(),
          _mesh(mesh),
          _other(other),
          _overlap(delta_),
          _level(level)
        {
        }

        ///CTOR from buffer
        template<typename BufferType_>
        Halo(const BufferType_& buffer, MeshType_& new_mesh, IndexType_ new_other_rank) :
          _halo_elements(),
          _halo_element_counterparts(),
          _mesh(new_mesh),
          _other(new_other_rank),
          _overlap(delta_),
          _level(level_)
        {
          for(IndexType_ i(0) ; i < buffer.size() ; ++i)
          {
            _halo_elements.push_back(buffer.get_element(i));
            _halo_element_counterparts.push_back(buffer.get_element_counterpart(i));
          }
        }

        ///DTOR
        ~Halo()
        {
        }

        ///Add correspondence of i
        void add_element_pair(IndexType_ i, IndexType_ j)
        {
          _halo_elements.push_back(i);
          _halo_element_counterparts.push_back(j);
        }

        ///public access functions:
        IndexType_ get_element_counterpart(IndexType_ index)
        {
          return _halo_element_counterparts.at(index);
        }

        IndexType_ get_element(IndexType_ index)
        {
          return _halo_elements.at(index);
        }

        IndexType_ size()
        {
          return IndexType_(_halo_elements.size());
        }

        MeshType_ & get_mesh()
        {
          return _mesh;
        }

        IndexType_ get_other()
        {
          return _other;
        }

        unsigned get_overlap()
        {
          return _overlap;
        }

        PolytopeLevels get_level()
        {
          return _level;
        }

        Halo& operator=(const Halo& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_halo_elements = rhs._halo_elements;
          this->_halo_element_counterparts = rhs._halo_element_counterparts;
          this->_mesh = rhs._mesh;
          this->_other = rhs._other;
          this->_overlap = rhs._overlap;
          this->_level = rhs._level;

          return *this;
        }

        ///implementation of Bufferable interface
        virtual BufferedData<StorageType_> buffer()
        {
          BufferedData<StorageType_> result;
          result.get().push_back(BufferedSharedArray<IndexType_>::create(3));
          result.get().push_back(BufferedSharedArray<IndexType_>::create(_halo_elements.size()));
          result.get().push_back(BufferedSharedArray<IndexType_>::create(_halo_element_counterparts.size()));

          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[0] = 3;
          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[1] = _halo_elements.size();
          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[2] = _halo_element_counterparts.size();

          return result;
        }

        virtual void to_buffer(BufferedData<StorageType_>& buffer)
        {
          for(IndexType_ i(0) ; i < _halo_elements.size() ; ++i)
          {
            (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(1).get())))[i] = _halo_elements.at(i);
          }

          for(IndexType_ i(0) ; i < _halo_element_counterparts.size() ; ++i)
          {
            (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(2).get())))[i] = _halo_element_counterparts.at(i);
          }
        }

        virtual void from_buffer(const BufferedData<StorageType_>& buffer)
        {
          _halo_elements.clear();
          _halo_element_counterparts.clear();

          for(IndexType_ i(0) ; i < (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(0).get())))[0] ; ++i)
          {
            _halo_elements.push_back( (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(1).get())))[i] );
          }

          for(IndexType_ i(0) ; i < (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(0).get())))[1] ; ++i)
          {
            _halo_element_counterparts.push_back( (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(2).get())))[i] );
          }
        }

        ///implementation of Communicateable interface
        void send_recv(BufferedData<StorageType_>& sendbuffers,
                       int destrank,
                       BufferedData<StorageType_>& recvbuffers,
                       int sourcerank)
        {
#ifndef FEAST_SERIAL_MODE
          for(IndexType_ i(0) ; i < sendbuffers.get().size() ; ++i)
          {
            Comm<Parallel>::send_recv(((BufferedSharedArray<IndexType_>*)sendbuffers.get().at(i).get())->get(),
                                      (*(BufferedSharedArray<IndexType_>*)((sendbuffers.get().at(0).get())))[i],
                                      destrank,
                                      ((BufferedSharedArray<IndexType_>*)recvbuffers.get().at(i).get())->get(),
                                      (*(BufferedSharedArray<IndexType_>*)((recvbuffers.get().at(0).get())))[i],
                                      sourcerank);
          }
#else
          for(IndexType_ i(0) ; i < sendbuffers.get().size() ; ++i)
          {
            Comm<Serial>::send_recv(((BufferedSharedArray<IndexType_>*)sendbuffers.get().at(i).get())->get(),
                                      (*(BufferedSharedArray<IndexType_>*)((sendbuffers.get().at(0).get())))[i],
                                      destrank,
                                      ((BufferedSharedArray<IndexType_>*)recvbuffers.get().at(i).get())->get(),
                                      (*(BufferedSharedArray<IndexType_>*)((recvbuffers.get().at(0).get())))[i],
                                      sourcerank);
          }
#endif
        }

      private:
        StorageType_<IndexType_, std::allocator<IndexType_> > _halo_elements;
        StorageType_<IndexType_, std::allocator<IndexType_> > _halo_element_counterparts;

        MeshType_ & _mesh;
        IndexType_ _other;

        unsigned _overlap;

        PolytopeLevels _level;
    };
  }
}


#endif

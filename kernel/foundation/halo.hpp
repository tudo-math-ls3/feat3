#pragma once
#ifndef KERNEL_FOUNDATION_HALO_HH
#define KERNEL_FOUNDATION_HALO_HH 1

#include<vector>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/buffer.hpp>
#include<kernel/archs.hpp>

using namespace FEAST::Archs;

namespace FEAST
{
  namespace Foundation
  {
    template<typename MeshType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IndexType_ = Index>
    class HaloBase :
      public Bufferable<BufferedData<StorageType_> >,
      public Communicateable<BufferedData<StorageType_>, com_send_receive>
    {
      public:
        ///type exports:
        typedef IndexType_ index_type_;
        typedef MeshType_ mesh_type_;
        typedef BufferedData<StorageType_> buffer_type_;

        ///CTOR
        HaloBase(MeshType_ & mesh, IndexType_ other = 0) : //TODO move to template
          _halo_elements(),
          _halo_element_counterparts(),
          _mesh(mesh),
          _other(other)
        {
        }

        ///CTOR from buffer
        template<typename BufferType_>
        HaloBase(const BufferType_& buffer, MeshType_& new_mesh, IndexType_ new_other_rank) :
          _halo_elements(),
          _halo_element_counterparts(),
          _mesh(new_mesh),
          _other(new_other_rank)
        {
          for(IndexType_ i(0) ; i < buffer.size() ; ++i)
          {
            _halo_elements.push_back(buffer.get_element(i));
            _halo_element_counterparts.push_back(buffer.get_element_counterpart(i));
          }
        }

        ///copy-CTOR
        HaloBase(const HaloBase& other) : //TODO move to template
          _halo_elements(other._halo_elements),
          _halo_element_counterparts(other._halo_element_counterparts),
          _mesh(other._mesh),
          _other(other._other)
        {
        }

        ///DTOR
        virtual ~HaloBase()
        {
        }

        ///Add correspondence of i
        virtual void add_element_pair(IndexType_ i, IndexType_ j)
        {
          _halo_elements.push_back(i);
          _halo_element_counterparts.push_back(j);
        }

        ///public access functions:
        virtual IndexType_ get_element_counterpart(IndexType_ index) const
        {
          return _halo_element_counterparts.at(index);
        }

        virtual IndexType_ get_element(IndexType_ index) const
        {
          return _halo_elements.at(index);
        }

        virtual IndexType_ size() const
        {
          return IndexType_(_halo_elements.size());
        }

        virtual const MeshType_ & get_mesh() const
        {
          return _mesh;
        }

        virtual MeshType_ & get_mesh()
        {
          return _mesh;
        }

        virtual IndexType_ get_other() const
        {
          return _other;
        }

        virtual unsigned get_overlap() const = 0;
        virtual PolytopeLevels get_level() const = 0;

        virtual StorageType_<IndexType_, std::allocator<IndexType_> >& get_elements()
        {
          return _halo_elements;
        }

        virtual const StorageType_<IndexType_, std::allocator<IndexType_> >& get_elements() const
        {
          return _halo_elements;
        }

        virtual StorageType_<IndexType_, std::allocator<IndexType_> >& get_element_counterparts()
        {
          return _halo_element_counterparts;
        }

        virtual const StorageType_<IndexType_, std::allocator<IndexType_> >& get_element_counterparts() const
        {
          return _halo_element_counterparts;
        }

        ///implementation of Bufferable interface
        virtual BufferedData<StorageType_> buffer(IndexType_ estimated_size_increase = 0)
        {
          BufferedData<StorageType_> result;
          result.get().push_back(BufferedSharedArray<IndexType_>::create(3));
          result.get().push_back(BufferedSharedArray<IndexType_>::create(_halo_elements.size() + estimated_size_increase));
          result.get().push_back(BufferedSharedArray<IndexType_>::create(_halo_element_counterparts.size() + estimated_size_increase));

          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[0] = 3;
          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[1] = _halo_elements.size() + estimated_size_increase;
          (*(BufferedSharedArray<IndexType_>*)((result.get().at(0).get())))[2] = _halo_element_counterparts.size() + estimated_size_increase;

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

          for(IndexType_ i(0) ; i < (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(0).get())))[1] ; ++i)
          {
            _halo_elements.push_back( (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(1).get())))[i] );
          }

          for(IndexType_ i(0) ; i < (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(0).get())))[2] ; ++i)
          {
            _halo_element_counterparts.push_back( (*(BufferedSharedArray<IndexType_>*)((buffer.get().at(2).get())))[i] );
          }
        }

        ///implementation of Communicateable interface
        virtual void send_recv(BufferedData<StorageType_>& sendbuffers,
                       int destrank,
                       BufferedData<StorageType_>& recvbuffers,
                       int sourcerank)
        {
#ifndef SERIAL
          Comm<Parallel>::send_recv(((BufferedSharedArray<IndexType_>*)sendbuffers.get().at(0).get())->get(),
              3,
              destrank,
              ((BufferedSharedArray<IndexType_>*)recvbuffers.get().at(0).get())->get(),
              3,
              sourcerank);

          for(IndexType_ i(1) ; i < 3 ; ++i)
          {
            Comm<Parallel>::send_recv(((BufferedSharedArray<IndexType_>*)sendbuffers.get().at(i).get())->get(),
                                      (*(BufferedSharedArray<IndexType_>*)((sendbuffers.get().at(0).get())))[i],
                                      destrank,
                                      ((BufferedSharedArray<IndexType_>*)recvbuffers.get().at(i).get())->get(),
                                      (*(BufferedSharedArray<IndexType_>*)((recvbuffers.get().at(0).get())))[i],
                                      sourcerank);
          }
#else
          Comm<Serial>::send_recv(((BufferedSharedArray<IndexType_>*)sendbuffers.get().at(0).get())->get(),
              3,
              destrank,
              ((BufferedSharedArray<IndexType_>*)recvbuffers.get().at(0).get())->get(),
              3,
              sourcerank);

          for(IndexType_ i(1) ; i < 3 ; ++i)
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

      protected:
        StorageType_<IndexType_, std::allocator<IndexType_> > _halo_elements;
        StorageType_<IndexType_, std::allocator<IndexType_> > _halo_element_counterparts;

        MeshType_ & _mesh;
        IndexType_ _other;

    };

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
      public HaloBase<MeshType_, StorageType_, IndexType_>
    {
      public:
        ///type exports:
        typedef typename HaloBase<MeshType_, StorageType_, IndexType_>::index_type_ index_type_;
        typedef typename HaloBase<MeshType_, StorageType_, IndexType_>::mesh_type_ mesh_type_;
        typedef typename HaloBase<MeshType_, StorageType_, IndexType_>::buffer_type_ buffer_type_;

        ///CTOR
        Halo(MeshType_ & mesh, IndexType_ other = 0) : //TODO move to template
          HaloBase<MeshType_, StorageType_, IndexType_>(mesh, other),
          _overlap(delta_),
          _level(level_)
        {
        }

        ///CTOR from buffer
        template<typename BufferType_>
        Halo(const BufferType_& buffer, MeshType_& new_mesh, IndexType_ new_other_rank) :
          HaloBase<MeshType_, StorageType_, IndexType_>(buffer, new_mesh, new_other_rank),
          _overlap(delta_),
          _level(level_)
        {
          for(IndexType_ i(0) ; i < buffer.size() ; ++i)
          {
            this->_halo_elements.push_back(buffer.get_element(i));
            this->_halo_element_counterparts.push_back(buffer.get_element_counterpart(i));
          }
        }

        ///copy-CTOR
        Halo(const Halo& other) : //TODO move to template
          HaloBase<MeshType_, StorageType_, IndexType_>(other),
          _overlap(other._overlap),
          _level(other._level)
        {
        }

        virtual unsigned get_overlap() const
        {
          return _overlap;
        }

        virtual PolytopeLevels get_level() const
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

      private:
        unsigned _overlap;

        PolytopeLevels _level;
    };


    template<typename MeshType_,
             template<typename, typename> class StorageType_,
             typename IndexType_>
    bool compare_other(const std::shared_ptr<HaloBase<MeshType_, StorageType_, IndexType_> >& l, const std::shared_ptr<HaloBase<MeshType_, StorageType_, IndexType_> >& r)
    {
      return (l->get_other() < r->get_other());
    }

  }
}


#endif

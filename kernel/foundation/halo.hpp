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
             template<typename, typename> class BufferStorageType_ = std::vector>
    class HaloBase :
      public Bufferable<BufferedData<BufferStorageType_> >,
      public Communicateable<BufferedData<BufferStorageType_>, com_send_receive>
    {
      public:
        ///type exports:
        typedef typename MeshType_::topology_type_::storage_type_ compound_storage_type_;
        typedef typename MeshType_::topology_type_::index_type_ index_type_;
        typedef MeshType_ mesh_type_;
        typedef BufferedData<BufferStorageType_> buffer_type_;

        ///CTOR
        HaloBase() :
          _halo_elements(),
          _mesh(nullptr),
          _other(0)
        {
        }

        ///CTOR
        HaloBase(MeshType_ & mesh, index_type_ other = 0) : //TODO move to template
          _halo_elements(),
          _mesh(&mesh),
          _other(other)
        {
        }

        ///CTOR from buffer
        template<typename BufferType_>
        HaloBase(const BufferType_& buffer, MeshType_& new_mesh, index_type_ new_other_rank) :
          _halo_elements(),
          _mesh(&new_mesh),
          _other(new_other_rank)
        {
          for(index_type_ i(0) ; i < buffer.size() ; ++i)
          {
            _halo_elements.push_back(buffer.get_element(i));
          }
        }

        ///copy-CTOR
        HaloBase(const HaloBase& other) : //TODO move to template
          _halo_elements(other._halo_elements),
          _mesh(other._mesh),
          _other(other._other)
        {
        }

        ///DTOR
        virtual ~HaloBase()
        {
        }

        ///Add correspondence of i
        virtual void push_back(const index_type_ i)
        {
          _halo_elements.push_back(i);
        }

        ///public access functions:

        virtual index_type_ get_element(const index_type_ index) const
        {
          return _halo_elements.at(index);
        }

        virtual void erase_element(const index_type_ index)
        {
          _halo_elements.erase(_halo_elements.begin() + (const typename compound_storage_type_::difference_type)index);
        }

        virtual index_type_ size() const
        {
          return index_type_(_halo_elements.size());
        }

        virtual const mesh_type_* get_mesh() const
        {
          return _mesh;
        }

        virtual mesh_type_* get_mesh()
        {
          return _mesh;
        }

        virtual void reset_mesh(mesh_type_* mp)
        {
          _mesh = mp;
        }

        virtual index_type_ get_other() const
        {
          return _other;
        }

        virtual void reset_other(index_type_ i)
        {
          _other = i;
        }

        virtual unsigned get_overlap() const = 0;
        virtual PolytopeLevels get_level() const = 0;

        virtual compound_storage_type_& get_elements()
        {
          return _halo_elements;
        }

        virtual const compound_storage_type_& get_elements() const
        {
          return _halo_elements;
        }

        ///implementation of Bufferable interface
        virtual BufferedData<BufferStorageType_> buffer(index_type_ estimated_size_increase = 0)
        {
          BufferedData<BufferStorageType_> result;
          result.get().push_back(BufferedSharedArray<index_type_>::create(3));
          result.get().push_back(BufferedSharedArray<index_type_>::create((index_type_)(_halo_elements.size()) + estimated_size_increase));

          (*(BufferedSharedArray<index_type_>*)((result.get().at(0).get())))[0] = 2;
          (*(BufferedSharedArray<index_type_>*)((result.get().at(0).get())))[1] = (index_type_)(_halo_elements.size()) + estimated_size_increase;

          return result;
        }

        virtual void to_buffer(BufferedData<BufferStorageType_>& buffer)
        {
          for(index_type_ i(0) ; i < _halo_elements.size() ; ++i)
          {
            (*(BufferedSharedArray<index_type_>*)((buffer.get().at(1).get())))[i] = _halo_elements.at(i);
          }
        }

        virtual void from_buffer(const BufferedData<BufferStorageType_>& buffer)
        {
          _halo_elements.clear();

          for(index_type_ i(0) ; i < (*(BufferedSharedArray<index_type_>*)((buffer.get().at(0).get())))[1] ; ++i)
          {
            _halo_elements.push_back( (*(BufferedSharedArray<index_type_>*)((buffer.get().at(1).get())))[i] );
          }
        }

        ///implementation of Communicateable interface
        virtual void send_recv(BufferedData<BufferStorageType_>& sendbuffers,
                       Index destrank,
                       BufferedData<BufferStorageType_>& recvbuffers,
                       Index sourcerank)
        {
          Comm::send_recv(((BufferedSharedArray<index_type_>*)sendbuffers.get().at(0).get())->get(),
              2,
              destrank,
              ((BufferedSharedArray<index_type_>*)recvbuffers.get().at(0).get())->get(),
              2,
              sourcerank);

          Comm::send_recv(((BufferedSharedArray<index_type_>*)sendbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<index_type_>*)((sendbuffers.get().at(0).get())))[1],
              destrank,
              ((BufferedSharedArray<index_type_>*)recvbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<index_type_>*)((recvbuffers.get().at(0).get())))[1],
              sourcerank);
        }

      protected:
        compound_storage_type_ _halo_elements;

        mesh_type_* _mesh;
        index_type_ _other;

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
     * \author Markus Geveler
     */
    template<unsigned delta_,
             typename Level_,
             typename MeshType_,
             template<typename, typename> class BufferStorageType_ = std::vector>
    class Halo :
      public HaloBase<MeshType_, BufferStorageType_>
    {
      public:
        ///type exports:
        typedef typename HaloBase<MeshType_, BufferStorageType_>::index_type_ index_type_;
        typedef typename HaloBase<MeshType_, BufferStorageType_>::mesh_type_ mesh_type_;
        typedef typename HaloBase<MeshType_, BufferStorageType_>::buffer_type_ buffer_type_;

        ///CTOR
        Halo() :
          HaloBase<MeshType_, BufferStorageType_>(),
          _overlap(delta_),
          _level(Level_::tag_value)
        {
        }

        ///CTOR
        Halo(MeshType_ & mesh, index_type_ other = 0) : //TODO move to template
          HaloBase<MeshType_, BufferStorageType_>(mesh, other),
          _overlap(delta_),
          _level(Level_::tag_value)
        {
        }

        ///CTOR from buffer
        template<typename BufferType_>
        Halo(const BufferType_& buffer, MeshType_& new_mesh, index_type_ new_other_rank) :
          HaloBase<MeshType_, BufferStorageType_>(buffer, new_mesh, new_other_rank),
          _overlap(delta_),
          _level(Level_::tag_value)
        {
          for(index_type_ i(0) ; i < buffer.size() ; ++i)
          {
            this->_halo_elements.push_back(buffer.get_element(i));
          }
        }

        ///copy-CTOR
        Halo(const Halo& other) : //TODO move to template
          HaloBase<MeshType_, BufferStorageType_>(other),
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
             template<typename, typename> class StorageType_>
    bool compare_other(const std::shared_ptr<HaloBase<MeshType_, StorageType_> >& l, const std::shared_ptr<HaloBase<MeshType_, StorageType_> >& r)
    {
      return (l->get_other() < r->get_other());
    }

  }
}
#endif

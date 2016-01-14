#pragma once
#ifndef KERNEL_FOUNDATION_HALO_HH
#define KERNEL_FOUNDATION_HALO_HH 1

#include<vector>
#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/buffer.hpp>
#include<kernel/archs.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template<typename MeshType_,
             typename WT_ = double,
             template<typename, typename> class BufferStorageType_ = std::vector>
    class HaloBase :
      public Bufferable<BufferedData<BufferStorageType_> >,
      public Communicateable<BufferedData<BufferStorageType_>, com_send_receive>
    {
      public:

        ///CommLink nested type
        template<typename WeightType_ = WT_>
        class CommLink
        {
          public:

            typedef WT_ weight_data_type_;

            CommLink() :
              _rank(typename MeshType_::index_type_(0)),
              _weight(WT_(0.5)),
              _complete(false),
              _global_id(0)
            {
            }

            CommLink(typename MeshType_::index_type_ rank, WeightType_ w) :
              _rank(rank),
              _weight(w),
              _complete(false),
              _global_id(0)
            {
            }

            CommLink(typename MeshType_::index_type_ rank, WeightType_ w, typename MeshType_::index_type_ id) :
              _rank(rank),
              _weight(w),
              _complete(true),
              _global_id(id)
            {
            }

            CommLink(const CommLink& other) :
              _rank(other._rank),
              _weight(other._weight),
              _complete(other._complete),
              _global_id(other._global_id)
            {
            }

            typename MeshType_::index_type_ get_rank() const
            {
              return _rank;
            }

            typename MeshType_::index_type_& get_rank()
            {
              return _rank;
            }

            WeightType_ weight() const
            {
              return _weight;
            }

            WeightType_& weight()
            {
              return _weight;
            }

            bool complete()
            {
              return _complete;
            }

            void set_complete()
            {
              _complete = true;
            }

            typename MeshType_::index_type_ get_global_id() const
            {
              return _global_id;
            }

            typename MeshType_::index_type_& get_global_id()
            {
              return _global_id;
            }

          private:
            typename MeshType_::index_type_ _rank;
            WeightType_ _weight;
            bool _complete;
            typename MeshType_::index_type_ _global_id;
        };

        ///type exports:
        typedef CommLink<WT_> comm_link_type_;
        typedef typename MeshType_::topology_type_::storage_type_ compound_storage_type_;
        typedef typename MeshType_::topology_type_::index_type_ index_type_;
        typedef MeshType_ mesh_type_;
        typedef BufferedData<BufferStorageType_> buffer_type_;

        ///CTOR
        HaloBase() :
          _halo_elements(),
          _mesh(nullptr),
          _cl()
        {
        }

        ///CTOR
        HaloBase(MeshType_ & mesh, index_type_ other = 0) : //TODO move to template
          _halo_elements(),
          _mesh(&mesh),
          _cl(other, WT_(0.5))
        {
        }

        ///CTOR from buffer
        template<typename BufferType_>
        HaloBase(const BufferType_& b, MeshType_& new_mesh, index_type_ new_other_rank, WT_ w = WT_(0.5)) :
          _halo_elements(),
          _mesh(&new_mesh),
          _cl(new_other_rank, w)
        {
          for(index_type_ i(0) ; i < b.size() ; ++i)
          {
            _halo_elements.push_back(b.get_element(i));
          }
        }

        ///copy-CTOR
        HaloBase(const HaloBase& other) : //TODO move to template
          _halo_elements(other._halo_elements),
          _mesh(other._mesh),
          _cl(other._cl)
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
          _halo_elements.erase(_halo_elements.begin() + (typename compound_storage_type_::difference_type)index);
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
          return _cl.get_rank();
        }

        virtual void reset_other(index_type_ i)
        {
          _cl.get_rank() = i;
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

        virtual void set_elements(compound_storage_type_&& data)
        {
          _halo_elements = data;
        }

        virtual Index get_global_id()
        {
          return _cl.get_global_id();
        }

        virtual void set_id(typename MeshType_::index_type_ id)
        {
          _cl.get_global_id() = id;
          _cl.set_complete();
        }

        ///implementation of Bufferable interface
        virtual BufferedData<BufferStorageType_> buffer(index_type_ estimated_size_increase = 0) override
        {
          BufferedData<BufferStorageType_> result;
          result.get().push_back(BufferedSharedArray<index_type_>::create(3));
          result.get().push_back(BufferedSharedArray<index_type_>::create((index_type_)(_halo_elements.size()) + estimated_size_increase));

          (*(BufferedSharedArray<index_type_>*)((result.get().at(0).get())))[0] = 2;
          (*(BufferedSharedArray<index_type_>*)((result.get().at(0).get())))[1] = (index_type_)(_halo_elements.size()) + estimated_size_increase;

          return result;
        }

        virtual void to_buffer(BufferedData<BufferStorageType_>& b) override
        {
          for(index_type_ i(0) ; i < _halo_elements.size() ; ++i)
          {
            (*(BufferedSharedArray<index_type_>*)((b.get().at(1).get())))[i] = _halo_elements.at(i);
          }
        }

        virtual void from_buffer(const BufferedData<BufferStorageType_>& b) override
        {
          _halo_elements.clear();

          for(index_type_ i(0) ; i < (*(BufferedSharedArray<index_type_>*)((b.get().at(0).get())))[1] ; ++i)
          {
            _halo_elements.push_back( (*(BufferedSharedArray<index_type_>*)((b.get().at(1).get())))[i] );
          }
        }

        ///implementation of Communicateable interface
        virtual void send_recv(BufferedData<BufferStorageType_>& sendbuffers,
                       typename MeshType_::index_type_ destrank,
                       BufferedData<BufferStorageType_>& recvbuffers,
                       typename MeshType_::index_type_ sourcerank) override
        {
          Status s1;
          Comm::send_recv(((BufferedSharedArray<index_type_>*)sendbuffers.get().at(0).get())->get(),
              2,
              destrank,
              ((BufferedSharedArray<index_type_>*)recvbuffers.get().at(0).get())->get(),
              2,
              sourcerank,
              s1);

          Status s2;
          Comm::send_recv(((BufferedSharedArray<index_type_>*)sendbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<index_type_>*)((sendbuffers.get().at(0).get())))[1],
              destrank,
              ((BufferedSharedArray<index_type_>*)recvbuffers.get().at(1).get())->get(),
              (*(BufferedSharedArray<index_type_>*)((recvbuffers.get().at(0).get())))[1],
              sourcerank,
              s2);
        }


        virtual comm_link_type_& comm_link()
        {
          return _cl;
        }

        comm_link_type_& comm_link() const
        {
          return _cl;
        }

        virtual const compound_storage_type_* elements() const
        {
          return &_halo_elements;
        }

      protected:
        compound_storage_type_ _halo_elements;

        mesh_type_* _mesh;

        comm_link_type_ _cl;
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
             typename WT_ = double,
             template<typename, typename> class BufferStorageType_ = std::vector>
    class Halo :
      public HaloBase<MeshType_, WT_, BufferStorageType_>
    {
      public:
        ///type exports:
        typedef typename HaloBase<MeshType_, WT_, BufferStorageType_>::template CommLink<WT_> comm_link_type_;
        typedef typename HaloBase<MeshType_, WT_, BufferStorageType_>::index_type_ index_type_;
        typedef typename HaloBase<MeshType_, WT_, BufferStorageType_>::mesh_type_ mesh_type_;
        typedef typename HaloBase<MeshType_, WT_, BufferStorageType_>::buffer_type_ buffer_type_;
        typedef Level_ level_;

        ///CTOR
        Halo() :
          HaloBase<MeshType_, WT_, BufferStorageType_>(),
          _overlap(delta_),
          _level(Level_::tag_value)
        {
        }

        ///CTOR
        Halo(MeshType_ & mesh, index_type_ other = 0) : //TODO move to template
          HaloBase<MeshType_, WT_, BufferStorageType_>(mesh, other),
          _overlap(delta_),
          _level(Level_::tag_value)
        {
        }

        ///CTOR from buffer
        template<typename BufferType_>
        Halo(const BufferType_& buffer, MeshType_& new_mesh, index_type_ new_other_rank) :
          HaloBase<MeshType_, WT_, BufferStorageType_>(buffer, new_mesh, new_other_rank),
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
          HaloBase<MeshType_, WT_, BufferStorageType_>(other),
          _overlap(other._overlap),
          _level(other._level)
        {
        }

        virtual unsigned get_overlap() const override
        {
          return _overlap;
        }

        virtual PolytopeLevels get_level() const override
        {
          return _level;
        }

        Halo& operator=(const Halo& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_halo_elements = rhs._halo_elements;
          this->_mesh = rhs._mesh;
          this->_cl = rhs._cl;
          this->_overlap = rhs._overlap;
          this->_level = rhs._level;

          return *this;
        }

      private:
        unsigned _overlap;

        PolytopeLevels _level;
    };

    struct HaloTags
    {
      template<typename MeshType_,
               typename WT_ = double,
               template<typename, typename> class ST_ = std::vector>
      static ST_<Index, std::allocator<Index> > value(ST_<std::shared_ptr<HaloBase<MeshType_, WT_, ST_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_, WT_, ST_> > > >& halos)
      {
        ST_<Index, std::allocator<Index> > res;
        for(auto& h_i : halos)
          res.push_back(h_i->get_global_id());

        return res;
      }
    };


    template<typename MeshType_,
             typename WT_ = double,
             template<typename, typename> class StorageType_ = std::vector>
    bool compare_other(const std::shared_ptr<HaloBase<MeshType_, WT_, StorageType_> >& l, const std::shared_ptr<HaloBase<MeshType_, WT_, StorageType_> >& r)
    {
      return (l->get_other() < r->get_other());
    }

  }
}
#endif

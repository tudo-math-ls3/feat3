#pragma once
#ifndef KERNEL_FOUNDATION_HALO_HH
#define KERNEL_FOUNDATION_HALO_HH 1

#include<vector>
#include<kernel/foundation/mesh.hpp>

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
             typename MeshType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IndexType_ = Index>
    class Halo
    {
      public:
        ///type exports:
        typedef IndexType_ index_type_;
        typedef MeshType_ mesh_type_;

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

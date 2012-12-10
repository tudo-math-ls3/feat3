#pragma once
#ifndef KERNEL_SPACE_DOF_MAPPING_COMMON_HPP
#define KERNEL_SPACE_DOF_MAPPING_COMMON_HPP 1

// includes, FEAST
#include <kernel/space/dof_mapping_base.hpp>
#include <kernel/util/adjactor.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Identity-Dof-Mapping base-class template
     *
     * \note This class is meant to be used as a base-class only, therefore its constructor is protected.
     *
     * \tparam Space_
     * The finite-element space that this dof-mapping is used by.
     *
     * \tparam shape_dim_
     * The dimension of the shape that his dof-mapping is defined on.
     *
     * \tparam dofs_per_cell_
     * The number of dofs per entity.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      int dofs_per_cell_ = 1>
    class DofMappingIdentity :
      public DofMappingBase<Space_>
    {
    public:
      /// adjactor image iterator type
      typedef Adjactor::IndexImageIterator ImageIterator;

    protected:
      /// constructor
      explicit DofMappingIdentity(const Space_& space) :
        DofMappingBase<Space_>(space)
      {
      }

    public:
      /** \copydoc DofMappingBase::get_num_local_dofs() */
      Index get_num_local_dofs() const
      {
        return Index(dofs_per_cell_);
      }

      /** \copydoc DofMappingBase::get_index() */
      Index get_index(Index local_dof_idx, Index /*contrib_idx*/ = 0) const
      {
        return Index(dofs_per_cell_) * this->_cell_index + local_dof_idx;
      }

      /** \copydoc DofMappingBase::image_begin() */
      ImageIterator image_begin(Index domain_node) const
      {
        return ImageIterator(Index(dofs_per_cell_) * domain_node);
      }

      /** \copydoc DofMappingBase::image_end() */
      ImageIterator image_end(Index domain_node) const
      {
        return ImageIterator(Index(dofs_per_cell_) * (domain_node + 1));
      }
    }; // class DofMappingIdentity<...>

    /**
     * \brief Single-Entity Dof-Mapping class template
     *
     * This class implements the Dof-Mapping interface for a finite-element space which has dofs only
     * in one specific entity dimension.
     *
     * \tparam Space_
     * The finite-element space that this dof-mapping is used by.
     *
     * \tparam shape_dim_
     * The dimension of the shape that his dof-mapping is defined on.
     *
     * \tparam codim_
     * The co-dimension of the shape that the dofs are assigned to.
     *
     * \tparam dofs_per_cell_
     * The number of dofs per entity.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      int codim_,
      int dofs_per_cell_ = 1>
    class DofMappingSingleEntity :
      public DofMappingBase<Space_>
    {
    public:
      /// space typedef
      typedef Space_ SpaceType;
      /// shape type
      typedef typename SpaceType::ShapeType ShapeType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = ShapeType::dimension,
        /// dof dimension
        dof_dim = shape_dim - codim_
      };

      // make sure the dof-dimension is less than the shape dimension
      static_assert((codim_ >= 0) && (codim_ <= shape_dim), "invalid co-dimension");

    protected:
      /// trafo type
      typedef typename SpaceType::TrafoType TrafoType;
      /// mesh type
      typedef typename TrafoType::MeshType MeshType;
      /// index-set type
      typedef typename MeshType::template IndexSet<shape_dim, dof_dim>::Type IndexSetType;

    public:
      /**
       * \brief ImageIterator implementation
       */
      class ImageIterator
      {
      private:
        /// a pointer to the index-set
        const IndexSetType* _index_set;
        /// the current cell index
        Index _cell;
        /// the current local dof index
        Index _loc_dof;

      public:
        explicit ImageIterator(const IndexSetType* index_set, Index cell, Index loc_dof) :
          _index_set(index_set),
          _cell(cell),
          _loc_dof(loc_dof)
        {
        }

        ImageIterator& operator++()
        {
          ++_loc_dof;
          return *this;
        }

        Index operator*() const
        {
          Index ldi_q = _loc_dof / dofs_per_cell_;
          Index ldi_r = _loc_dof % dofs_per_cell_;
          return Index(dofs_per_cell_) * ((*_index_set)[_cell][ldi_q]) + ldi_r;
        }

        bool operator!=(const ImageIterator& other) const
        {
          return (_index_set != other._index_set) || (_cell != other._cell) || (_loc_dof != other._loc_dof);
        }
      }; // class DofMappingSingleEntity<...>::ImageIterator

    protected:
      /// dofs-at-cell index set reference
      const IndexSetType& _index_set;

    public:
      /** \copydoc DofMappingBase::DofMappingBase() */
      explicit DofMappingSingleEntity(const Space_& space) :
        DofMappingBase<Space_>(space),
        _index_set(space.get_trafo().get_mesh().template get_index_set<shape_dim, dof_dim>())
      {
      }

      /** \copydoc DofMappingBase::get_num_local_dofs() */
      Index get_num_local_dofs() const
      {
        return Index(IndexSetType::num_indices) * Index(dofs_per_cell_);
      }

      /** \copydoc DofMappingBase::get_index() */
      Index get_index(Index local_dof_idx, Index /*contrib_idx*/ = 0) const
      {
        Index ldi_q = local_dof_idx / dofs_per_cell_;
        Index ldi_r = local_dof_idx % dofs_per_cell_;
        return Index(dofs_per_cell_) * _index_set(this->_cell_index, ldi_q) + ldi_r;
        //return Index(dofs_per_cell_) * _index_set[this->_cell_index][ldi_q] + ldi_r;
      }

      /** \copydoc DofMappingBase::image_begin() */
      ImageIterator image_begin(Index domain_node) const
      {
        return ImageIterator(&_index_set, domain_node, 0);
      }

      /** \copydoc DofMappingBase::image_end() */
      ImageIterator image_end(Index domain_node) const
      {
        return ImageIterator(&_index_set, domain_node, get_num_local_dofs());
      }
    };

    /**
     * \brief Partial specialisation of DofMappingSingleEntity for co-dimension zero
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      int dofs_per_cell_>
    class DofMappingSingleEntity<Space_, 0, dofs_per_cell_> :
      public DofMappingIdentity<Space_, dofs_per_cell_>
    {
    public:
      /** \copydoc DofMappingBase::DofMappingBase() */
      explicit DofMappingSingleEntity(const Space_& space) :
        DofMappingIdentity<Space_, dofs_per_cell_>(space)
      {
      }
    };
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_MAPPING_COMMON_HPP

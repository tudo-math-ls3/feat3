// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_DOF_MAPPING_COMMON_HPP
#define KERNEL_SPACE_DOF_MAPPING_COMMON_HPP 1

// includes, FEAT
#include <kernel/space/dof_mapping_base.hpp>
#include <kernel/util/assertion.hpp>

namespace FEAT
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
    protected:
      /// constructor
      explicit DofMappingIdentity(const Space_& space) :
        DofMappingBase<Space_>(space)
      {
      }

    public:
      /** \copydoc DofMappingBase::get_num_local_dofs() */
      int get_num_local_dofs() const
      {
        return dofs_per_cell_;
      }

      /** \copydoc DofMappingBase::get_num_global_dofs() */
      Index get_num_global_dofs() const
      {
        return this->_space.get_trafo().get_mesh().get_num_entities(Space_::shape_dim) * Index(dofs_per_cell_);
      }

      /** \copydoc DofMappingBase::get_index() */
      Index get_index(int local_dof_idx) const
      {
        return Index(dofs_per_cell_) * this->_cell_index + Index(local_dof_idx);
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

      /// shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// dof dimension
      static constexpr int dof_dim = shape_dim - codim_;

      // make sure the dof-dimension is less than the shape dimension
      static_assert((codim_ >= 0) && (codim_ <= shape_dim), "invalid co-dimension");

    protected:
      /// trafo type
      typedef typename SpaceType::TrafoType TrafoType;
      /// mesh type
      typedef typename TrafoType::MeshType MeshType;
      /// index-set type
      typedef typename MeshType::template IndexSet<shape_dim, dof_dim>::Type IndexSetType;

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
      int get_num_local_dofs() const
      {
        return IndexSetType::num_indices * dofs_per_cell_;
      }

      /** \copydoc DofMappingBase::get_num_global_dofs() */
      Index get_num_global_dofs() const
      {
        return this->_space.get_trafo().get_mesh().get_num_entities(dof_dim) * Index(dofs_per_cell_);
      }

      /** \copydoc DofMappingBase::get_index() */
      Index get_index(int local_dof_idx) const
      {
        int ldi_q = local_dof_idx / dofs_per_cell_;
        int ldi_r = local_dof_idx % dofs_per_cell_;
        return Index(dofs_per_cell_) * _index_set(this->_cell_index, ldi_q) + Index(ldi_r);
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

    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        template<typename, int> class Traits_,
        typename Tag_,
        int cell_dim_ = Shape_::dimension,
        int shape_dim_ = Shape_::dimension>
      struct UniformDofMappingHelper;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Uniform Dof-Mapping class template
     *
     * This class implements the Dof-Mapping interface for a uniform finite-element space, which has
     * a fixed number of dofs for each entity dimension.
     *
     * \tparam Space_
     * The finite-element space that this dof-mapping is used by.
     *
     * \tparam Traits_
     * A traits class template defining the number of dofs per entity dimension.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      template<typename Tag_, int dim_> class DofTraits_,
      typename DofTag_>
    class DofMappingUniform
      : public DofMappingBase<Space_>
    {
    public:
      /// base-class typedef
      typedef DofMappingBase<Space_> BaseClass;
      /// space typedef
      typedef Space_ SpaceType;
      /// shape typedef
      typedef typename SpaceType::ShapeType ShapeType;

      // total number of local dofs
      static constexpr int dof_count = Intern::UniformDofMappingHelper<ShapeType, DofTraits_, DofTag_>::dof_count;

    private:
      /// dof index vector
      Index dof_idx[dof_count];

    public:
      explicit DofMappingUniform(const SpaceType& space) :
        BaseClass(space)
      {
      }

      /** \copydoc DofMappingBase::prepare() */
      void prepare(Index cell_index)
      {
        BaseClass::prepare(cell_index);
        Index count(0);
        Intern::UniformDofMappingHelper<ShapeType, DofTraits_, DofTag_>
          ::assemble(dof_idx, this->_space.get_mesh(), cell_index, count);
        XASSERTM(count == Index(dof_count), "dof-count mismatch");
      }

      /** \copydoc DofMappingBase::get_num_local_dofs() */
      int get_num_local_dofs() const
      {
        return dof_count;
      }

      /** \copydoc DofMappingBase::get_num_global_dofs() */
      Index get_num_global_dofs() const
      {
        return Intern::UniformDofMappingHelper<ShapeType, DofTraits_, DofTag_>::global_count(this->_space.get_mesh());
      }

      /** \copydoc DofMappingBase::get_index() */
      Index get_index(int local_dof_idx) const
      {
        ASSERTM((local_dof_idx >= 0) && (local_dof_idx < dof_count), "local dof-index out-of-range");
        return dof_idx[local_dof_idx];
      }
    };

    /// \cond internal
    namespace Intern
    {
      template<typename Shape_, template<typename, int> class Traits_, typename Tag_, int cell_dim_, int shape_dim_>
      struct UniformDofMappingHelper
      {
        static_assert(cell_dim_ < shape_dim_, "invalid cell dimension");

        static constexpr int cell_count = Shape::FaceTraits<Shape_, cell_dim_>::count;
        // number of dofs per cell
        static constexpr int dofs_per_cell = Traits_<Tag_, cell_dim_>::count;
        // number of dofs for all cells of this dimension
        static constexpr int dof_cell_count = cell_count * dofs_per_cell;
        // total dof count
        static constexpr int dof_count = UniformDofMappingHelper<Shape_, Traits_, Tag_, cell_dim_-1>::dof_count + dof_cell_count;

        template<typename MeshType_>
        static Index assemble(Index dof_idx[], const MeshType_& mesh, Index cell, Index& k)
        {
          Index offs(UniformDofMappingHelper<Shape_, Traits_, Tag_, cell_dim_-1>::assemble(dof_idx, mesh, cell, k));

          const auto& index_set(mesh.template get_index_set<shape_dim_, cell_dim_>());
          for(int i(0); i < cell_count; ++i)
          {
            for(int j(0); j < dofs_per_cell; ++j, ++k)
            {
              dof_idx[k] = offs + index_set(cell,i) * Index(dofs_per_cell) + Index(j);
            }
          }
          return offs + Index(dofs_per_cell) * mesh.get_num_entities(cell_dim_);
        }

        template<typename MeshType_>
        static Index global_count(const MeshType_& mesh)
        {
          return UniformDofMappingHelper<Shape_, Traits_, Tag_, cell_dim_-1>::global_count(mesh) +
            Index(Traits_<Tag_, cell_dim_>::count) * mesh.get_num_entities(cell_dim_);
        }
      };

      // specialisation for cell_dim_ = shape_dim_
      template<typename Shape_, template<typename, int> class Traits_, typename Tag_, int cell_dim_>
      struct UniformDofMappingHelper<Shape_, Traits_, Tag_, cell_dim_, cell_dim_>
      {
        // number of dofs per cell
        static constexpr int dofs_per_cell = Traits_<Tag_, cell_dim_>::count;
        // number of dofs for all cells of this dimension
        static constexpr int dof_cell_count = dofs_per_cell;
        // total dof count
        static constexpr int dof_count = UniformDofMappingHelper<Shape_, Traits_, Tag_, cell_dim_-1>::dof_count + dof_cell_count;

        template<typename MeshType_>
        static Index assemble(Index dof_idx[], const MeshType_& mesh, Index cell, Index& k)
        {
          Index offs(UniformDofMappingHelper<Shape_, Traits_, Tag_, cell_dim_-1>::assemble(dof_idx, mesh, cell, k));

          for(int j(0); j < dofs_per_cell; ++j, ++k)
          {
            dof_idx[k] = offs + cell * dofs_per_cell + Index(j);
          }
          return offs + Index(dofs_per_cell) * mesh.get_num_entities(cell_dim_);
        }

        template<typename MeshType_>
        static Index global_count(const MeshType_& mesh)
        {
          return UniformDofMappingHelper<Shape_, Traits_, Tag_, cell_dim_-1>::global_count(mesh) +
            Index(Traits_<Tag_, cell_dim_>::count) * mesh.get_num_entities(cell_dim_);
        }
      };

      // specialisation for cell_dim_ = 0
      template<typename Shape_, template<typename, int> class Traits_, typename Tag_, int shape_dim_>
      struct UniformDofMappingHelper<Shape_, Traits_, Tag_, 0, shape_dim_>
      {
        static constexpr int cell_count = Shape::FaceTraits<Shape_, 0>::count;
        // number of dofs per cell
        static constexpr int dofs_per_cell = Traits_<Tag_, 0>::count;
        // number of dofs for all cells of this dimension
        static constexpr int dof_cell_count = cell_count * dofs_per_cell;
        // total dof count
        static constexpr int dof_count = dof_cell_count;

        template<typename MeshType_>
        static Index assemble(Index dof_idx[], const MeshType_& mesh, Index cell, Index& k)
        {
          const auto& index_set(mesh.template get_index_set<shape_dim_, 0>());

          k = 0;
          for(int i(0); i < cell_count; ++i)
          {
            for(int j(0); j < dofs_per_cell; ++j, ++k)
            {
              dof_idx[k] = index_set(cell,i) * dofs_per_cell + Index(j);
            }
          }
          return Index(dofs_per_cell) * mesh.get_num_entities(0);
        }

        template<typename MeshType_>
        static Index global_count(const MeshType_& mesh)
        {
          return Index(Traits_<Tag_, 0>::count) * mesh.get_num_entities(0);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DOF_MAPPING_COMMON_HPP

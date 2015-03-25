#pragma once
#ifndef KERNEL_SPACE_DOF_ASSIGNMENT_BASE_HPP
#define KERNEL_SPACE_DOF_ASSIGNMENT_BASE_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Finite-Element Dof-Assignment base-class template.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      int shape_dim_,
      typename DataType_>
    class DofAssignmentBase
    {
    public:
      /// space typedef
      typedef Space_ SpaceType;

    protected:
      /// space reference
      const SpaceType& _space;
      /// currently active cell index
      Index _cell_index;

      /// constructor
      explicit DofAssignmentBase(const SpaceType& space) :
        _space(space),
        _cell_index(~Index(0))
      {
      }

    public:
      void prepare(Index cell_index)
      {
        _cell_index = cell_index;
      }

      void finish()
      {
        _cell_index = ~Index(0);
      }

      Index get_max_assigned_dofs() const
      {
        return 1;
      }

#ifdef DOXYGEN
      Index get_num_assigned_dofs() const;
#endif // DOXYGEN

      Index get_max_contribs() const
      {
        return 1;
      }

      Index get_num_contribs(Index /*assign_idx*/) const
      {
        return 1;
      }

      Index get_derive_order(Index /*assign_idx*/) const
      {
        return 0;
      }

#ifdef DOXYGEN
      Index get_index(Index assign_idx, Index contrib_idx = 0) const;
#endif // DOXYGEN

      DataType_ get_weight(Index /*assign_idx*/, Index /*contrib_idx*/ = 0) const
      {
        return DataType_(1.0);
      }
    }; // class DofAssignmentBase

    /// \cond internal
    namespace Intern
    {
      template<typename Tag_, template<typename, int> class Traits_, int cell_dim_>
      struct UniformDofAssignHelper
      {
        template<typename Mesh_>
        static Index dof_offset(const Mesh_& mesh)
        {
          return UniformDofAssignHelper<Tag_, Traits_, cell_dim_ - 1>::dof_offset(mesh)
            + mesh.get_num_entities(cell_dim_-1) * Index(Traits_<Tag_, cell_dim_ - 1>::count);
        }
      };

      template<typename Tag_, template<typename, int> class Traits_>
      struct UniformDofAssignHelper<Tag_, Traits_, 0>
      {
        template<typename Mesh_>
        static Index dof_offset(const Mesh_&)
        {
          return Index(0);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Uniform Dof-Assignment class template.
     *
     * \tparam Space_
     * The finite element space this dof-assignment is used by.
     *
     * \tparam shape_dim_
     * The dimension of the shape that this dof-assignment refers to.
     *
     * \tparam DataType_
     * The data-type to be used by the dof-assignment.
     *
     * \tparam DofTraits_
     * A dof-traits class template that defines the number of dofs per dimension.
     *
     * \tparam DofTag_
     * A tag class that is passed as a first parameter to the DofTraits_ class template.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      int shape_dim_,
      typename DataType_,
      template<typename Tag_, int dim_> class DofTraits_,
      typename DofTag_>
    class DofAssignmentUniform
      : public DofAssignmentBase<Space_, shape_dim_, DataType_>
    {
    public:
      typedef DofAssignmentBase<Space_, shape_dim_, DataType_> BaseClass;
      typedef Space_ SpaceType;
      typedef typename SpaceType::ShapeType ShapeType;

      // number of dofs for each cell
      static constexpr int dofs_per_cell = DofTraits_<DofTag_, shape_dim_>::count;

    protected:
      /// the offset of the first dof of this cell dimension
      Index _dof_offset;

    public:
      explicit DofAssignmentUniform(const SpaceType& space) :
        BaseClass(space),
        _dof_offset(Intern::UniformDofAssignHelper<DofTag_, DofTraits_, shape_dim_>::dof_offset(space.get_mesh()))
      {
      }

      /** \copydoc DofAssignmentBase::get_max_assigned_dofs() */
      Index get_max_assigned_dofs() const
      {
        return Index(dofs_per_cell);
      }

      /** \copydoc DofAssignmentBase::get_num_assigned_dofs() */
      Index get_num_assigned_dofs() const
      {
        return Index(dofs_per_cell);
      }

      /** \copydoc DofAssignmentBase::get_derive_order() */
      Index get_derive_order(Index assign_idx)
      {
        return DofTraits_<DofTag_, shape_dim_>::derive_order(assign_idx);
      }

      /** \copydoc DofAssignmentBase::get_index() */
      Index get_index(Index assign_idx, Index DOXY(contrib_idx) = 0) const
      {
        return _dof_offset + Index(dofs_per_cell) * this->_cell_index + assign_idx;
      }
    };
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_ASSIGNMENT_BASE_HPP

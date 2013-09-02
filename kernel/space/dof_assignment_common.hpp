#pragma once
#ifndef KERNEL_SPACE_DOF_ASSIGNMENT_COMMON_HPP
#define KERNEL_SPACE_DOF_ASSIGNMENT_COMMON_HPP 1

// includes, FEAST
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/adjacency/adjactor.hpp>

namespace FEAST
{
  namespace Space
  {
    template<
      typename Space_,
      int shape_dim_,
      typename DataType_>
    class DofAssignmentNull :
      public DofAssignmentBase<Space_, shape_dim_, DataType_>
    {
      typedef DofAssignmentBase<Space_, shape_dim_, DataType_> BaseClass;

    public:
      explicit DofAssignmentNull(const Space_& space) :
        BaseClass(space)
      {
      }

      /** \copydoc DofAssignmentBase::get_max_assigned_dofs() */
      Index get_max_assigned_dofs() const
      {
        return 0;
      }

      /** \copydoc DofAssignmentBase::get_num_assigned_dofs() */
      Index get_num_assigned_dofs() const
      {
        return 0;
      }


      Index get_max_contribs() const
      {
        return 0;
      }

      Index get_num_contribs(Index /*assign_idx*/) const
      {
        throw InternalError("invalid call of DofAssignmentNull::get_num_contribs()");
      }

      /** \copydoc DofAssignmentBase::get_index() */
      Index get_index(Index /*assign_idx*/, Index /*contrib_idx*/ = 0) const
      {
        throw InternalError("invalid call of DofAssignmentNull::get_index()");
      }

      /** \copydoc DofAssignmentBase::get_weight() */
      DataType_ get_weight(Index /*assign_idx*/, Index /*contrib_idx*/ = 0) const
      {
        throw InternalError("invalid call of DofAssignmentNull::get_weight()");
      }
    };

    template<
      typename Space_,
      int shape_dim_,
      typename DataType_,
      int dofs_per_cell_ = 1>
    class DofAssignmentIdentity :
      public DofAssignmentBase<Space_, shape_dim_, DataType_>
    {
      typedef DofAssignmentBase<Space_, shape_dim_, DataType_> BaseClass;


    public:
      explicit DofAssignmentIdentity(const Space_& space) :
        BaseClass(space)
      {
      }

      /** \copydoc DofAssignmentBase::get_max_assigned_dofs() */
      Index get_max_assigned_dofs() const
      {
        return Index(dofs_per_cell_);
      }

      /** \copydoc DofAssignmentBase::get_num_assigned_dofs() */
      Index get_num_assigned_dofs() const
      {
        return Index(dofs_per_cell_);
      }

      /** \copydoc DofAssignmentBase::get_index() */
      Index get_index(Index assign_idx, Index DOXY(contrib_idx) = 0) const
      {
        return Index(dofs_per_cell_) * this->_cell_index + assign_idx;
      }
    }; // class DofAssignmentIdentity

    template<
      typename Space_,
      int shape_dim_,
      typename DataType_,
      int dof_dim_,
      int dofs_per_cell_ = 1>
    class DofAssignmentSingleEntity :
      public DofAssignmentNull<Space_, shape_dim_, DataType_>
    {
      typedef DofAssignmentNull<Space_, shape_dim_, DataType_> BaseClass;

    public:
      explicit DofAssignmentSingleEntity(const Space_& space) :
        BaseClass(space)
      {
      }
    };

    template<
      typename Space_,
      int shape_dim_,
      typename DataType_,
      int dofs_per_cell_>
    class DofAssignmentSingleEntity<Space_, shape_dim_, DataType_, shape_dim_, dofs_per_cell_> :
      public DofAssignmentIdentity<Space_, shape_dim_, DataType_>
    {
      typedef DofAssignmentIdentity<Space_, shape_dim_, DataType_> BaseClass;

    public:
      explicit DofAssignmentSingleEntity(const Space_& space) :
        BaseClass(space)
      {
      }
    };
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_ASSIGNMENT_COMMON_HPP

#pragma once
#ifndef KERNEL_SPACE_DOF_ASSIGNMENT_COMMON_HPP
#define KERNEL_SPACE_DOF_ASSIGNMENT_COMMON_HPP 1

// includes, FEAT
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/adjacency/adjactor.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Null Dof-Assignment class template
     *
     * This class defines a dof-assignment which has no dofs for the selected shape dimension.
     *
     * \tparam Space_
     * The finite element space that this dof-assignment is used by.
     *
     * \tparam shape_dim_
     * The dimension of the shape whose dof-assignment is to be considered.
     *
     * \tparam DataType_
     * The data type used by the dof-assignment.
     *
     * \author Peter Zajac
     */
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
      int get_max_assigned_dofs() const
      {
        return 0;
      }

      /** \copydoc DofAssignmentBase::get_num_assigned_dofs() */
      int get_num_assigned_dofs() const
      {
        return 0;
      }

      /** \copydoc DofAssignmentBase::get_index() */
      Index get_index(int /*assign_idx*/) const
      {
        throw InternalError("invalid call of DofAssignmentNull::get_index()");
      }
    }; // class DofAssignmentNull<...>

    /**
     * \brief Identity Dof-Assignment class template
     *
     * This class defines a dof-assignment which has no dofs for the selected shape dimension.
     *
     * \tparam Space_
     * The finite element space that this dof-assignment is used by.
     *
     * \tparam shape_dim_
     * The dimension of the shape whose dof-assignment is to be considered.
     *
     * \tparam DataType_
     * The data type used by the dof-assignment.
     *
     * \author Peter Zajac
     */
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
      int get_max_assigned_dofs() const
      {
        return dofs_per_cell_;
      }

      /** \copydoc DofAssignmentBase::get_num_assigned_dofs() */
      int get_num_assigned_dofs() const
      {
        return dofs_per_cell_;
      }

      /** \copydoc DofAssignmentBase::get_index() */
      Index get_index(int assign_idx) const
      {
        return Index(dofs_per_cell_) * this->_cell_index + Index(assign_idx);
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
} // namespace FEAT

#endif // KERNEL_SPACE_DOF_ASSIGNMENT_COMMON_HPP

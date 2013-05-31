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

#ifdef DOXYGEN
      Index get_index(Index assign_idx, Index contrib_idx = 0) const;
#endif // DOXYGEN

      DataType_ get_weight(Index assign_idx, Index contrib_idx = 0) const
      {
        return DataType_(1.0);
      }
    }; // class DofAssignmentBase
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_ASSIGNMENT_BASE_HPP

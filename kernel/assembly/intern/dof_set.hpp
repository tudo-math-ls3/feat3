#pragma once
#ifndef KERNEL_ASSEMBLY_INTERN_DOF_SET_HPP
#define KERNEL_ASSEMBLY_INTERN_DOF_SET_HPP 1

// includes, FEAST
#include <kernel/assembly/base.hpp>

// includes, system
#include <set>

namespace FEAST
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Space_,
        typename CellSet_,
        int shape_dim_>
      struct DofSetAdder
      {
        static void apply(std::set<Index>& idx, const Space_& space, const CellSet_& cell_set)
        {
          // fetch the target set for this dimension
          const typename CellSet_::template TargetSet<shape_dim_>::Type&
            target_set(cell_set.template get_target_set<shape_dim_>());

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);

          // loop over all target indices
          const Index num_entities = target_set.get_num_entities();
          for(Index i(0); i < num_entities; ++i)
          {
            dof_assign.prepare(target_set[i]);
            const Index num_assign(dof_assign.get_num_assigned_dofs());
            for(Index j(0); j < num_assign; ++j)
            {
              const Index num_contribs(dof_assign.get_num_contribs(j));
              for(Index k(0); k < num_contribs; ++k)
              {
                idx.insert(dof_assign.get_index(j, k));
              }
            }
            dof_assign.finish();
          }
        }
      };

      template<
        typename Space_,
        typename CellSet_,
        int shape_dim_ = Space_::shape_dim>
      struct DofSetAddWrapper
      {
        static void apply(std::set<Index>& idx, const Space_& space, const CellSet_& cell_set)
        {
          DofSetAddWrapper<Space_, CellSet_, shape_dim_ - 1>::apply(idx, space, cell_set);
          DofSetAdder<Space_, CellSet_, shape_dim_>::apply(idx, space, cell_set);
        }
      };

      template<
        typename Space_,
        typename CellSet_>
      struct DofSetAddWrapper<Space_, CellSet_, 0>
      {
        static void apply(std::set<Index>& idx, const Space_& space, const CellSet_& cell_set)
        {
          DofSetAdder<Space_, CellSet_, 0>::apply(idx, space, cell_set);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_INTERN_DOF_SET_HPP

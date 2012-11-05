#pragma once
#ifndef KERNEL_FINITE_ELEMENT_DOF_SUBSET_HPP
#define KERNEL_FINITE_ELEMENT_DOF_SUBSET_HPP 1

// includes, FEAST
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/shape.hpp>

namespace FEAST
{
  namespace Space
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Space_,
        typename CellSet_,
        int shape_dim_>
      struct DofSubsetHelper
      {
        static Index count(Index& contribs, const Space_& space, const CellSet_& cell_set)
        {
          // fetch the target set for this dimension
          const typename CellSet_::template TargetSet<shape_dim_>::Type&
            target_set(cell_set.template get_target_set<shape_dim_>());
          if(target_set.get_num_entities() <= 0)
            return 0;

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);
          if(dof_assign.get_max_assigned_dofs() <= 0)
            return 0;

          // loop over all target indices
          Index num_entities = target_set.get_num_entities();
          Index count(0);
          for(Index i(0); i < num_entities; ++i)
          {
            dof_assign.prepare(target_set[i]);
            Index num_assign(dof_assign.get_num_assigned_dofs());
            for(Index j(0); j < num_assign; ++j)
            {
              contribs += dof_assign.get_num_contribs(j);
            }
            count += num_assign;
            dof_assign.finish();
          }

          return count;
        }

        static Index fill(Index* ptr, Index idx[], Index offset, const Space_& space, const CellSet_& cell_set)
        {
          // fetch the target set for this dimension
          const typename CellSet_::template TargetSet<shape_dim_>::Type&
            target_set(cell_set.template get_target_set<shape_dim_>());
          if(target_set.get_num_entities() <= 0)
            return offset;

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);
          if(dof_assign.get_max_assigned_dofs() <= 0)
            return offset;

          // loop over all target indices
          Index num_entities = target_set.get_num_entities();
          for(Index i(0); i < num_entities; ++i)
          {
            dof_assign.prepare(target_set[i]);
            Index num_assign(dof_assign.get_num_assigned_dofs());
            for(Index j(0); j < num_assign; ++j, ++ptr)
            {
              *ptr = offset;
              Index num_contribs(dof_assign.get_num_contribs(j));
              for(Index k(0); k < num_contribs; ++k, ++offset)
              {
                idx[offset] = dof_assign.get_index(j, k);
              }
            }
            *ptr = offset;
            dof_assign.finish();
          }

          return offset;
        }
      };

      template<
        typename Space_,
        typename CellSet_,
        int shape_dim_ = Space_::shape_dim>
      struct DofSubsetHelpWrapper
      {
        static Index count(Index& contribs, const Space_& space, const CellSet_& cell_set)
        {
          // recursive call
          return DofSubsetHelpWrapper<Space_, CellSet_, shape_dim_ - 1>::count(contribs, space, cell_set) +
            DofSubsetHelper<Space_, CellSet_, shape_dim_>::count(contribs, space, cell_set);
        }

        static Index fill(Index* ptr, Index idx[], const Space_& space, const CellSet_& cell_set)
        {
          Index offset =  DofSubsetHelpWrapper<Space_, CellSet_, shape_dim_ - 1>::fill(ptr, idx, space, cell_set);
          return DofSubsetHelper<Space_, CellSet_, shape_dim_>::fill(ptr, idx, offset, space, cell_set);
        }
      };

      template<
        typename Space_,
        typename CellSet_>
      struct DofSubsetHelpWrapper<Space_, CellSet_, 0>
      {
        static Index count(Index& contribs, const Space_& space, const CellSet_& cell_set)
        {
          return DofSubsetHelper<Space_, CellSet_, 0>::count(contribs, space, cell_set);
        }

        static Index fill(Index* ptr, Index idx[], const Space_& space, const CellSet_& cell_set)
        {
          return DofSubsetHelper<Space_, CellSet_, 0>::fill(ptr, idx, 0, space, cell_set);
        }
      };
    } // namespace Intern
    /// \endcond

    class DofSubset
    {
    public:
      template<
        typename Space_,
        typename CellSet_>
      static Graph assemble(
        const Space_& space,
        const CellSet_& cell_set)
      {
        Index contribs = 0;
        Index count = Intern::DofSubsetHelpWrapper<Space_, CellSet_>::count(contribs, space, cell_set);
        Graph graph(count, space.get_num_dofs(), contribs);
        Index* ptr = graph.get_domain_ptr();
        Index* idx = graph.get_image_idx();
        Intern::DofSubsetHelpWrapper<Space_, CellSet_>::fill(ptr, idx, space, cell_set);
        return graph;
      }
    };
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_FINITE_ELEMENT_DOF_SUBSET_HPP

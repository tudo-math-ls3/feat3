#pragma once
#ifndef KERNEL_SPACE_DOF_MIRROR_HPP
#define KERNEL_SPACE_DOF_MIRROR_HPP 1

// includes, FEAST
#include <kernel/adjacency/graph.hpp>

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
      struct DofMirrorHelper
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
      struct DofMirrorHelpWrapper
      {
        static Index count(Index& contribs, const Space_& space, const CellSet_& cell_set)
        {
          // recursive call
          return DofMirrorHelpWrapper<Space_, CellSet_, shape_dim_ - 1>::count(contribs, space, cell_set) +
            DofMirrorHelper<Space_, CellSet_, shape_dim_>::count(contribs, space, cell_set);
        }

        static Index fill(Index* ptr, Index idx[], const Space_& space, const CellSet_& cell_set)
        {
          Index offset =  DofMirrorHelpWrapper<Space_, CellSet_, shape_dim_ - 1>::fill(ptr, idx, space, cell_set);
          return DofMirrorHelper<Space_, CellSet_, shape_dim_>::fill(ptr, idx, offset, space, cell_set);
        }
      };

      template<
        typename Space_,
        typename CellSet_>
      struct DofMirrorHelpWrapper<Space_, CellSet_, 0>
      {
        static Index count(Index& contribs, const Space_& space, const CellSet_& cell_set)
        {
          return DofMirrorHelper<Space_, CellSet_, 0>::count(contribs, space, cell_set);
        }

        static Index fill(Index* ptr, Index idx[], const Space_& space, const CellSet_& cell_set)
        {
          return DofMirrorHelper<Space_, CellSet_, 0>::fill(ptr, idx, 0, space, cell_set);
        }
      };

      template<
        typename Space_,
        int shape_dim_ = Space_::shape_dim>
      struct MaxDofContrib
      {
        static Index value(const Space_& space)
        {
          Index max_contrib(MaxDofContrib<Space_, shape_dim_-1>::value(space));
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);
          return std::max(max_contrib, dof_assign.get_max_contribs());
        }
      };

      template<typename Space_>
      struct MaxDofContrib<Space_, 0>
      {
        static Index value(const Space_& space)
        {
          typename Space_::template DofAssignment<0>::Type dof_assign(space);
          return dof_assign.get_max_contribs();
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Dof-Mirror assembler class template.
     *
     * \author Peter Zajac
     */
    class DofMirror
    {
    public:
      /**
       * \brief Assembles the Dof-Mirror adjacency graph.
       */
      template<
        typename Space_,
        typename CellSet_>
      static Adjacency::Graph assemble(
        const Space_& space,
        const CellSet_& cell_set)
      {
        // ensure that the space has not more than 1 Dof contribution...
        if(Intern::MaxDofContrib<Space_>::value(space) != 1)
          throw InternalError("Cannot compute Dof-Mirror graph: multiple DOF contributions!");

        Index contribs = 0;
        Index count = Intern::DofMirrorHelpWrapper<Space_, CellSet_>::count(contribs, space, cell_set);
        Adjacency::Graph graph(count, space.get_num_dofs(), contribs);
        Index* ptr = graph.get_domain_ptr();
        Index* idx = graph.get_image_idx();
        Intern::DofMirrorHelpWrapper<Space_, CellSet_>::fill(ptr, idx, space, cell_set);
        return graph;
      }
    }; // class DofMirror

  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_MIRROR_HPP

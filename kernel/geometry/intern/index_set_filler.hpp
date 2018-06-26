#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_INDEX_SET_FILLER_HPP
#define KERNEL_GEOMETRY_INTERN_INDEX_SET_FILLER_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Fills a MeshPart's IndexSetHolder by using TargetSet information and the parent's IndexSetHolder
       *
       * \tparam dim
       * Shape dimension of the IndexSet that gets filled.
       *
       * This only fills the vertex@shape information and leaves the rest to the RedundantIndexSetBuilder
       *
       * \author Jordi Paul
       */
      template<int dim>
      struct IndexSetFiller
      {
        template<typename IndexSetHolderType, typename TargetSetHolderType, typename ParentIndexSetHolderType>
        static void fill_ish(IndexSetHolderType& ish, const TargetSetHolderType& tsh, const ParentIndexSetHolderType& parent_ish)
        {
          // Recurse down
          IndexSetFiller<dim-1>::fill_ish(ish, tsh, parent_ish);

          // The MeshPart's vertex@shape set
          auto& index_set(ish.template get_index_set<dim, 0>());
          // The Mesh's vertex@shape set
          auto& index_set_parent(parent_ish.template get_index_set<dim, 0>());

          auto& target_set_vertex(tsh.template get_target_set<0>());
          auto& target_set_dim(tsh.template get_target_set<dim>());

          // For every vertex in the parent, this will contain its index in the MeshPart
          std::vector<Index> inverse_target_map(index_set_parent.get_index_bound());

          // Set every index to something out of range to catch errors
          for(Index i(0); i < index_set_parent.get_index_bound(); ++i)
            inverse_target_map[i] = target_set_vertex.get_num_entities() + Index(1);

          for(Index i(0); i < target_set_vertex.get_num_entities(); ++i)
            inverse_target_map[target_set_vertex[i]] = i;

          // Now we can just iterate over the shapes of the MeshPart
          for(Index cell(0); cell < target_set_dim.get_num_entities(); ++cell)
          {
            // For the shape cell, get its index in the parent. Then get the local vertex' index in the parent and
            // map that back to the MeshPart with the inverse_target_map. Ez!
            for(int i(0); i < index_set_parent.get_num_indices(); ++i)
              index_set[cell][i] = inverse_target_map[index_set_parent[target_set_dim[cell]][i]];
          }
        }
      }; // struct IndexSetFiller<int>

      /**
       * \brief Explicit specialisation as end of template recursion.
       */
      template<>
      struct IndexSetFiller<0>
      {
        template<typename IndexSetHolderType, typename TargetSetHolderType, typename ParentIndexSetHolderType>
        static void fill_ish(IndexSetHolderType& DOXY(ish), const TargetSetHolderType& DOXY(tsh), const ParentIndexSetHolderType& DOXY(parent_ish))
        {
        }
      }; // struct IndexSetFiller<0>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_INDEX_SET_FILLER_HPP

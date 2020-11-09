// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_TARGET_SET_COMPUTER_HPP
#define KERNEL_GEOMETRY_INTERN_TARGET_SET_COMPUTER_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/target_set.hpp>

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
       * \brief Wrapper class for filling TargetSetHolder objects
       *
       * \tparam end_dim_
       * Dimension to stop the template recursion at.
       *
       * \tparam current_dim_
       * Dimension for which a new container is filled.
       *
       * \author Jordi Paul
       *
       */
      template<int end_dim_, int current_dim_>
      struct TargetSetComputer
      {
        /**
         * \brief Fills a TargetSetHolder from bottom to top using information from an IndexSetHolder
         *
         * The IndexSetHolder contains the topology of the parent mesh that the TargetSetHolder refers to.
         *
         * \tparam TargetSetHolderType
         * Type of the TargetSetHolder to be filled.
         *
         * \tparam IndexSetHolderType
         * Type containing mesh topology information.
         *
         */
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void bottom_to_top(TargetSetHolderType& _target_set_holder, const IndexSetHolderType& parent_ish)
        {
          // Template recursion: Call the lower dimensional version first
          TargetSetComputer<end_dim_, current_dim_-1>::template bottom_to_top(_target_set_holder, parent_ish);

          typedef typename IndexSetHolderType::template IndexSet<current_dim_, current_dim_-1>::Type ParentIndexSetType;

          // TargetSet to fill
          TargetSet& my_target_set(_target_set_holder.template get_target_set<current_dim_>());

          // Only do something if the target set is still empty
          if(my_target_set.get_num_entities() == 0)
          {
            // Lower dimensional target set known to exist
            TargetSet& ts_below(_target_set_holder.template get_target_set<current_dim_-1>());

            // Indexset for entities of dimension current_dim_-1 at entities of dimension current_dim_
            const ParentIndexSetType& is_parent_below(parent_ish.template get_index_set<current_dim_, current_dim_-1>());
            // IndexSet for storing the indices of the MeshPart entities lying at entities of the parent
            ParentIndexSetType lower_parent_to_mesh_part(Index(is_parent_below.get_num_indices()));

            // For every entity of the parent, this will save whether it is referenced by the TargetSet
            TargetSet lower_origin_set(is_parent_below.get_index_bound());
            // and this is the marker for that
            Index marker(is_parent_below.get_index_bound());
            for(Index i(0); i < lower_origin_set.get_num_entities(); ++i)
              lower_origin_set[i] = 0;
            // Set the values to something bogus for every index present
            for(Index i(0); i < ts_below.get_num_entities(); ++i)
              lower_origin_set[ts_below[i]] = marker;

            // Count the number of entities that get added
            Index num_entities_current_dim(0);

            // Temporary TargetSet initialized with the maximum possible size (= number of entities in the parent mesh,
            // which is the case if the MeshPart contains all entities of the parent mesh)
            TargetSet tmp_target_set(is_parent_below.get_num_entities());

            // Check all entities of dimension current_dim_ from the parent
            for(Index i(0); i < is_parent_below.get_num_entities(); ++i)
            {
              bool is_in_mesh_part(true);

              // Check if all entities of dimension current_dim_-1 at entity i are referenced by the TargetSet
              for(int j(0); j < ParentIndexSetType::num_indices; ++j)
              {
                // This is the index in the parent
                Index parent_index(is_parent_below[i][j]);
                if(lower_origin_set[parent_index] != marker)
                  is_in_mesh_part = false;
              }

              // If all subshapes belonged to the MeshPart, create the new parent mapping
              if(is_in_mesh_part)
              {
                tmp_target_set[num_entities_current_dim] = i;
                num_entities_current_dim++;
              }
            }

            // tmp_target_set possibly has the wrong size, so manually create a correctly sized TargetSet
            TargetSet new_target_set(num_entities_current_dim);
            for(Index i(0); i < new_target_set.get_num_entities(); ++i)
              new_target_set[i] = tmp_target_set[i];

            // Update the information in the TargetSetHolder
            my_target_set = std::move(new_target_set);
          }
        } // void bottom_to_top<typename, typename>()

        /**
         * \brief Fills a TargetSetHolder from top to bottom using information from an IndexSetHolder
         *
         * The IndexSetHolder contains the topology of the parent mesh that the TargetSetHolder refers to.
         *
         * \tparam TargetSetHolderType
         * Type of the TargetSetHolder to be filled.
         *
         * \tparam IndexSetHolderType
         * Type containing mesh topology information.
         *
         */
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void top_to_bottom(TargetSetHolderType& _target_set_holder, const IndexSetHolderType& parent_ish)
        {
          // Call higher dimensional version first
          TargetSetComputer<end_dim_, current_dim_+1>::template top_to_bottom(_target_set_holder, parent_ish);

          typedef typename IndexSetHolderType::template IndexSet<current_dim_+1, current_dim_>::Type ParentIndexSetType;

          // TargetSet to fill
          TargetSet& my_target_set(_target_set_holder.template get_target_set<current_dim_>());

          // Only do something if the target set is still empty
          if(my_target_set.get_num_entities() == 0)
          {
            // Higher dimensional target set known to exist
            TargetSet& ts_above(_target_set_holder.template get_target_set<current_dim_+1>());

            // Indexset for entities of dimension current_dim_ at entities of dimension current_dim_+1
            const ParentIndexSetType& is_parent_above(parent_ish.template get_index_set<current_dim_+1, current_dim_>());
            // IndexSet for storing the indices of the MeshPart entities lying at entities of the parent
            ParentIndexSetType upper_parent_to_mesh_part(Index(is_parent_above.get_num_indices()));

            // For every entity of current_dim_ of the parent, this will save whether it belongs to an entity of
            // dimensin current_dim_+1 that is referenced through the TargetSetHolder
            TargetSet upper_origin_set(is_parent_above.get_index_bound());
            // And this is the marker for that
            Index marker(is_parent_above.get_index_bound());
            for(Index i(0); i < upper_origin_set.get_num_entities(); ++i)
              upper_origin_set[i] = 0;

            for(Index i(0); i < ts_above.get_num_entities(); ++i)
            {
              // Index of the current entity in the parent
              Index parent_index(ts_above[i]);
              for(int j(0); j < is_parent_above.get_num_indices(); ++j)
              {
                // This is i.e. the global number of local edge j at face parent_index in the parent
                Index parent_sub_entity_index(is_parent_above(parent_index, j));
                upper_origin_set[parent_sub_entity_index] = marker;
              }
            }

            // Temporary TargetSet initialized with the maximum possible size (= number of entities in the parent mesh,
            // which is the case if the MeshPart contains all entities of the parent mesh)
            TargetSet tmp_target_set(is_parent_above.get_index_bound());
            // Count the number of entities that get added
            Index num_entities_current_dim(0);
            for(Index i(0); i < upper_origin_set.get_num_entities(); ++i)
            {
              if(upper_origin_set[i] == marker)
              {
                tmp_target_set[num_entities_current_dim] = i;
                num_entities_current_dim++;
              }
            }

            // tmp_target_set possibly has the wrong size, so manually create a correctly sized TargetSet
            TargetSet new_target_set(num_entities_current_dim);
            for(Index i(0); i < new_target_set.get_num_entities(); ++i)
              new_target_set[i] = tmp_target_set[i];

            // Update the information in the TargetSetHolder
            my_target_set = std::move(new_target_set);
          }
        }
      }; // struct TargetSetComputer<Index, Index>

      /**
       * \brief Specialization as end of template recursion
       */
      template<int end_dim_>
      struct TargetSetComputer<end_dim_, end_dim_>
      {
        /// \brief End of template recursion: Nothing to do
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void bottom_to_top(TargetSetHolderType& , const IndexSetHolderType&)
        {
        }

        /// \brief End of template recursion: Nothing to do
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void top_to_bottom(TargetSetHolderType& , const IndexSetHolderType&)
        {
        }
      }; // struct TargetSetComputer<int>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_TARGET_SET_COMPUTER_HPP

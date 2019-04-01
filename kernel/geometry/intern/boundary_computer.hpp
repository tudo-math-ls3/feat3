// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_BOUNDARY_COMPUTER_HPP
#define KERNEL_GEOMETRY_INTERN_BOUNDARY_COMPUTER_HPP 1

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
      template<
        typename Shape_,
        // Note: The shape_dim_ parameter *must* always be equal to Shape_::dimension !
        int shape_dim_ = Shape_::dimension,
        int cell_dim_ = shape_dim_>
      class BoundaryFaceComputer :
        public BoundaryFaceComputer<Shape_, shape_dim_, cell_dim_-1>
      {
      public:
        typedef BoundaryFaceComputer<Shape_, shape_dim_, cell_dim_-1> BaseClass;

      protected:
        /// a vector of all (cell_dim_-1)-dimensional boundary face indices
        std::vector<Index> _faces;

        template<typename ParentIndexSetHolder_>
        void compute(const ParentIndexSetHolder_& index_set_holder, const std::vector<Index>& facets)
        {
          // In this function, we compute all (cell_dim_-1)-dimensional faces of the mesh, which
          // are adjacent to at least one of the boundary facets given in the facets vector.

          const auto& face_index_set = index_set_holder.template get_index_set<shape_dim_-1, cell_dim_-1>();

          // allocate a temporary vector; this stores the number of cells adjacent to each facet
          std::vector<Index> caf(face_index_set.get_index_bound(), Index(0));

          // loop over all boundary facets
          for(Index i(0); i < Index(facets.size()); ++i)
          {
            for(int j(0); j < face_index_set.get_num_indices(); ++j)
              ++caf[face_index_set(facets[i], j)];
          }

          // count the number of boundary facets
          Index count(0);
          for(Index i(0); i < Index(caf.size()); ++i)
          {
            if(caf[i] > Index(0))
              ++count;
          }

          // allocate boundary facet sets
          _faces.reserve(count);
          for(Index i(0); i < Index(caf.size()); ++i)
          {
            if(caf[i] > Index(0))
              _faces.push_back(i);
          }

          // call base-class
          BaseClass::compute(index_set_holder, facets);
        }

      public:
        Index get_num_entities(int dim)
        {
          if(dim+1 == cell_dim_)
            return Index(_faces.size());
          else if(dim+1 <= cell_dim_)
            return BaseClass::get_num_entities(dim);
          else
            return Index(0);
        }

        void fill_target_sets(TargetSetHolder<Shape_>& target_set_holder)
        {
          TargetSet& target_set = target_set_holder.template get_target_set<cell_dim_-1>();
          XASSERTM(target_set.get_num_entities() == Index(_faces.size()), "invalid target set size");

          for(Index i(0); i < target_set.get_num_entities(); ++i)
            target_set[i] = _faces[i];

          BaseClass::fill_target_sets(target_set_holder);
        }
      }; // class BoundaryFaceComputer<...>

      // Specialisation for cell_dim_ = shape_dim_-1
      template<typename Shape_, int shape_dim_>
      class BoundaryFaceComputer<Shape_, shape_dim_, shape_dim_> :
        public BoundaryFaceComputer<Shape_, shape_dim_, shape_dim_-1>
      {
        typedef BoundaryFaceComputer<Shape_, shape_dim_, shape_dim_-1> BaseClass;

      protected:
        /// a vector of all boundary facet indices
        std::vector<Index> _faces;

      public:
        template<typename ParentIndexSetHolder_>
        explicit BoundaryFaceComputer(const ParentIndexSetHolder_& index_set_holder)
        {
          // Phase 1:
          // Compute all facets, i.e. (n-1)-dimensional faces, which reside on the boundary.
          // We do this by computing how many cells are adjacent to a particular facet.
          // If the number of cells is exactly 1, the facet is a boundary facet.

          const auto& face_index_set = index_set_holder.template get_index_set<shape_dim_, shape_dim_-1>();

          // allocate a temporary vector; this stores the number of cells adjacent to each facet
          std::vector<Index> caf(face_index_set.get_index_bound(), Index(0));

          // loop over all cells of the mesh
          for(Index i(0); i < face_index_set.get_num_entities(); ++i)
          {
            // loop over all facets adjacent to the current cell and increments its counter by one
            for(int j(0); j < face_index_set.get_num_indices(); ++j)
              ++caf[face_index_set(i,j)];
          }

          // count the number of boundary facets
          Index count(0);
          for(Index i(0); i < Index(caf.size()); ++i)
          {
            // for each conformal mesh, a facet must be adjacent to either 1 (boundary facet) or 2 (inner facet) cells.
            ASSERTM((caf[i] > 0) && (caf[i] < 3), "invalid number of cells at facet");
            if(caf[i] == Index(1))
              ++count;
          }

          // allocate boundary facet sets
          _faces.reserve(count);
          for(Index i(0); i < Index(caf.size()); ++i)
          {
            if(caf[i] == Index(1))
              _faces.push_back(i);
          }

          // Phase 2:
          // Compute all lower-dimensional faces adjacent to any of our boundary facets.
          BaseClass::compute(index_set_holder, _faces);
        }

        Index get_num_entities(int dim)
        {
          if(dim+1 == shape_dim_)
            return Index(_faces.size());
          else if(dim+1 <= shape_dim_)
            return BaseClass::get_num_entities(dim);
          else
            return Index(0);
        }

        void fill_target_sets(TargetSetHolder<Shape_>& target_set_holder)
        {
          TargetSet& target_set = target_set_holder.template get_target_set<shape_dim_-1>();
          XASSERTM(target_set.get_num_entities() == Index(_faces.size()), "invalid target set size");

          for(Index i(0); i < target_set.get_num_entities(); ++i)
            target_set[i] = _faces[i];

          BaseClass::fill_target_sets(target_set_holder);
        }
      }; // class BoundaryFaceComputer<...,shape_dim_,shape_dim_-1>

      // This class is required to terminate the inheritance recursion of the class.
      template<typename Shape_, int shape_dim_>
      class BoundaryFaceComputer<Shape_, shape_dim_, 0>
      {
      protected:
        template<typename ParentIndexSetHolder_>
        void compute(const ParentIndexSetHolder_&, const std::vector<Index>&)
        {
          // dummy
        }

      public:
        Index get_num_entities(int)
        {
          return Index(0);
        }

        void fill_target_sets(TargetSetHolder<Shape_>&)
        {
          // dummy
        }
      }; // class BoundaryFaceComputer<...,0>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
#endif // KERNEL_GEOMETRY_INTERN_BOUNDARY_COMPUTER_HPP

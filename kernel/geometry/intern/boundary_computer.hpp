// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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
        typedef std::array<std::vector<int>, std::size_t(shape_dim_)> BndMasksType;
        typedef std::array<std::vector<Index>, std::size_t(shape_dim_)> FaceIdxType;

      protected:
        template<typename ParentIndexSetHolder_>
        void _compute_mask(const ParentIndexSetHolder_& index_set_holder, BndMasksType& bnd_masks, const std::vector<Index>& facets)
        {
          // In this function, we compute all (cell_dim_-1)-dimensional faces of the mesh, which
          // are adjacent to at least one of the boundary facets given in the facets vector.

          // get the mask for this face dimension
          std::vector<int>& mask = bnd_masks.at(std::size_t(cell_dim_-1));

          const auto& face_index_set = index_set_holder.template get_index_set<shape_dim_-1, cell_dim_-1>();

          // allocate a temporary vector; this stores the number of cells adjacent to each facet
          mask.clear();
          mask.resize(face_index_set.get_index_bound(), 0);

          // loop over all boundary facets
          for(Index i(0); i < Index(facets.size()); ++i)
          {
            for(int j(0); j < face_index_set.get_num_indices(); ++j)
              ++mask[face_index_set(facets[i], j)];
          }

          BaseClass::_compute_mask(index_set_holder, bnd_masks, facets);
        }

        void _compute_faces(const BndMasksType& bnd_masks, FaceIdxType& face_idx)
        {
          // get the mask for this face dimension
          const std::vector<int>& mask = bnd_masks.at(std::size_t(cell_dim_-1));
          std::vector<Index>& faces = face_idx.at(std::size_t(cell_dim_-1));

          // count the number of boundary faces
          Index count(0);
          for(Index i(0); i < Index(mask.size()); ++i)
          {
            if(mask[i] != 0)
              ++count;
          }

          // allocate boundary face sets
          faces.reserve(count);
          for(Index i(0); i < Index(mask.size()); ++i)
          {
            if(mask[i] != 0)
              faces.push_back(i);
          }

          // call base-class
          BaseClass::_compute_faces(bnd_masks, face_idx);
        }

      public:
        void fill_target_sets(TargetSetHolder<Shape_>& target_set_holder, const FaceIdxType& face_idx)
        {
          const std::vector<Index>& faces = face_idx.at(std::size_t(cell_dim_-1)) ;

          TargetSet& target_set = target_set_holder.template get_target_set<cell_dim_-1>();
          XASSERTM(target_set.get_num_entities() == Index(faces.size()), "invalid target set size");

          for(Index i(0); i < target_set.get_num_entities(); ++i)
            target_set[i] = faces[i];

          BaseClass::fill_target_sets(target_set_holder, face_idx);
        }

        Index build_halo_buffer(std::vector<int>& buffer, const BndMasksType& bnd_masks, const TargetSetHolder<Shape_>& halo_trg_holder)
        {
          Index offset = BaseClass::build_halo_buffer(buffer, bnd_masks, halo_trg_holder);
          const TargetSet& target_set = halo_trg_holder.template get_target_set<cell_dim_-1>();
          const std::vector<int>& mask = bnd_masks.at(std::size_t(cell_dim_-1));
          for(Index i(0); i < target_set.get_num_entities(); ++i, ++offset)
            buffer[offset] = mask[target_set[i]];
          return offset;
        }

        Index mask_halo_buffer(BndMasksType& bnd_masks, const std::vector<int>& buffer, const TargetSetHolder<Shape_>& halo_trg_holder)
        {
          Index offset = BaseClass::mask_halo_buffer(bnd_masks, buffer, halo_trg_holder);
          const TargetSet& target_set = halo_trg_holder.template get_target_set<cell_dim_-1>();
          std::vector<int>& mask = bnd_masks.at(std::size_t(cell_dim_-1));
          for(Index i(0); i < target_set.get_num_entities(); ++i, ++offset)
          {
            if(buffer[offset] != 0)
              mask[target_set[i]] = 1;
          }
          return offset;
        }
      }; // class BoundaryFaceComputer<...>

      // Specialization for cell_dim_ = shape_dim_-1
      template<typename Shape_, int shape_dim_>
      class BoundaryFaceComputer<Shape_, shape_dim_, shape_dim_> :
        public BoundaryFaceComputer<Shape_, shape_dim_, shape_dim_-1>
      {
        typedef BoundaryFaceComputer<Shape_, shape_dim_, shape_dim_-1> BaseClass;
        typedef std::array<std::vector<int>, std::size_t(shape_dim_)> BndMasksType;
        typedef std::array<std::vector<Index>, std::size_t(shape_dim_)> FaceIdxType;

      public:
        template<typename ParentIndexSetHolder_>
        void compute_all(const ParentIndexSetHolder_& index_set_holder, BndMasksType& bnd_masks, FaceIdxType& face_idx)
        {
          // get the mask for this face dimension
          std::vector<int>& mask = bnd_masks.at(std::size_t(shape_dim_-1));
          std::vector<Index>& faces = face_idx.at(std::size_t(shape_dim_-1));

          // Phase 1:
          // Compute all facets, i.e. (n-1)-dimensional faces, which reside on the boundary.
          // We do this by computing how many cells are adjacent to a particular facet.
          // If the number of cells is exactly 1, the facet is a boundary facet.

          const auto& face_index_set = index_set_holder.template get_index_set<shape_dim_, shape_dim_-1>();

          // allocate a temporary vector; this stores the number of cells adjacent to each facet
          mask.clear();
          mask.resize(face_index_set.get_index_bound(), 0);

          // loop over all cells of the mesh
          for(Index i(0); i < face_index_set.get_num_entities(); ++i)
          {
            // loop over all facets adjacent to the current cell and increments its counter by one
            for(int j(0); j < face_index_set.get_num_indices(); ++j)
              ++mask[face_index_set(i,j)];
          }

          // count the number of boundary facets
          Index count(0);
          for(Index i(0); i < Index(mask.size()); ++i)
          {
            // for each conformal mesh, a facet must be adjacent to either 1 (boundary facet) or 2 (inner facet) cells.
            ASSERTM((mask[i] > 0) && (mask[i] < 3), "invalid number of cells at facet");
            if(mask[i] == 1)
              ++count;
          }

          // allocate boundary facet sets
          faces.reserve(count);
          for(Index i(0); i < Index(mask.size()); ++i)
          {
            if(mask[i] == 1)
              faces.push_back(i);
          }

          // Phase 2:
          // Compute all lower-dimensional faces adjacent to any of our boundary facets.
          BaseClass::_compute_mask(index_set_holder, bnd_masks, face_idx.back());
          BaseClass::_compute_faces(bnd_masks, face_idx);
        }

        template<typename ParentIndexSetHolder_>
        void compute_masks(const ParentIndexSetHolder_& index_set_holder, BndMasksType& bnd_masks, FaceIdxType& face_idx, const std::vector<int>& facet_mask)
        {
          // get the mask for this face dimension
          std::vector<int>& mask = bnd_masks.at(std::size_t(shape_dim_-1));
          std::vector<Index>& faces = face_idx.at(std::size_t(shape_dim_-1));

          // Phase 1:
          // Compute all facets, i.e. (n-1)-dimensional faces, which reside on the boundary.
          // We do this by computing how many cells are adjacent to a particular facet.
          // If the number of cells is exactly 1, the facet is a boundary facet.

          const auto& face_index_set = index_set_holder.template get_index_set<shape_dim_, shape_dim_-1>();

          // allocate a temporary vector; this stores the number of cells adjacent to each facet
          mask.clear();
          mask.resize(face_index_set.get_index_bound(), 0);

          // loop over all cells of the mesh
          for(Index i(0); i < face_index_set.get_num_entities(); ++i)
          {
            // loop over all facets adjacent to the current cell and increments its counter by one
            for(int j(0); j < face_index_set.get_num_indices(); ++j)
              ++mask[face_index_set(i,j)];
          }

          // reset all masked facets to 0, so that they won't be part of the created boundary mesh part
          XASSERTM(facet_mask.size() == mask.size(), "invalid mask vector size");
          for(std::size_t i(0); i < mask.size(); ++i)
            mask[i] = (facet_mask[i] != 0 ? 0 : mask[i]);

          // count the number of boundary facets
          Index count(0);
          for(Index i(0); i < Index(mask.size()); ++i)
          {
            // for each conformal mesh, a facet must be adjacent to either 1 (boundary facet) or 2 (inner facet) cells.
            ASSERTM((facet_mask[i] != 0) || ((mask[i] > 0) && (mask[i] < 3)), "invalid number of cells at facet");
            if(mask[i] == 1)
              ++count;
          }

          // allocate boundary face sets
          faces.reserve(count);
          for(Index i(0); i < Index(mask.size()); ++i)
          {
            if(mask[i] == 1)
              faces.push_back(i);
          }

          BaseClass::_compute_mask(index_set_holder, bnd_masks, faces);
        }

        void compute_faces(const BndMasksType& bnd_masks, FaceIdxType& face_idx)
        {
          // Compute all lower-dimensional faces adjacent to any of our boundary facets.
          BaseClass::_compute_faces(bnd_masks, face_idx);
        }

        void fill_target_sets(TargetSetHolder<Shape_>& target_set_holder, const FaceIdxType& face_idx)
        {
          const std::vector<Index>& faces = face_idx.at(std::size_t(shape_dim_-1));

          TargetSet& target_set = target_set_holder.template get_target_set<shape_dim_-1>();
          XASSERTM(target_set.get_num_entities() == Index(faces.size()), "invalid target set size");

          for(Index i(0); i < target_set.get_num_entities(); ++i)
            target_set[i] = faces[i];

          BaseClass::fill_target_sets(target_set_holder, face_idx);
        }

        Index build_halo_buffer(std::vector<int>& buffer, const BndMasksType& bnd_masks, const TargetSetHolder<Shape_>& halo_trg_holder)
        {
          Index offset = BaseClass::build_halo_buffer(buffer, bnd_masks, halo_trg_holder);
          const TargetSet& target_set = halo_trg_holder.template get_target_set<shape_dim_-1>();
          const std::vector<int>& mask = bnd_masks.at(std::size_t(shape_dim_-1));
          for(Index i(0); i < target_set.get_num_entities(); ++i, ++offset)
            buffer[offset] = mask[target_set[i]];
          return offset;
        }

        Index mask_halo_buffer(BndMasksType& bnd_masks, const std::vector<int>& buffer, const TargetSetHolder<Shape_>& halo_trg_holder)
        {
          Index offset = BaseClass::mask_halo_buffer(bnd_masks, buffer, halo_trg_holder);
          const TargetSet& target_set = halo_trg_holder.template get_target_set<shape_dim_-1>();
          std::vector<int>& mask = bnd_masks.at(std::size_t(shape_dim_-1));
          for(Index i(0); i < target_set.get_num_entities(); ++i, ++offset)
          {
            if(buffer[offset] == 1)
              mask[target_set[i]] = 1;
          }
          return offset;
        }
      }; // class BoundaryFaceComputer<...,shape_dim_,shape_dim_-1>

      // This class is required to terminate the inheritance recursion of the class.
      template<typename Shape_, int shape_dim_>
      class BoundaryFaceComputer<Shape_, shape_dim_, 0>
      {
      protected:
        typedef std::array<std::vector<int>, std::size_t(shape_dim_)> BndMasksType;
        typedef std::array<std::vector<Index>, std::size_t(shape_dim_)> FaceIdxType;

        template<typename ParentIndexSetHolder_>
        void _compute_mask(const ParentIndexSetHolder_&, BndMasksType&, const std::vector<Index>&)
        {
          // dummy
        }

        void _compute_faces(const BndMasksType&, FaceIdxType&)
        {
          // dummy
        }

      public:
        void fill_target_sets(TargetSetHolder<Shape_>&, const FaceIdxType&)
        {
          // dummy
        }

        Index build_halo_buffer(std::vector<int>&, const BndMasksType&, const TargetSetHolder<Shape_>&)
        {
          return 0u;
        }

        Index mask_halo_buffer(BndMasksType&, const std::vector<int>&, const TargetSetHolder<Shape_>&)
        {
          return 0u;
        }
      }; // class BoundaryFaceComputer<...,0>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

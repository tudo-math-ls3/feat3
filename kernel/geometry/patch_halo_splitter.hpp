// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_PATCH_HALO_SPLITTER_HPP
#define KERNEL_GEOMETRY_PATCH_HALO_SPLITTER_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/adjacency/graph.hpp>

#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_, int shape_dim_ = Shape_::dimension>
      class PatchInvMapWrapper;

      template<typename Shape_>
      class PatchHaloSplitPart;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Base-Mesh Patch Halo splitter
     */
    template<typename Mesh_>
    class PatchHaloSplitter;

    /**
     * \brief Base-Mesh Patch Halo splitter implementation for ConformalMesh
     *
     * This class is used to split the base-mesh halos in the case of fully recursive partitioning
     * to obtain the inter-patch halos of the child patches.
     *
     * \author Peter Zajac
     */
    template<typename Shape_, int num_coords_, typename Coord_>
    class PatchHaloSplitter<ConformalMesh<Shape_, num_coords_, Coord_>> :
      public Geometry::Factory<Geometry::MeshPart<Geometry::ConformalMesh<Shape_, num_coords_, Coord_>>>
    {
    public:
      /// Our shape type
      typedef Shape_ ShapeType;
      /// The shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// The mesh type
      typedef Geometry::ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
      /// The MeshPart type
      typedef Geometry::MeshPart<MeshType> MeshPartType;

    protected:
      typedef Intern::PatchHaloSplitPart<ShapeType> SplitHaloType;
      /// The BaseMesh we refer to
      const MeshType& _base_mesh;
      /// The MeshPart identifying this patch
      const MeshPartType& _patch_mesh_part;
      /// inverse patch mapping
      std::vector<std::vector<Index>> _patch_inv_map;
      /// split halos
      std::map<int, std::shared_ptr<SplitHaloType>> _split_halos;

      /// intersected halo indices
      std::vector<std::vector<Index>> _isect_halo_idx;

    public:
      /**
       * \brief Creates the halo splitter
       *
       * \param[in] base_mesh
       * A \resident reference to the base-mesh whose halos are to be split.
       *
       * \param[in] patch_mesh_part
       * A \resident reference to the mesh-part that defines the patch of the partitioned base-mesh.
       */
      explicit PatchHaloSplitter(const MeshType& base_mesh, const MeshPartType& patch_mesh_part) :
        _base_mesh(base_mesh),
        _patch_mesh_part(patch_mesh_part),
        _patch_inv_map(std::size_t(shape_dim+1)),
        _isect_halo_idx(std::size_t(shape_dim+1))
      {
        // allocate inverse patch mapping
        _patch_inv_map.resize(std::size_t(shape_dim+1));
        for(int dim(0); dim <= shape_dim; ++dim)
        {
          _patch_inv_map.at(std::size_t(dim)).resize(_base_mesh.get_num_entities(dim), ~Index(0));
        }
        // build inverse patch mapping
        Intern::PatchInvMapWrapper<ShapeType>::build(_patch_inv_map, patch_mesh_part.get_target_set_holder());
      }

      /**
       * \brief Adds a base-mesh halo to the splitter and returns the required buffer size.
       *
       * \param[in] halo_rank
       * The rank of the neighbor process that is adjacent to this halo wrt the parent partitioning.
       *
       * \param[in] halo
       * A \transient reference to the halo mesh-part that is to be split.
       *
       * \returns
       * The size of the buffer that is required to store the serialized split halo.
       * If the halo and the patch are disjoint, then the split halo is empty and this function returns 0.
       */
      std::size_t add_halo(int halo_rank, const MeshPartType& halo)
      {
        // create new split halo object
        auto sh = std::make_shared<SplitHaloType>();

        // If the following function returns false, then this patch is not adjacent to that halo
        if(!sh->split(_patch_inv_map, halo.get_target_set_holder()))
          return std::size_t(0);

        // okay, store halo and return split size
        return _split_halos.emplace(halo_rank, sh).first->second->size();
      }

      /**
       * \brief Serializes a split halo into a buffer.
       *
       * \param[in] halo_rank
       * The rank of the halo whose split data is to be serialized.
       *
       * \param[in] child_rank
       * The rank of the neighbor process that is adjacent to this halo wrt to the child partitioning.
       *
       * \returns
       * A vector containing the serialized split halo data.
       */
      std::vector<Index> serialize_split_halo(int halo_rank, int child_rank) const
      {
        const auto it = _split_halos.find(halo_rank);
        XASSERT(it != _split_halos.end());

        std::vector<Index> buffer;
        it->second->serialize(buffer, child_rank);
        return buffer;
      }

      /**
       * \brief Intersects a foreign split halo with the split halo of this patch
       *
       * \param[in] halo_rank
       * The rank of the halo whose split data is to be intersected with the foreign split halo.
       *
       * \param[in] buffer
       * The buffer vector containing the serialized split halo data of the potential neighbor child processes.
       *
       * \param[in] buffer_offset
       * The offset of the current neighbor split halo data within the buffer vector.
       *
       * \returns
       * \c true, if the split halo of this patch and the current neighbor process intersect, or
       * \c false, if the split halos are disjoint.
       */
      bool intersect_split_halo(int halo_rank, const std::vector<Index>& buffer, const Index buffer_offset)
      {
        // clear all intersection indices
        for(auto& x : _isect_halo_idx)
          x.clear();

        const auto it = _split_halos.find(halo_rank);
        if(it == _split_halos.end())
          return false;

        // try to intersect
        return it->second->intersect(_isect_halo_idx, buffer, std::size_t(buffer_offset));
      }

      /* *************************************************************************************** */
      /* F A C T O R Y   I N T E R F A C E   I M P L E M E N T A T I O N                         */
      /* *************************************************************************************** */

      virtual Index get_num_entities(int dim) override
      {
        return Index(_isect_halo_idx.at(std::size_t(dim)).size());
      }

      virtual void fill_attribute_sets(typename MeshPartType::AttributeSetContainer& /*attribute_container*/) override
      {
        // no attributes for halo mesh-parts
      }

      virtual void fill_index_sets(std::unique_ptr<typename MeshPartType::IndexSetHolderType>& index_set_holder) override
      {
        XASSERT(index_set_holder.get() == nullptr);
        // no index sets for halo mesh-parts
      }

      virtual void fill_target_sets(typename MeshPartType::TargetSetHolderType& target_set_holder) override
      {
        Intern::PatchInvMapWrapper<ShapeType>::fill(target_set_holder, _isect_halo_idx);
      }
    }; // class PatchHaloSplitter<ConformalMesh<...>>

    /// \cond internal
    namespace Intern
    {
      //template<int shape_dim_>
      class PatchInvMap
      {
      public:
        static void build(std::vector<Index>& pim, const TargetSet& ts)
        {
          for(Index i(0); i < ts.get_num_entities(); ++i)
            pim[ts[i]] = i;
        }

        static bool split(std::vector<Index>& shi, std::vector<Index>& spi, const std::vector<Index>& pim, const TargetSet& ts)
        {
          shi.reserve(ts.get_num_entities());
          spi.reserve(ts.get_num_entities());

          for(Index i(0); i < ts.get_num_entities(); ++i)
          {
            const Index patch_idx = pim[ts[i]];
            if(patch_idx != ~Index(0))
            {
              spi.push_back(patch_idx);
              shi.push_back(i);
            }
          }

          return !spi.empty();
        }

        static void fill(TargetSet& trg, const std::vector<Index>& idx)
        {
          const Index n = trg.get_num_entities();
          XASSERT(n == Index(idx.size()));
          for(Index i(0); i < n; ++i)
            trg[i] = idx[std::size_t(i)];
        }
      }; // class PatchInvMap

      template<typename Shape_, int shape_dim_>
      class PatchInvMapWrapper
      {
      public:
        static_assert(shape_dim_ > 0, "invalid shape dimension");
        static void build(std::vector<std::vector<Index>>& pim, const TargetSetHolder<Shape_>& tsh)
        {
          PatchInvMapWrapper<Shape_, shape_dim_-1>::build(pim, tsh);
          PatchInvMap::build(pim.at(std::size_t(shape_dim_)), tsh.template get_target_set<shape_dim_>());
        }

        static bool split(
          std::vector<std::vector<Index>>& shi, std::vector<std::vector<Index>>& spi,
          const std::vector<std::vector<Index>>& pim, const TargetSetHolder<Shape_>& tsh)
        {
          bool have = PatchInvMapWrapper<Shape_, shape_dim_-1>::split(shi, spi, pim, tsh);
          const std::size_t k = std::size_t(shape_dim_);
          return PatchInvMap::split(shi.at(k), spi.at(k), pim.at(k), tsh.template get_target_set<shape_dim_>()) || have;
        }

        static void fill(TargetSetHolder<Shape_>& tsh, const std::vector<std::vector<Index>>& idx)
        {
          PatchInvMapWrapper<Shape_, shape_dim_-1>::fill(tsh, idx);
          PatchInvMap::fill(tsh.template get_target_set<shape_dim_>(), idx.at(std::size_t(shape_dim_)));
        }
      };

      template<typename Shape_>
      class PatchInvMapWrapper<Shape_, 0>
      {
      public:
        static void build(std::vector<std::vector<Index>>& pim, const TargetSetHolder<Shape_>& tsh)
        {
          PatchInvMap::build(pim.at(std::size_t(0)), tsh.template get_target_set<0>());
        }
        static bool split(
          std::vector<std::vector<Index>>& shi, std::vector<std::vector<Index>>& spi,
          const std::vector<std::vector<Index>>& pim, const TargetSetHolder<Shape_>& tsh)
        {
          const std::size_t k = std::size_t(0);
          return PatchInvMap::split(shi.at(k), spi.at(k), pim.at(k), tsh.template get_target_set<0>());
        }

        static void fill(TargetSetHolder<Shape_>& tsh, const std::vector<std::vector<Index>>& idx)
        {
          PatchInvMap::fill(tsh.template get_target_set<0>(), idx.at(std::size_t(0)));
        }
      };

      template<typename Shape_>
      class PatchHaloSplitPart
      {
      protected:
        static constexpr int shape_dim = Shape_::dimension;

        /// split halo indices w.r.t the base mesh halo
        std::vector<std::vector<Index>> _halo_idx;
        /// split halo indices w.r.t the patch mesh part
        std::vector<std::vector<Index>> _patch_idx;

      public:
        explicit PatchHaloSplitPart()
        {
          std::size_t sdim = std::size_t(shape_dim+1);
          _halo_idx.resize(sdim);
          _patch_idx.resize(sdim);
        }

        bool split(const std::vector<std::vector<Index>>& pim, const TargetSetHolder<Shape_>& tsh)
        {
          return PatchInvMapWrapper<Shape_>::split(_halo_idx, _patch_idx, pim, tsh);
        }

        std::size_t size() const
        {
          std::size_t s = std::size_t(0);
          for(const auto& x : _halo_idx)
            s += x.size();

          // s == 0 ==> patch and halo are disjoint
          return s == std::size_t(0) ? s : ++s + _halo_idx.size();
        }

        void serialize(std::vector<Index>& buffer, int child_rank)
        {
          const std::size_t s = this->size();
          if(this->size() == std::size_t(0))
            return;

          if(buffer.capacity() < buffer.size() + s)
            buffer.reserve(buffer.size() + s);

          buffer.push_back(Index(child_rank));

          for(std::size_t i(0); i < _halo_idx.size(); ++i)
            buffer.push_back(Index(_halo_idx.at(i).size()));

          for(std::size_t i(0); i < _halo_idx.size(); ++i)
          {
            const std::vector<Index>& v = _halo_idx.at(i);
            for(std::size_t j(0); j < v.size(); ++j)
              buffer.push_back(v[j]);
          }
        }

        bool intersect(std::vector<std::vector<Index>>& isect, const std::vector<Index>& other, const std::size_t offset) const
        {
          // offset of first index within other halo
          std::size_t off = offset + std::size_t(shape_dim+2);

          bool ret = false;

          // compute offsets from sizes
          for(std::size_t i(0); i <= std::size_t(shape_dim); ++i)
          {
            // get index counts
            const std::size_t n_1 = _halo_idx.at(i).size();
            const std::size_t n_2 = other.at(offset+i+1u);

            // get the minimum of both
            const std::size_t n_min = Math::min(n_1, n_2);
            if(n_min == std::size_t(0))
            {
              off += n_2;
              continue;
            }

            // get my halo indices
            const Index* idx_1 = _halo_idx.at(i).data();
            const Index* idx_2 = &other.data()[off];
            const Index* idx_p = _patch_idx.at(i).data();

            // intersect idx_1 and idx_2
            std::vector<Index>& idx = isect.at(i);
            idx.reserve(n_min);
            for(std::size_t j_1(0), j_2(0); (j_1 < n_1) && (j_2 < n_2); )
            {
              if(idx_1[j_1] < idx_2[j_2])
                ++j_1;
              else if(idx_1[j_1] > idx_2[j_2])
                ++j_2;
              else
              {
                // match; store halo index w.r.t. patch
                idx.push_back(idx_p[j_1]);
                ++j_1;
                ++j_2;
                ret = true;
              }
            }

            off += n_2;
          }

          // okay
          return ret;
        }
      };
    } // namespace Intern
    /// \endcond

  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PATCH_HALO_SPLITTER_HPP

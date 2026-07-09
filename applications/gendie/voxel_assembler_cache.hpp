// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once
#include <kernel/voxel_assembly/burgers_velo_material_assembler.hpp>

#include <deque>
#include <map>
#include <memory>
#include <optional>

namespace Gendie
{
  namespace Intern
  {
    /**
     * \brief Per-level cache of voxel Burgers velocity-material assemblers.
     *
     * Constructing a \c VoxelBurgersVeloMaterialAssembler is expensive and mesh-bound: it builds a
     * LagrangeDataHandler (node/cell-to-dof/sorter buffers, connectivity, coloring ingest) plus CUDA
     * launch tuning. Everything it derives depends only on the finite element \c space and the element
     * \c coloring, both of which are fixed for a gendie run, so the object can be built once per level
     * and reused across nonlinear iterations; only the cheap scalar material parameters have to be
     * re-set per assembly call.
     *
     * The gendie assembler classes are templated only on the domain, but the voxel assembler type
     * additionally depends on the space, data type and index type - and \c _assemble_matrices is
     * instantiated for multiple matrix precisions (float and double, see SteadyFlowSolverFactory).
     * A single fixed-type member would therefore not compile, so the cache is type-erased and keyed
     * per concrete assembler type. The key is the address of a per-type \c static tag object, which
     * gives a unique, RTTI-free identity for each instantiated type.
     */
    class VoxelAssemblerCache
    {
      template<typename T_>
      struct TypeTag { static inline const char id{}; };

      std::map<const void*, std::shared_ptr<void>> _caches;

    public:
      /**
       * \brief Returns the cached voxel assembler for the given level, constructing it on first use.
       *
       * \param[in] level   The multigrid level index (0 = finest).
       * \param[in] space   The finite element space of that level.
       * \param[in] coloring The element coloring of that level.
       */
      template<typename Space_, typename DT_, typename IT_, typename Coloring_>
      FEAT::VoxelAssembly::VoxelBurgersVeloMaterialAssembler<Space_, DT_, IT_>&
      get(std::size_t level, const Space_& space, const Coloring_& coloring)
      {
        typedef FEAT::VoxelAssembly::VoxelBurgersVeloMaterialAssembler<Space_, DT_, IT_> AsmType;
        typedef std::deque<std::optional<AsmType>> CacheType;

        const void* key = &TypeTag<AsmType>::id;
        auto it = _caches.find(key);
        std::shared_ptr<CacheType> cache;
        if(it == _caches.end())
        {
          cache = std::make_shared<CacheType>();
          _caches.emplace(key, cache);
        }
        else
        {
          cache = std::static_pointer_cast<CacheType>(it->second);
        }

        if(cache->size() <= level)
          cache->resize(level + std::size_t(1));

        auto& slot = (*cache)[level];
        if(!slot.has_value())
          slot.emplace(space, coloring); // move-only type: emplace constructs it in place

        return *slot;
      }

      /// Releases all cached assemblers (freeing their mesh/coloring/CUDA resources).
      void clear()
      {
        _caches.clear();
      }
    };
  } // namespace Intern
} // namespace Gendie

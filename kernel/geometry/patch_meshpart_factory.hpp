#pragma once
#ifndef KERNEL_GEOMETRY_PATCH_MESHPART_FACTORY_HPP
#define KERNEL_GEOMETRY_PATCH_MESHPART_FACTORY_HPP 1

// includes, FEAT
#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/mesh_part.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief  Factory for building a PatchBaseMesh from a given partitioning
     *
     * \tparam MeshType_
     * Type of the BaseMesh
     *
     * \author Jordi Paul
     */
    template<typename MeshType_>
    class PatchMeshPartFactory : public Factory< MeshPart<MeshType_> >
    {
      public:
        /// The shape type of the mesh
        typedef typename MeshType_::ShapeType ShapeType;
        /// Mesh typedef
        typedef MeshPart<MeshType_> MeshPartType;
        /// Data type for attributes
        typedef typename MeshType_::VertexSetType::CoordType AttributeDataType;
        /// Mesh attribute holder type
        typedef typename MeshPartType::MeshAttributeContainer MeshAttributeContainer;
        /// Index set holder type
        typedef typename MeshPartType::IndexSetHolderType IndexSetHolderType;
        /// Target set holder type
        typedef typename MeshPartType::TargetSetHolderType TargetSetHolderType;
        /// The shape dimension
        static constexpr int shape_dim = ShapeType::dimension;

      private:
        /// Number of entities for each shape dimension
        Index _num_entities[shape_dim + 1];
        /// This will hold all BaseMesh cell numbers of cells present in the patch
        std::vector<Index> _cells_patch;

      public:
        /**
         * \brief Builds a factory from a given partitioning
         *
         * \param[in] my_rank
         * Identifier for this patch. Usually this is the MPI rank
         *
         * \param[in] ranks_at_cells
         * For all cells in the BaseMesh, this contains the patches it is present in (this means that true overlap
         * is supported)
         *
         */
        explicit PatchMeshPartFactory(Index my_rank, const Adjacency::Graph& ranks_at_cells)
        {
          // The patch initially has no entities
          for(int i(0); i < shape_dim + 1; ++i)
            _num_entities[i] = 0;

          // Number of cells in the BaseMesh
          Index ncells_base(ranks_at_cells.get_num_nodes_domain());

          // For each cell in the BaseMesh, find out if it belongs to this patch
          for(Index cell(0); cell < ncells_base; ++cell)
          {
            for(auto it(ranks_at_cells.image_begin(cell)); it != ranks_at_cells.image_end(cell); ++it)
            {
              if(*it == my_rank)
              {
                _cells_patch.push_back(cell);
                continue;
              }
            }
          }

          // Now we know the number of cells
          _num_entities[shape_dim] = Index(_cells_patch.size());
        }

        /// virtual destructor
        virtual ~PatchMeshPartFactory()
        {
          _cells_patch.clear();
        }

        /**
         * \brief Returns the number of entities.
         *
         * \param[in] dim
         * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= shape_dim.
         *
         * \returns
         * The number of entities of dimension \p dim.
         */
        virtual Index get_num_entities(int dim) override
        {
          XASSERT(dim >= 0);
          XASSERT(dim <= shape_dim);

          return _num_entities[dim];
        }

        /**
         * \brief Pretends to fill the attribute sets
         *
         * \param[in,out] attribute_set_holder
         * The attribute set holder whose attribute sets are to be filled.
         *
         * This MeshPart does not have any attributes.
         */
        virtual void fill_attribute_sets(MeshAttributeContainer& DOXY(attribute_set_holder)) override
        {
          // Do nothing
          return;
        }

        /**
         * \brief Pretends to fill the index sets (=topology)
         *
         * \param[in,out] index_set_holder
         * The index set holder whose index sets are to be filled.
         *
         * This MeshPart does not have a topology.
         */
        virtual void fill_index_sets(IndexSetHolderType*& DOXY(index_set_holder)) override
        {
          // Do nothing
          return;
        }

        /**
         * \brief Fills the target sets.
         *
         * \param[in,out] target_set_holder
         * The target set holder whose target sets are to be filled.
         *
         * This fills just the TargetSet corresponding to the cell dimension.
         */
        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          Index ncells_patch(get_num_entities(shape_dim));

          TargetSet new_target_set(ncells_patch);
          for(Index k(0); k < ncells_patch; ++k)
            new_target_set[k] = _cells_patch[k];

          target_set_holder.template get_target_set<shape_dim>() = std::move(new_target_set);

          return;
        }

    }; // class PatchMeshPartFactory<MeshPart<...>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PATCH_MESHPART_FACTORY_HPP

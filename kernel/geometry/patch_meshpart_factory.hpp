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
     * \author Peter Zajac
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
        /// This will hold all BaseMesh cell numbers of cells present in the patch
        std::vector<Index> _cells_patch;

      public:
        /**
         * \brief Builds a factory from a given partitioning
         *
         * \param[in] my_rank
         * Identifier for this patch. Usually this is the MPI rank
         *
         * \param[in] elems_at_rank
         * The graph representing the elements-at-rank adjacency.
         */
        explicit PatchMeshPartFactory(Index my_rank, const Adjacency::Graph& elems_at_rank)
        {
          XASSERT(my_rank < elems_at_rank.get_num_nodes_domain());

          // extract the elements of our rank
          _cells_patch.reserve(elems_at_rank.degree(my_rank));
          for(auto it = elems_at_rank.image_begin(my_rank); it != elems_at_rank.image_end(my_rank); ++it)
            _cells_patch.push_back(*it);
        }

        /// virtual destructor
        virtual ~PatchMeshPartFactory()
        {
        }

        /**
         * \brief Checks whether the patch is empty.
         */
        bool empty() const
        {
          return _cells_patch.empty();
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
          return (dim < shape_dim ? Index(0) : Index(_cells_patch.size()));
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
          TargetSet& target_set = target_set_holder.template get_target_set<shape_dim>();
          for(Index k(0); k < Index(_cells_patch.size()); ++k)
            target_set[k] = _cells_patch[k];
        }
    }; // class PatchMeshPartFactory<MeshPart<...>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PATCH_MESHPART_FACTORY_HPP

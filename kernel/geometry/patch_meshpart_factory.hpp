#pragma once
#ifndef KERNEL_GEOMETRY_PATCH_MESHPART_FACTORY_HPP
#define KERNEL_GEOMETRY_PATCH_MESHPART_FACTORY_HPP 1

// includes, FEAST
#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/mesh_part.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
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
        typedef MeshAttributeHolder<ShapeType, AttributeDataType> AttributeHolderType;
        /// index set holder type
        typedef typename MeshPartType::IndexSetHolderType IndexSetHolderType;
        /// target set holder type
        typedef typename MeshPartType::TargetSetHolderType TargetSetHolderType;

        static constexpr int shape_dim = ShapeType::dimension;

      private:
        /// number of entities for fine mesh
        Index _num_entities[shape_dim + 1];
        std::vector<Index> _cells_patch;

      public:
        explicit PatchMeshPartFactory(Index my_rank, const MeshType_& DOXY(base_mesh_),
          const Adjacency::Graph& ranks_at_cells)
        {
          for(int i(0); i < shape_dim + 1; ++i)
            _num_entities[i] = 0;

          Index ncells_base(ranks_at_cells.get_num_nodes_domain());

          for(Index cell(0); cell < ncells_base; ++cell)
          {
            Adjacency::Graph::ImageIterator it(ranks_at_cells.image_begin(cell));
            Adjacency::Graph::ImageIterator jt(ranks_at_cells.image_end(cell));
            for(; it != jt; ++it)
            {
              if(*it == my_rank)
              {
                _cells_patch.push_back(cell);
                continue;
              }
            }
          }

          _num_entities[shape_dim] = Index(_cells_patch.size());
        }

        /// virtual destructor
        virtual ~PatchMeshPartFactory()
        {
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
          ASSERT_(dim >= 0); //, "There are no entities of dim "+stringify(dim)+" < 0!");
          ASSERT_(dim <= shape_dim); //, "There are no entities of dim "+stringify(dim)+" > "
              //+stringify(shape_dim)+" = ShapeType::dimension");

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
        virtual void fill_attribute_sets(AttributeHolderType& DOXY(attribute_set_holder)) override
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

        /**
         * \brief Returns the name of the MeshPart this factory will construct
         *
         * \returns The name
         */
        virtual String get_identifier() const override
        {
          return "root_meshpart";
        }

        /**
         * \brief Returns the name of the parent mesh
         *
         * \returns The name of the parent mesh the constructed MeshPart will refer to
         */
        virtual String get_parent_identifier() const override
        {
          return "root";
        };
    }; // class Factory<MeshPart<...>>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_PATCH_MESHPART_FACTORY_HPP

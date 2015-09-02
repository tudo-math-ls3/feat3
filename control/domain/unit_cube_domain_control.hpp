#pragma once
#ifndef CONTROL_DOMAIN_UNIT_CUBE_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_UNIT_CUBE_DOMAIN_CONTROL_HPP 1

#include <control/domain/domain_control.hpp>
#include <kernel/geometry/unit_cube_patch_generator.hpp>

namespace FEAST
{
  namespace Control
  {
    namespace Domain
    {
      template<typename Mesh_>
      class UnitCubeDomainControl :
        public DomainControl<Mesh_>
      {
      public:
        typedef DomainControl<Mesh_> BaseClass;
        typedef DomainLayer<Mesh_> LayerType;
        typedef DomainLevel<Mesh_> LevelType;
        typedef Mesh_ MeshType;
        typedef Geometry::MeshAtlas<MeshType> AtlasType;

      public:
        explicit UnitCubeDomainControl(int rank, int nprocs, int lvl_max, int lvl_min = -1) :
          BaseClass()
        {
          typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

          // communication data structures
          std::vector<Index> ranks, ctags;

          // create root mesh node
          MeshNodeType* mesh_node = nullptr;
          int lvl = Geometry::UnitCubePatchGenerator<MeshType>::create(rank, Math::max(nprocs,1), mesh_node, ranks, ctags);

          // push layer
          this->_layers.push_back(new LayerType(std::move(ranks), std::move(ctags)));

          // adjust lvl_max and lvl_min
          lvl_max = Math::max(lvl_max, 0);
          if(lvl_min < 0)
            lvl_min = Math::max(lvl_max + lvl_min + 1, 0);
          else
            lvl_min = Math::min(lvl_min, lvl_max);

          // refine up to desired minimum level
          for(; lvl < lvl_min; ++lvl)
          {
            MeshNodeType* coarse_node = mesh_node;
            mesh_node = coarse_node->refine();
            delete coarse_node;
          }

          // add coarse mesh node
          this->_levels.push_back(new LevelType(lvl, mesh_node));

          // refine up to desired maximum level
          for(; lvl < lvl_max;)
          {
            MeshNodeType* coarse_node = mesh_node;
            mesh_node = coarse_node->refine();
            this->_levels.push_back(new LevelType(++lvl, mesh_node));
          }
        }
      }; // class UnitCubeDomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAST

#endif // CONTROL_DOMAIN_UNIT_CUBE_DOMAIN_CONTROL_HPP

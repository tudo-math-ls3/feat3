#pragma once
#ifndef CONTROL_DOMAIN_MESH_STREAMER_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_MESH_STREAMER_DOMAIN_CONTROL_HPP 1

#include <control/domain/domain_control.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_node.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Domain
    {
      template<typename Mesh_>
      class MeshStreamerDomainControl :
        public DomainControl<Mesh_>
      {
      public:
        typedef DomainControl<Mesh_> BaseClass;
        typedef DomainLayer<Mesh_> LayerType;
        typedef DomainLevel<Mesh_> LevelType;
        typedef Mesh_ MeshType;
        typedef Geometry::MeshAtlas<MeshType> AtlasType;

      public:
        explicit MeshStreamerDomainControl(Geometry::MeshFileReader& mesh_file_reader, int lvl_max, int lvl_min = -1) :
          BaseClass()
        {
          this->_create(mesh_file_reader, lvl_max, lvl_min);
        }


        explicit MeshStreamerDomainControl(const String& filename, int lvl_max, int lvl_min = -1) :
          BaseClass()
        {
          // create input file stream
          std::ifstream ifs(filename, std::ios_base::in);
          if(!ifs.is_open() || !ifs.good())
            throw InternalError(__func__,__FILE__,__LINE__,"Could not open file "+filename);

          Geometry::MeshFileReader mesh_file_reader(ifs);
          this->_create(mesh_file_reader, lvl_max, lvl_min);
        }

      protected:
        void _create(Geometry::MeshFileReader& mesh_file_reader, int lvl_max, int lvl_min)
        {
          typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

          // Communication data structures
          std::vector<Index> ranks, ctags;

          // Push layer
          this->_layers.push_back(new LayerType(std::move(ranks), std::move(ctags)));

          // Create mesh atlas
          this->_atlas = new AtlasType();
          // Create temporary root mesh node
          MeshNodeType* mesh_node = new MeshNodeType(nullptr, this->_atlas);

          // Parse mesh_node and atlas
          mesh_file_reader.parse(*mesh_node, *(this->_atlas));

          // adjust lvl_max and lvl_min
          lvl_max = Math::max(lvl_max, 0);
          if(lvl_min < 0)
            lvl_min = Math::max(lvl_max + lvl_min + 1, 0);
          else
            lvl_min = Math::min(lvl_min, lvl_max);

          // refine up to desired minimum level
          int lvl(0);
          for(; lvl < lvl_min; ++lvl)
          {
            MeshNodeType* coarse_node = mesh_node;
            mesh_node = coarse_node->refine();
            delete coarse_node;
          }

          // add coarse mesh node
          this->_levels.push_back(new LevelType(lvl_min, mesh_node));

          // refine up to desired maximum level
          for(; lvl < lvl_max; ++lvl)
          {
            MeshNodeType* coarse_node = mesh_node;
            mesh_node = coarse_node->refine();
            this->_levels.push_back(new LevelType(lvl+1, mesh_node));
          }
        }
      }; // class MeshStreamerDomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_MESH_STREAMER_DOMAIN_CONTROL_HPP

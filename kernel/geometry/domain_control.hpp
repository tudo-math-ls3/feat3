#pragma once
#ifndef KERNEL_GEOMETRY_DOMAIN_CONTROL_HPP
#define KERNEL_GEOMETRY_DOMAIN_CONTROL_HPP 1

#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/unit_cube_patch_generator.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>

#include <vector>
#include <deque>

namespace FEAST
{
  namespace Geometry
  {
    template<typename Mesh_>
    class DomainLayer
    {
    protected:
      std::vector<Index> _ranks;
      std::vector<Index> _ctags;

    public:
      explicit DomainLayer(std::vector<Index>&& ranks, std::vector<Index>&& ctags) :
        _ranks(std::forward<std::vector<Index>>(ranks)),
        _ctags(std::forward<std::vector<Index>>(ctags))
      {
        ASSERT_(_ranks.size() == _ctags.size());
      }

      virtual ~DomainLayer()
      {
      }

      Index size() const
      {
        return Index(_ranks.size());
      }

      Index get_rank(Index i) const
      {
        return _ranks.at(std::size_t(i));
      }

      Index get_ctag(Index i) const
      {
        return _ctags.at(std::size_t(i));
      }

      const std::vector<Index>& get_ranks() const
      {
        return _ranks;
      }

      const std::vector<Index>& get_ctags() const
      {
        return _ctags;
      }
    };

    template<typename Mesh_>
    class DomainLevel
    {
    public:
      typedef Mesh_ MeshType;
      typedef Geometry::MeshPart<Mesh_> PartType;
      typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

    protected:
      int _level_index;
      MeshNodeType* _mesh_node;

    public:
      explicit DomainLevel() :
        _level_index(0),
        _mesh_node(nullptr)
      {
      }

      explicit DomainLevel(int lvl_index, MeshNodeType* node) :
        _level_index(lvl_index),
        _mesh_node(node)
      {
      }

      virtual ~DomainLevel()
      {
        if(_mesh_node != nullptr)
        {
          delete _mesh_node;
        }
      }

      int get_level_index() const
      {
        return _level_index;
      }

      MeshNodeType* get_mesh_node()
      {
        return _mesh_node;
      }

      const MeshNodeType* get_mesh_node() const
      {
        return _mesh_node;
      }

      MeshType& get_mesh()
      {
        return *_mesh_node->get_mesh();
      }

      const MeshType& get_mesh() const
      {
        return *_mesh_node->get_mesh();
      }

      const PartType* find_halo_part(Index rank) const
      {
        // try to fetch the corresponding mesh part node
        return this->_mesh_node->find_mesh_part(String("halo:") + stringify(rank));
      }
    };

    template<typename Mesh_>
    class DomainControl
    {
    public:
      typedef DomainLayer<Mesh_> LayerType;
      typedef DomainLevel<Mesh_> LevelType;
      typedef Mesh_ MeshType;
      typedef Geometry::MeshAtlas<MeshType> AtlasType;

    protected:
      AtlasType* _atlas;
      std::deque<LayerType*> _layers;
      std::deque<LevelType*> _levels;

    public:
      DomainControl() :
        _atlas(nullptr)
      {
      }

      virtual ~DomainControl()
      {
        // delete layers
        while(!_layers.empty())
        {
          if(_layers.back() != nullptr)
            delete _layers.back();
          _layers.pop_back();
        }
        // delete levels
        while(!_levels.empty())
        {
          if(_levels.back() != nullptr)
            delete _levels.back();
          _levels.pop_back();
        }
        // delete atlas
        if(_atlas != nullptr)
        {
          delete _atlas;
        }
      }

      Index num_levels() const
      {
        return Index(_levels.size());
      }

      Index num_layers() const
      {
        return Index(_layers.size());
      }

      int min_level_index() const
      {
        return _levels.front()->get_level_index();
      }

      int max_level_index() const
      {
        return _levels.back()->get_level_index();
      }

      std::deque<LevelType*> get_levels()
      {
        return _levels;
      }

      const std::deque<LevelType*> get_levels() const
      {
        return _levels;
      }

      std::deque<LayerType*> get_layers()
      {
        return _layers;
      }

      const std::deque<LayerType*> get_layers() const
      {
        return _layers;
      }

      LevelType& get_level(Index i)
      {
        return *_levels.at(std::size_t(i));
      }

      const LevelType& get_level(Index i) const
      {
        return *_levels.at(std::size_t(i));
      }

      LayerType& get_layer(Index i)
      {
        return *_layers.at(std::size_t(i));
      }

      const LayerType& get_layer(Index i) const
      {
        return *_layers.at(std::size_t(i));
      }

      AtlasType* get_atlas()
      {
        return _atlas;
      }

      const AtlasType* get_atlas() const
      {
        return _atlas;
      }
    }; // class DomainControl<...>


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
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_DOMAIN_CONTROL_HPP

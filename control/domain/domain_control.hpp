#pragma once
#ifndef CONTROL_DOMAIN_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_DOMAIN_CONTROL_HPP 1

#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_quality_heuristic.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>

#include <vector>
#include <deque>

namespace FEAT
{
  namespace Control
  {
    /**
     * \brief Domain control namespace
     */
    namespace Domain
    {
      template<typename Mesh_>
      class DomainLayer
      {
      protected:
        Dist::Comm* _comm;
        std::vector<Index> _ranks;
        std::vector<Index> _ctags;

      public:
        explicit DomainLayer(Dist::Comm* comm_, std::vector<Index>&& ranks, std::vector<Index>&& ctags) :
          _comm(comm_),
          _ranks(std::forward<std::vector<Index>>(ranks)),
          _ctags(std::forward<std::vector<Index>>(ctags))
        {
          XASSERT(_ranks.size() == _ctags.size());
        }

        virtual ~DomainLayer()
        {
        }

        std::size_t bytes() const
        {
          return (_ranks.size() + _ctags.size()) * sizeof(Index);
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

        const Dist::Comm* get_comm() const
        {
          return _comm;
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

        std::size_t bytes() const
        {
          return _mesh_node->bytes();
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
          return this->_mesh_node->find_mesh_part(String("_halo:") + stringify(rank));
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
        typedef typename MeshType::CoordType CoordType;

      protected:
        /// the main communicator
        Dist::Comm* _comm;

        AtlasType* _atlas;
        std::deque<LayerType*> _layers;
        std::deque<LevelType*> _levels;

      public:
        DomainControl() :
          _comm(nullptr),
          _atlas(nullptr)
        {
        }

        explicit DomainControl(Dist::Comm* comm_) :
          _comm(comm_),
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

        static String name()
        {
          return "DomainControl<"+MeshType::name()+">";
        }

        virtual void print() const
        {
          Index pad_width(30);
          String msg;

          msg = String("num_levels").pad_back(pad_width, '.') + String(": ") + stringify(num_levels());
          _comm->print(msg);

          const auto& my_mesh = _levels.back()->get_mesh();
          Index ncells(my_mesh.get_num_entities(MeshType::shape_dim));
          _comm->allreduce(&ncells, &ncells, std::size_t(1), Dist::op_sum);

          msg = String("Cells on level "+stringify(_levels.back()->get_level_index())).pad_back(pad_width, '.')
            + String(": ") + stringify(ncells);
          _comm->print(msg);
        }

        std::size_t bytes() const
        {
          std::size_t s(0);
          if(_atlas != nullptr)
            s += _atlas->bytes();
          for(auto it = _layers.begin(); it != _layers.end(); ++it)
            s += (*it)->bytes();
          for(auto it = _levels.begin(); it != _levels.end(); ++it)
            s += (*it)->bytes();
          return s;
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

        std::deque<LevelType*>& get_levels()
        {
          return _levels;
        }

        const std::deque<LevelType*>& get_levels() const
        {
          return _levels;
        }

        std::deque<LayerType*>& get_layers()
        {
          return _layers;
        }

        const std::deque<LayerType*>& get_layers() const
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

        /**
         * \brief Computes mesh quality heuristics
         *
         * \param[out] edge_angle
         * The worst angle between two edges. Keep in mind that in 3d, this can be >0 even for deteriorated cells.
         *
         * \param[out] qi_min
         * The minimum quality indicator over all cells.
         *
         * \param[out] qi_mean
         * The mean quality indicator overa all cells.
         *
         * \param[out] edge_angle_cellwise
         * For debugging or visualisation purposes, this can receive the worst edge angle for every cell.
         *
         * \param[out] qi_cellwise
         * For debugging or visualisation purposes, this can receive the quality indicator for every cell.
         *
         * \param[in] lvl_index
         * Index of the level to compute everything for. Defaults to the maximum level.
         *
         */
        void compute_mesh_quality(CoordType& edge_angle, CoordType& qi_min, CoordType& qi_mean,
        CoordType* edge_angle_cellwise = nullptr, CoordType* qi_cellwise = nullptr,
        int lvl_index = -1) const
        {
          // max_level_index cannot be called for the default argument, so we do it here
          if(lvl_index == -1)
          {
            lvl_index = max_level_index();
          }
          ASSERT(lvl_index >= min_level_index());
          ASSERT(lvl_index <= max_level_index());

          CoordType qi_sum(0);

          for(const auto& it: _levels)
          {
            if(it->get_level_index() == lvl_index)
            {
              const auto& my_mesh = it->get_mesh();

              Index ncells(my_mesh.get_num_entities(MeshType::shape_dim));
              _comm->allreduce(&ncells, &ncells, std::size_t(1), Dist::op_sum);

              Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(qi_min, qi_sum,
              my_mesh.template get_index_set<MeshType::shape_dim, 0>(), my_mesh.get_vertex_set(), qi_cellwise);

              _comm->allreduce(&qi_min, &qi_min, std::size_t(1), Dist::op_min);
              _comm->allreduce(&qi_sum, &qi_sum, std::size_t(1), Dist::op_sum);

              edge_angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
              my_mesh.template get_index_set<MeshType::shape_dim, 0>(), my_mesh.get_vertex_set(), edge_angle_cellwise);
              _comm->allreduce(&edge_angle, &edge_angle, std::size_t(1), Dist::op_min);

              qi_mean = qi_sum/CoordType(ncells);

              return;
            }
          }

          // We should never get to this point
          throw InternalError(__func__,__FILE__,__LINE__,
          "Could not find level with index "+stringify(lvl_index)+"!\n");

        }
      }; // class DomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_DOMAIN_CONTROL_HPP

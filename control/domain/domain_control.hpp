#pragma once
#ifndef CONTROL_DOMAIN_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_DOMAIN_CONTROL_HPP 1

#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_quality_heuristic.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>

#include <control/domain/domain_level.hpp>

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
      class DomainLayer
      {
      protected:
        /// the communicator for this level
        Dist::Comm _comm;

        /// the layer index
        int _layer_index;

        /// the ranks of the neighbour processes in this layer
        std::vector<int> _neighbour_ranks;
        /// the ranks of the child processes w.r.t. the previous layer
        std::vector<int> _child_ranks;
        /// the ranks of the parent processes w.r.t. the next layer
        std::vector<int> _parent_ranks;

      public:
        explicit DomainLayer(Dist::Comm&& comm_, int lyr_idx) :
          _comm(std::forward<Dist::Comm>(comm_)),
          _layer_index(lyr_idx),
          _neighbour_ranks(),
          _child_ranks(),
          _parent_ranks()
        {
        }

        DomainLayer(DomainLayer&& other) :
          _comm(std::forward<Dist::Comm>(other._comm)),
          _layer_index(other._layer_index),
          _neighbour_ranks(std::forward<std::vector<int>>(other._neighbour_ranks)),
          _child_ranks(std::forward<std::vector<int>>(other._child_ranks)),
          _parent_ranks(std::forward<std::vector<int>>(other._parent_ranks))
        {
        }

        virtual ~DomainLayer()
        {
        }

        std::size_t bytes() const
        {
          return (_neighbour_ranks.size() + _child_ranks.size() + _parent_ranks.size()) * sizeof(int);
        }

        const Dist::Comm& comm() const
        {
          return _comm;
        }

        const Dist::Comm* comm_ptr() const
        {
          return &_comm;
        }

        bool is_child() const
        {
          return !_parent_ranks.empty();
        }

        bool is_parent() const
        {
          return !_child_ranks.empty();
        }

        bool is_ghost() const
        {
          return is_child() && !is_parent();
        }

        int get_layer_index() const
        {
          return _layer_index;
        }

        Index neighbour_count() const
        {
          return Index(_neighbour_ranks.size());
        }

        void push_neighbour(int neighbour_rank_)
        {
          _neighbour_ranks.push_back(neighbour_rank_);
        }

        int neighbour_rank(Index i) const
        {
          return _neighbour_ranks.at(std::size_t(i));
        }

        void set_neighbour_ranks(const std::vector<int>& neighbours)
        {
          _neighbour_ranks = neighbours;
        }

        const std::vector<int>& get_neighbour_ranks() const
        {
          return _neighbour_ranks;
        }

        Index child_count() const
        {
          return Index(_child_ranks.size());
        }

        void push_child(int child_rank_)
        {
          _child_ranks.push_back(child_rank_);
        }

        int child_rank(Index i) const
        {
          return _child_ranks.at(std::size_t(i));
        }

        void set_child_ranks(const std::vector<int>& children)
        {
          _child_ranks = children;
        }

        const std::vector<int>& get_child_ranks() const
        {
          return _child_ranks;
        }

        Index parent_count() const
        {
          return Index(_parent_ranks.size());
        }

        void push_parent(int parent_rank_)
        {
          _parent_ranks.push_back(parent_rank_);
        }

        int parent_rank(Index i) const
        {
          return _parent_ranks.at(std::size_t(i));
        }

        void set_parent_ranks(const std::vector<int>& parents)
        {
          _parent_ranks = parents;
        }

        const std::vector<int>& get_parent_ranks() const
        {
          return _parent_ranks;
        }
      }; // class DomainLayer

      template<typename DomLvl_>
      class VirtualLevel
      {
      public:
        typedef DomainLayer LayerType;
        typedef DomLvl_ LevelType;

      protected:
        std::shared_ptr<LevelType> _level;
        std::shared_ptr<LevelType> _level_child;
        std::shared_ptr<LevelType> _level_parent;
        std::shared_ptr<LayerType> _layer;
        std::shared_ptr<LayerType> _layer_child;
        std::shared_ptr<LayerType> _layer_parent;

      public:
        explicit VirtualLevel(std::shared_ptr<LevelType> level_, std::shared_ptr<LayerType> layer_) :
          _level(level_),
          _level_child(),
          _level_parent(),
          _layer(layer_),
          _layer_child(),
          _layer_parent()
        {
        }

        explicit VirtualLevel(
          std::shared_ptr<LevelType> level_child_,
          std::shared_ptr<LayerType> layer_child_,
          std::shared_ptr<LevelType> level_parent_,
          std::shared_ptr<LayerType> layer_parent_) :
          _level(level_parent_ ? level_parent_ : level_child_),
          _level_child(level_child_),
          _level_parent(level_parent_),
          _layer(layer_parent_ ? layer_parent_ : layer_child_),
          _layer_child(layer_child_),
          _layer_parent(layer_parent_)
        {
          XASSERT(bool(_level_child));
          XASSERT(bool(_layer_child));
          XASSERT(_layer_child->is_child());

          if(_level_parent)
          {
            XASSERT(bool(_layer_parent));
            XASSERT(_layer_parent->is_parent());
            XASSERT(_level_child->get_level_index() == _level_parent->get_level_index());
          }
          else
          {
            XASSERT(!bool(_layer_parent));
          }
        }

        virtual ~VirtualLevel()
        {
        }

        bool is_child() const
        {
          return bool(_level_child);
        }

        bool is_parent() const
        {
          return bool(_level_parent);
        }

        bool is_ghost() const
        {
          return is_child() && !is_parent();
        }

        LevelType& operator*()
        {
          return *_level;
        }

        const LevelType& operator*() const
        {
          return *_level;
        }

        LevelType* operator->()
        {
          return _level.get();
        }

        const LevelType* operator->() const
        {
          return _level.get();
        }

        LevelType& level()
        {
          return *_level;
        }

        const LevelType& level() const
        {
          return *_level;
        }

        LayerType& layer()
        {
          return *_layer;
        }

        const LayerType& layer() const
        {
          return *_layer;
        }

        LevelType& level_c()
        {
          XASSERT(bool(_level_child));
          return *_level_child;
        }

        const LevelType& level_c() const
        {
          XASSERT(bool(_level_child));
          return *_level_child;
        }

        LayerType& layer_c()
        {
          XASSERT(bool(_layer_child));
          return *_layer_child;
        }

        const LayerType& layer_c() const
        {
          XASSERT(bool(_layer_child));
          return *_layer_child;
        }

        LevelType& level_p()
        {
          XASSERT(bool(_level_parent));
          return *_level_parent;
        }

        const LevelType& level_p() const
        {
          XASSERT(bool(_level_parent));
          return *_level_parent;
        }

        LayerType& layer_p()
        {
          XASSERT(bool(_layer_parent));
          return *_layer_parent;
        }

        const LayerType& layer_p() const
        {
          XASSERT(bool(_layer_parent));
          return *_layer_parent;
        }
      }; // class VirtualLevel<...>

      template<typename DomLvl_>
      class DomainControl
      {
      public:
        typedef DomLvl_ LevelType;
        typedef DomainLayer LayerType;
        typedef VirtualLevel<LevelType> VirtLevelType;

        typedef typename LevelType::MeshType MeshType;
        typedef Geometry::MeshAtlas<MeshType> AtlasType;
        typedef typename MeshType::CoordType CoordType;

      protected:
        const Dist::Comm& _comm;
        AtlasType _atlas;
        std::deque<std::shared_ptr<LayerType>> _layers;
        std::deque<std::deque<std::shared_ptr<LevelType>>> _layer_levels;
        std::deque<VirtLevelType> _virt_levels;
        std::size_t _virt_size;

      public:
        explicit DomainControl(const Dist::Comm& comm_) :
          _comm(comm_),
          _virt_size(0u)
        {
          _layers.push_back(std::make_shared<LayerType>(_comm.comm_dup(), 0));
          _layer_levels.resize(std::size_t(1));
        }

        virtual ~DomainControl()
        {
        }

        static String name()
        {
          return "DomainControl<"+MeshType::name()+">";
        }

        virtual void print() const
        {
          Index pad_width(30);
          String msg;

          XASSERT(!_layers.empty());

          msg = String("num_levels").pad_back(pad_width, '.') + String(": ") + stringify(size_physical());
          _comm.print(msg);

          const auto& my_mesh = front()->get_mesh();
          Index ncells(my_mesh.get_num_entities(MeshType::shape_dim));
          _comm.allreduce(&ncells, &ncells, std::size_t(1), Dist::op_sum);

          msg = String("Cells on level "+stringify(max_level_index()).pad_back(pad_width, '.'))
            + String(": ") + stringify(ncells);
          _comm.print(msg);
        }

        const Dist::Comm& comm() const
        {
          return _comm;
        }

        std::size_t bytes() const
        {
          std::size_t s = _atlas.bytes();;
          for(auto it = _layers.begin(); it != _layers.end(); ++it)
            s += (*it)->bytes();
          return s;
        }

        Index num_layers() const
        {
          return Index(_layers.size());
        }

        int min_level_index() const
        {
          return _virt_levels.back()->get_level_index();
        }

        int max_level_index() const
        {
          return _virt_levels.front()->get_level_index();
        }

        void push_level_front(int layer_index, std::shared_ptr<LevelType> level)
        {
          _layer_levels.at(std::size_t(layer_index)).push_front(level);
        }

        void push_level_back(int layer_index, std::shared_ptr<LevelType> level)
        {
          _layer_levels.at(std::size_t(layer_index)).push_back(level);
        }

        bool has_ghost() const
        {
          /// \todo does this really work?
          return _virt_size < _virt_levels.size();
          /*if(_virt_levels.empty())
            return false;
          else
            return _virt_levels.back().is_ghost();*/
        }

        std::size_t size_virtual() const
        {
          return _virt_size;
        }

        std::size_t size_physical() const
        {
          if(_virt_levels.empty())
            return std::size_t(0);
          else if(_virt_levels.back().is_ghost())
            return _virt_levels.size() - std::size_t(1);
          else
            return _virt_levels.size();
        }

        VirtLevelType& at(std::size_t i)
        {
          return _virt_levels.at(i);
        }

        const VirtLevelType& at(std::size_t i) const
        {
          return _virt_levels.at(i);
        }

        VirtLevelType& front()
        {
          return _virt_levels.front();
        }

        const VirtLevelType& front() const
        {
          return _virt_levels.front();
        }

        VirtLevelType& back()
        {
          return _virt_levels.back();
        }

        const VirtLevelType& back() const
        {
          return _virt_levels.back();
        }

        AtlasType& get_atlas()
        {
          return _atlas;
        }

        const AtlasType& get_atlas() const
        {
          return _atlas;
        }

        void compile_virtual_levels()
        {
          _virt_levels.clear();

          // push the finest level
          _virt_levels.push_back(VirtLevelType(_layer_levels.front().front(), _layers.front()));

          // loop over all layers, front to back
          for(std::size_t ilay(0); ilay < _layers.size(); ++ilay)
          {
            auto& layer = _layers.at(ilay);
            auto& laylevs = _layer_levels.at(ilay);

            // the front level has already been added
            if(laylevs.size() <= std::size_t(1))
              continue;

            // push all inner layer levels
            for(std::size_t ilev(1); (ilev+1) < laylevs.size(); ++ilev)
              _virt_levels.push_back(VirtLevelType(laylevs.at(ilev), layer));

            // push the last level
            if((ilay+1) < _layers.size())
            {
              // that's an overlapping virtual mesh level with two physical mesh levels
              _virt_levels.push_back(VirtLevelType(laylevs.back(), layer, _layer_levels.at(ilay+1).front(), _layers.at(ilay+1)));
            }
            else if(layer->is_child())
            {
              // that's a "ghost" virtual mesh level with only one physical mesh level on this process
              _virt_levels.push_back(VirtLevelType(laylevs.back(), layer, nullptr, nullptr));
            }
            else
            {
              // that's a standard virtual mesh level with exactly one physical mesh level
              _virt_levels.push_back(VirtLevelType(laylevs.back(), layer));
            }
          }

          // get virtual size of this process
          std::size_t my_virt_size = _virt_levels.size();

          // get global maximum of all virtual sizes
          _comm.allreduce(&my_virt_size, &_virt_size, std::size_t(1), Dist::op_max);
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
          XASSERT(lvl_index >= min_level_index());
          XASSERT(lvl_index <= max_level_index());
          XASSERT(!_layers.empty());

          CoordType qi_sum(0);

          for(const auto& it: _virt_levels)
          {
            if(it->get_level_index() == lvl_index)
            {
              const auto& my_mesh = it->get_mesh();

              Index ncells(my_mesh.get_num_entities(MeshType::shape_dim));
              _comm.allreduce(&ncells, &ncells, std::size_t(1), Dist::op_sum);

              Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(qi_min, qi_sum,
              my_mesh.template get_index_set<MeshType::shape_dim, 0>(), my_mesh.get_vertex_set(), qi_cellwise);

              _comm.allreduce(&qi_min, &qi_min, std::size_t(1), Dist::op_min);
              _comm.allreduce(&qi_sum, &qi_sum, std::size_t(1), Dist::op_sum);

              edge_angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
              my_mesh.template get_index_set<MeshType::shape_dim, 0>(), my_mesh.get_vertex_set(), edge_angle_cellwise);
              _comm.allreduce(&edge_angle, &edge_angle, std::size_t(1), Dist::op_min);

              qi_mean = qi_sum/CoordType(ncells);

              return;
            }
          }

          // We should never get to this point
          throw InternalError(__func__,__FILE__,__LINE__,
          "Could not find level with index "+stringify(lvl_index)+"!\n");

        }

        void dump_layers()
        {
          String msg;
          msg += "(" + stringify(_layers.size()) + "):";
          std::size_t nc = std::size_t(1);
          for(std::size_t i(0); i < _layers.size(); ++i)
          {
            const auto& lyr = *_layers.at(i);
            std::size_t np = std::size_t(Math::ilog10(lyr.comm().size()));
            if(i > std::size_t(0))
              msg += " |";
            msg += " [" + stringify(lyr.comm().rank()).pad_front(np) + "]";
            for(const auto& cr : lyr.get_child_ranks())
              msg += " " + stringify(cr).pad_front(nc);
            if(lyr.is_child())
            {
              msg += " {";
              for(const auto& pr : lyr.get_parent_ranks())
                msg += " " + stringify(pr).pad_front(np);
              msg += "}";
            }
            nc = np;
          }
          _comm.allprint(msg);
        }

        void dump_layer_levels()
        {
          String msg;
          for(auto it = _layer_levels.begin(); it != _layer_levels.end(); ++it)
          {
            if(it != _layer_levels.begin())
              msg += " |";
            for(auto jt = it->begin(); jt != it->end(); ++jt)
              msg += " " + stringify((*jt)->get_level_index());
          }
          _comm.allprint(msg);
        }

        void dump_virt_levels()
        {
          String msg;
          for(auto it = _virt_levels.begin(); it != _virt_levels.end(); ++it)
          {
            std::size_t np = std::size_t(Math::ilog10((*it).layer().comm().size()));
            if((*it).is_child())
              np = std::size_t(Math::ilog10((*it).layer_c().comm().size()));
            msg += " " + stringify((*it)->get_level_index());
            msg += "[" + stringify((*it).layer().comm().size()).pad_front(np) + "]";
            if((*it).is_child())
            {
              msg += " {";
              for(const auto& pr : (*it).layer_c().get_parent_ranks())
                msg += " " + stringify(pr).pad_front(np);
              if((*it).is_parent())
                msg += ":" + stringify((*it).layer_p().comm().rank()).pad_front(np) + "}";
              else
                msg += String(np+1, ' ') + "}";
            }
          }
          _comm.allprint(msg);
        }
      }; // class DomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_DOMAIN_CONTROL_HPP

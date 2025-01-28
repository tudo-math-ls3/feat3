// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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
      /**
       * \brief Domain Layer class
       *
       * This class implements a class for domain layers, i.e. objects which represent
       * a partitioning of the available process set in a domain decomposition based
       * parallel application.
       *
       * In the most basic case, a domain layer consists of two important parts:
       * - A Dist::Comm communicator which represents the set of all MPI processes
       *   that participate in the decomposition of a domain.
       * - A set of ranks that represent the neighbors of this process within
       *   the communicator.
       *
       * In a simple parallel application, there exists just a single domain layer,
       * whose communicator is the world communicator (i.e. the set of all available
       * processes). In more complex scenarios, an application may contain a list of
       * domain layers, where each domain layer represents a subset of the processes
       * contained in the previous domain layer, thus forming a "process hierarchy".
       *
       * In such a process hierarchy, a whole group of "siblings" in a particular
       * domain layer is reduced onto a single parent process in the next domain layer.
       *
       * In case of multiple domain layers forming a hierarchy, each layer (except
       * for the last one) also contains two additional objects:
       * - A Dist::Comm communicator which represents the set of all sibling processes
       *   which are reduced to a single parent in the next layer.
       * - The rank of parent process within the sibling communicator.
       *
       * \author Peter Zajac
       */
      class DomainLayer
      {
      protected:
        /// the layer index
        int _layer_index;
        /// the communicator for this layer
        Dist::Comm _comm;
        /// the ranks of the neighbor processes in this layer
        std::vector<int> _neighbor_ranks;

        /// the sibling communicator
        Dist::Comm _sibling_comm;
        /// the rank of the parent in the sibling comm
        int _parent_rank;

      public:
        explicit DomainLayer(Dist::Comm&& comm_, int lyr_idx) :
          _layer_index(lyr_idx),
          _comm(std::forward<Dist::Comm>(comm_)),
          _neighbor_ranks(),
          _sibling_comm(),
          _parent_rank(-1)
        {
        }

        DomainLayer(DomainLayer&& other) :
          _layer_index(other._layer_index),
          _comm(std::forward<Dist::Comm>(other._comm)),
          _neighbor_ranks(std::forward<std::vector<int>>(other._neighbor_ranks)),
          _sibling_comm(std::forward<Dist::Comm>(other._sibling_comm)),
          _parent_rank(other._parent_rank)
        {
        }

        virtual ~DomainLayer()
        {
        }

        void set_parent(Dist::Comm&& sibling_comm_, int parent_rank)
        {
          _sibling_comm = std::forward<Dist::Comm>(sibling_comm_);
          _parent_rank = parent_rank;
        }

        std::size_t bytes() const
        {
          return _neighbor_ranks.size() * sizeof(int);
        }

        const Dist::Comm& comm() const
        {
          return _comm;
        }

        const Dist::Comm* comm_ptr() const
        {
          return &_comm;
        }

        const Dist::Comm* sibling_comm_ptr() const
        {
          return &_sibling_comm;
        }

        bool is_child() const
        {
          return (_sibling_comm.size() > 0);
        }

        bool is_parent() const
        {
          return (_sibling_comm.rank() == _parent_rank);
        }

        bool is_ghost() const
        {
          return is_child() && !is_parent();
        }

        int get_layer_index() const
        {
          return _layer_index;
        }

        int get_parent_rank() const
        {
          return _parent_rank;
        }

        Index child_count() const
        {
          return Index(_sibling_comm.size());
        }

        Index neighbor_count() const
        {
          return Index(_neighbor_ranks.size());
        }

        void push_neighbor(int neighbor_rank_)
        {
          _neighbor_ranks.push_back(neighbor_rank_);
        }

        int neighbor_rank(Index i) const
        {
          return _neighbor_ranks.at(std::size_t(i));
        }

        void set_neighbor_ranks(const std::vector<int>& neighbors)
        {
          _neighbor_ranks = neighbors;
        }

        const std::vector<int>& get_neighbor_ranks() const
        {
          return _neighbor_ranks;
        }
      }; // class DomainLayer

      /**
       * \brief Virtual Level Lambda type enumeration
       *
       * This enumeration is used by the VirtualLevel::apply_lambda() function, see its
       * documentation for details.
       */
      enum class VirtualLevelLambda
      {
        /// indicates that the level is a normal level
        normal,
        /// indicates that the level is a parent level
        parent,
        /// indicates that the level is a child level
        child,
        /// indicates that the level is a base level
        base
      };

      /**
       * \brief Virtual Domain Level class template
       *
       * This class acts as a wrapper around the actual domain levels, which contain the
       * actual mesh, trafo and finite element space objects. The purpose of this class
       * it to correctly manage domain levels which lay on the interface between two
       * domain layers, as in this case, some special treatment is required for the
       * implementation of a recursively partitioned mesh hierarchy.
       *
       * \tparam DomLvl_
       * The class representing a particular domain level; this class should derive from the
       * FEAT::Control::Domain::DomainLevel class template.
       *
       * \author Peter Zajac
       */
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
        std::shared_ptr<LevelType> _level_base;
        std::shared_ptr<LayerType> _layer;
        std::shared_ptr<LayerType> _layer_child;
        std::shared_ptr<LayerType> _layer_parent;
        bool _has_base;

      public:
        explicit VirtualLevel(std::shared_ptr<LevelType> level_, std::shared_ptr<LayerType> layer_) :
          _level(level_),
          _level_child(),
          _level_parent(),
          _level_base(),
          _layer(layer_),
          _layer_child(),
          _layer_parent(),
          _has_base(false)
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
          _level_base(),
          _layer(layer_parent_ ? layer_parent_ : layer_child_),
          _layer_child(layer_child_),
          _layer_parent(layer_parent_),
          _has_base(false)
        {
          XASSERT(bool(_level_child));
          XASSERT(bool(_layer_child));
          XASSERT(_layer_child->is_child());

          if(_level_parent)
          {
            XASSERT(bool(_layer_parent));
            XASSERT(_layer_child->is_parent());
            XASSERT(_level_child->get_level_index() == _level_parent->get_level_index());
          }
          else
          {
            XASSERT(!bool(_layer_parent));
          }
        }

        VirtualLevel(VirtualLevel&&) = default;
        VirtualLevel& operator=(VirtualLevel&&) = default;
        VirtualLevel(const VirtualLevel&) = delete;
        VirtualLevel& operator=(const VirtualLevel&) = delete;
        virtual ~VirtualLevel() = default;

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

        bool has_base() const
        {
          return _has_base;
        }

        void set_base(std::shared_ptr<LevelType> level_base_)
        {
          _level_base = level_base_;
          _has_base = true;
        }

        void set_base()
        {
          _has_base = true;
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

        LevelType& level_b()
        {
          XASSERT(bool(_level_base));
          return *_level_base;
        }

        const LevelType& level_b() const
        {
          XASSERT(bool(_level_base));
          return *_level_base;
        }

        /**
         * \brief Applies a lambda expression onto every level object in this virtual level
         *
         * \param[in] lambda
         * The lambda expression that has to be applied onto each level object
         *
         * The lambda expression is expected to take 3 parameters, where the first one is
         * a LevelType& reference, the second one is a LayerType& reference, which represents
         * the layer corresponding to the level object, and the third parameter is a
         * VirtualLevelLambda enumeration value which specifies whether the lambda is being
         * applied onto a normal level, a parent level, a child level or a base level.
         */
        template<typename Lambda_>
        void apply_lambda(Lambda_&& lambda)
        {
          if(_level_parent)
            lambda(*_level_parent, *_layer_parent, VirtualLevelLambda::parent);
          if(_level_child)
            lambda(*_level_child, *_layer_child, VirtualLevelLambda::child);
          else
            lambda(*_level, *_layer, VirtualLevelLambda::normal);
          if(_level_base)
            lambda(*_level_base, *_layer, VirtualLevelLambda::base);
        }

        std::size_t bytes() const
        {
          std::size_t b(0u);
          if(_level_base)
            b += _level_base->bytes();
          if(is_parent())
            return b + this->_level_child->bytes() + this->_layer_parent->bytes();
          else if(is_child())
            return b + this->_level_child->bytes();
          else
            return b + this->_level->bytes();
        }
      }; // class VirtualLevel<...>

      /**
       * \brief Domain control base-class template
       *
       * \author Peter Zajac
       */
      template<typename DomLvl_>
      class DomainControl
      {
      public:
        typedef DomLvl_ LevelType;
        typedef DomainLayer LayerType;
        typedef VirtualLevel<LevelType> VirtLevelType;

        typedef typename LevelType::ShapeType ShapeType;
        typedef typename LevelType::MeshType MeshType;
        typedef typename LevelType::MeshNodeType MeshNodeType;
        typedef Geometry::MeshAtlas<MeshType> AtlasType;

      protected:
        const Dist::Comm& _comm;
        AtlasType _atlas;
        std::deque<std::shared_ptr<LayerType>> _layers;
        std::deque<std::shared_ptr<LevelType>> _base_levels;
        std::deque<std::deque<std::shared_ptr<LevelType>>> _layer_levels;
        std::deque<VirtLevelType> _virt_levels;
        std::size_t _num_global_layers;
        std::size_t _virt_size;

        /// keep base-mesh levels on root process?
        bool _keep_base_levels;

      public:
        explicit DomainControl(const Dist::Comm& comm_) :
          _comm(comm_),
          _num_global_layers(0u),
          _virt_size(0u),
          _keep_base_levels(false)
        {
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
          std::size_t s = _atlas.bytes();
          for(const auto& lyr : _layers)
            s += lyr->bytes();
          for(const auto& lyr_lvl : _layer_levels)
            for(const auto& lvl : lyr_lvl)
              s += lvl->bytes();
          return s;
        }

        std::size_t num_local_layers() const
        {
          return _layers.size();
        }

        std::size_t num_global_layers() const
        {
          return _num_global_layers;
        }

        void push_layer(std::shared_ptr<LayerType> layer)
        {
          // make sure that the layer index is valid
          if(!_layers.empty())
          {
            XASSERT((_layers.back()->get_layer_index()+1) == layer->get_layer_index());
          }
          _layers.push_back(layer);
          _layer_levels.push_back(std::deque<std::shared_ptr<LevelType>>());
        }

        LayerType& front_layer()
        {
          return *_layers.front();
        }

        const LayerType& front_layer() const
        {
          return *_layers.front();
        }

        LayerType& back_layer()
        {
          return *_layers.back();
        }

        const LayerType& back_layer() const
        {
          return *_layers.back();
        }

        int min_level_index() const
        {
          return _virt_levels.back()->get_level_index();
        }

        int max_level_index() const
        {
          return _virt_levels.front()->get_level_index();
        }

        int med_level_index() const
        {
          return (_layer_levels.size() < std::size_t(2) ? -1 : _layer_levels.front().back()->get_level_index());
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

        /**
         * \brief Instructs the domain controller to keep the base-mesh levels after partitioning
         *
         * Calling this function requests the domain controller to keep the base-mesh levels and
         * refine them even after the partitioning process. This functionality is required to
         * assemble a Global::Splitter object, which allows to serialize and deserialize vectors
         * independent of the partitioning used to create the vector.
         *
         * \attention
         * The functionality requested by this function is \b NOT scalable for large process
         * counts, i.e. the root process \b WILL run out of memory in large-scale simulations.
         * This functionality should only be used for moderate process counts.
         */
        void keep_base_levels()
        {
          this->_keep_base_levels = true;
        }

        void compile_virtual_levels()
        {
          // get number of layers of this process
          std::size_t my_num_layers = _layers.size();

          // get global maximum of all layer sizes
          _comm.allreduce(&my_num_layers, &_num_global_layers, std::size_t(1), Dist::op_max);

          // clear virtual levels
          _virt_levels.clear();

          // push the finest level
          _virt_levels.push_back(VirtLevelType(_layer_levels.front().front(), _layers.front()));

          // do we have a base-level? If so, then set the base for the
          if(_keep_base_levels)
          {
            if(!_base_levels.empty())
              _virt_levels.front().set_base(_base_levels.front());
            else
              _virt_levels.front().set_base();
          }

          // loop over all layers, front to back
          for(std::size_t ilay(0); ilay < _layers.size(); ++ilay)
          {
            // get the current layer
            std::shared_ptr<LayerType>& layer = _layers.at(ilay);

            // get the deque of the current layer's levels
            std::deque<std::shared_ptr<LevelType>>& laylevs = _layer_levels.at(ilay);

            // the front level has already been added
            if(laylevs.size() <= std::size_t(1))
              continue;

            // push all inner layer levels
            for(std::size_t ilev(1); (ilev+1) < laylevs.size(); ++ilev)
            {
              _virt_levels.push_back(VirtLevelType(laylevs.at(ilev), layer));

              // push base-levels? (only for front layer)
              if(_keep_base_levels && (ilay == std::size_t(0)))
              {
                if(!_base_levels.empty())
                  _virt_levels.back().set_base(_base_levels.at(ilev));
                else
                  _virt_levels.back().set_base();
              }
            }

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

        void add_trafo_mesh_part_charts()
        {
          for(auto& bl : _base_levels)
            bl->add_trafo_mesh_part_charts();
          for(auto& ll : _layer_levels)
            for(auto& l : ll)
              l->add_trafo_mesh_part_charts();
          /*for(auto it = _virt_levels.begin(); it != _virt_levels.end(); ++it)
          {
            if(it->is_parent())
              it->level_p().add_trafo_mesh_part_charts();
            if(it->is_child())
              it->level_c().add_trafo_mesh_part_charts();
            else
              it->level().add_trafo_mesh_part_charts();
          }*/
        }

        String dump_layers() const
        {
          String msg;
          msg += "(" + stringify(_layers.size()) + "):";
          for(std::size_t i(0); i < _layers.size(); ++i)
          {
            const auto& lyr = *_layers.at(i);
            std::size_t np = std::size_t(Math::ilog10(lyr.comm().size()));
            if(i > std::size_t(0))
              msg += " |";
            msg += " [" + stringify(lyr.comm().rank()).pad_front(np) + "]";
            if(lyr.is_child())
            {
              std::size_t ns = std::size_t(Math::ilog10(lyr.comm().size()));
              msg += " " + stringify(lyr.sibling_comm_ptr()->rank()).pad_front(ns);
              msg += " {" + stringify(lyr.get_parent_rank()).pad_front(np);
              msg += lyr.is_parent() ? "*}" : " }";
            }
          }
          return msg;
        }

        String dump_layer_levels() const
        {
          String msg;
          for(auto it = _layer_levels.begin(); it != _layer_levels.end(); ++it)
          {
            if(it != _layer_levels.begin())
              msg += " |";
            for(auto jt = it->begin(); jt != it->end(); ++jt)
              msg += " " + stringify((*jt)->get_level_index());
          }
          return msg;
        }

        String dump_virt_levels() const
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
              msg += " { " + stringify((*it).layer_c().get_parent_rank());
              if((*it).is_parent())
                msg += ":" + stringify((*it).layer_p().comm().rank()).pad_front(np) + "}";
              else
                msg += String(np+1, ' ') + "}";
            }
          }
          return msg;
        }
      }; // class DomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

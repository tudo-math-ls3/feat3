#pragma once
#ifndef CONTROL_DOMAIN_UNIT_CUBE_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_UNIT_CUBE_DOMAIN_CONTROL_HPP 1

#include <kernel/geometry/unit_cube_patch_generator.hpp>
#include <control/domain/domain_control.hpp>

#include <deque>
#include <vector>

namespace FEAT
{
  namespace Control
  {
    namespace Domain
    {
      template<typename DomainLevel_>
      class UnitCubeDomainControl :
        public DomainControl<DomainLevel_>
      {
      public:
        typedef DomainControl<DomainLevel_> BaseClass;
        typedef typename BaseClass::LayerType LayerType;
        typedef typename BaseClass::LevelType LevelType;
        typedef typename BaseClass::MeshType MeshType;
        typedef typename BaseClass::AtlasType AtlasType;
        typedef typename BaseClass::MeshNodeType MeshNodeType;

      public:
        explicit UnitCubeDomainControl(const Dist::Comm& comm_, int lvl_max, int lvl_min = -1) :
          BaseClass(comm_)
        {
          int rank = comm_.rank();
          int nprocs = comm_.size();

          this->_layers.push_back(std::make_shared<LayerType>(comm_.comm_dup(), 0));
          this->_layer_levels.resize(std::size_t(1));

          // communication data structures
          std::vector<int> ranks;

          // create root mesh node
          std::shared_ptr<MeshNodeType> mesh_node;
          int lvl = (int)Geometry::UnitCubePatchGenerator<MeshType>::create(rank, Math::max(nprocs,1), mesh_node, ranks);

          // set front layer neighbour ranks
          this->_layers.front()->set_neighbour_ranks(ranks);

          // adjust lvl_max and lvl_min
          lvl_max = Math::max(lvl_max, 0);
          if(lvl_min < 0)
            lvl_min = Math::max(lvl_max + lvl_min + 1, 0);
          else
            lvl_min = Math::min(lvl_min, lvl_max);

          // refine up to desired minimum level
          for(; lvl < lvl_min; ++lvl)
          {
            auto coarse_node = mesh_node;
            mesh_node = std::shared_ptr<MeshNodeType>(coarse_node->refine());
          }

          // add coarse mesh node
          this->_layer_levels.front().push_front(std::make_shared<LevelType>(lvl, mesh_node));

          // refine up to desired maximum level
          for(; lvl < lvl_max;)
          {
            auto coarse_node = mesh_node;
            mesh_node = std::shared_ptr<MeshNodeType>(coarse_node->refine());
            this->_layer_levels.front().push_front(std::make_shared<LevelType>(++lvl, mesh_node));
          }

          // compile the virtual level list
          this->compile_virtual_levels();
        }
      }; // class UnitCubeDomainControl<...>


      template<typename DomainLevel_>
      class HierarchUnitCubeDomainControl :
        public DomainControl<DomainLevel_>
      {
      public:
        typedef DomainControl<DomainLevel_> BaseClass;
        typedef typename BaseClass::LayerType LayerType;
        typedef typename BaseClass::LevelType LevelType;
        typedef typename BaseClass::MeshType MeshType;
        typedef typename BaseClass::AtlasType AtlasType;
        typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
        typedef Geometry::MeshPart<MeshType> MeshPartType;

        static_assert(MeshType::shape_dim == 2, "HierarchUnitCubeDomainControl works only for 2D meshes");

      public:
        explicit HierarchUnitCubeDomainControl(const Dist::Comm& comm_, const std::deque<int>& lvls) :
          BaseClass(comm_)
        {
          _create(lvls);
        }

        explicit HierarchUnitCubeDomainControl(const Dist::Comm& comm_, const std::vector<String>& lvls) :
          BaseClass(comm_)
        {
          std::deque<int> ilvls;
          XASSERT(!lvls.empty());

          ilvls.resize(lvls.size());
          for(std::size_t i(0); i < lvls.size(); ++i)
          {
            if(!lvls.at(i).parse(ilvls.at(i)))
            {
              comm_.print(std::cerr, "ERROR: failed to parse '" + lvls.at(i) + "' as level");
              FEAT::Runtime::abort();
            }
          }
          if(ilvls.size() < std::size_t(2))
            ilvls.push_back(0);

          _create(ilvls);
        }

      protected:
        static int _ilog4(int x)
        {
          int r = 0;
          for(; x > 1; ++r)
          {
            if((x % 4) > 0)
              return -1;
            x /= 4;
          }
          return r;
        }

        // converts a rank from 2-level to lexi ordering
        static int _2lvl2lexi(int rank, int size)
        {
          const int ls = _ilog4(size);
          int ix(0), iy(0);
          for(int i(0); i < ls; ++i)
          {
            ix |= (rank >> (i+0)) & (1 << i);
            iy |= (rank >> (i+1)) & (1 << i);
          }
          return iy*(1 << ls) + ix;
        }

        // converts a vector of ranks from lexi to 2-level ordering
        static void _lexi22lvl(std::vector<int>& ranks, std::map<int,int>& reranks, int size)
        {
          const int ls = _ilog4(size);
          for(auto& rank : ranks)
          {
            const int ix = rank % (1 << ls);
            const int iy = rank / (1 << ls);
            int rr = 0;
            for(int i(0); i < ls; ++i)
            {
              rr |= (ix & (1 << i)) << (i+0);
              rr |= (iy & (1 << i)) << (i+1);
            }
            reranks.insert(std::make_pair(rank, rr));
            rank = rr;
          }
        }

        void _create(const std::deque<int>& lvls)
        {
          const int log4n = _ilog4(this->_comm.size());
          XASSERTM(log4n >= 0, "number of processes must be a power of 4");

          XASSERT(!lvls.empty());

          // create layers for this process
          _create_layers(log4n, int(lvls.size())-1);

          this->_layer_levels.resize(this->_layers.size());

          // loop over all layers
          for(std::size_t i(0); i < this->_layers.size(); ++i)
          {
            auto& layer = *this->_layers.at(i);
            auto& laylevs = this->_layer_levels.at(i);

            const int csize = layer.comm().size();
            const int crank = _2lvl2lexi(layer.comm().rank(), csize);

            std::vector<int> ranks;
            std::map<int,int> reranks; // neighbour rank map: lexi -> 2lvl

            // create root mesh node
            std::shared_ptr<MeshNodeType> mesh_node;
            const int base_lvl = (int)Geometry::UnitCubePatchGenerator<MeshType>::create(crank, csize, mesh_node, ranks);

            // translate neighbour ranks from lexi to 2lvl
            _lexi22lvl(ranks, reranks, csize);

            // rename halos
            mesh_node->rename_halos(reranks);

            // set layer neighbours
            layer.set_neighbour_ranks(ranks);

            // create domain level
            laylevs.push_front(std::make_shared<LevelType>(base_lvl, mesh_node));

            // refine once and create child mesh parts
            if(i > std::size_t(0))
            {
              std::shared_ptr<MeshNodeType> ref_node = std::shared_ptr<MeshNodeType>(laylevs.front()->get_mesh_node()->refine());
              this->_create_child_meshparts(layer, ref_node);
              laylevs.push_front(std::make_shared<LevelType>(base_lvl+1, ref_node));
            }

            const int head_lvl = laylevs.front()->get_level_index();
            const int fin_lvl = lvls.at(i);
            const int crs_lvl = lvls.at(i+1);

            // refine up the desired fine level
            for(int l(head_lvl); l < fin_lvl; ++l)
            {
              std::shared_ptr<MeshNodeType> ref_node = std::shared_ptr<MeshNodeType>(laylevs.front()->get_mesh_node()->refine());
              laylevs.push_front(std::make_shared<LevelType>(l+1, ref_node));
            }

            // drop all levels below desired coarse level
            for(int l(base_lvl); l < crs_lvl; ++l)
            {
              laylevs.pop_back();
            }
          }

          // compile the virtual level list
          this->compile_virtual_levels();
        }

        void _create_layers(int /*log4n*/, int nlayers)
        {
          // create main layer
          this->_layers.push_back(std::make_shared<LayerType>(this->_comm.comm_dup(), 0));

          // the world comm layer is already pushed
          for(int i(0); (i+1) < nlayers; ++i)
          {
            // get child layer
            DomainLayer& child = *this->_layers.back();

            // get the child comm
            const Dist::Comm& comm_c = child.comm();

            // my rank in child comm
            const int rank = comm_c.rank();

            // create sibling comm
            Dist::Comm comm_s = comm_c.comm_create_range_incl(4, rank & ~3);
            XASSERT(comm_s.size() == 4);

            // set sibling comm in child layer
            child.set_parent(std::move(comm_s), 0); // rank 0 is parent

            // create parent comm
            Dist::Comm comm_p = comm_c.comm_create_range_incl(comm_c.size()/4, 0, 4);
            if(comm_p.is_null())
              break;

            // finally, add to our layer deque
            this->_layers.push_back(std::make_shared<DomainLayer>(std::move(comm_p), i+1));
          }
        }

        void _create_child_meshparts(const DomainLayer& /*layer*/, std::shared_ptr<MeshNodeType> mesh_node)
        {
          Index num_elems = mesh_node->get_mesh()->get_num_elements();
          XASSERT(num_elems == 4);

          const Index num_entities[] = {4, 4, 1};

          const Index ivert[4][4] =
          {
            { 0, 4, 6, 8},
            { 4, 1, 8, 7},
            { 6, 8, 2, 5},
            { 8, 7, 5, 3}
          };
          const Index iedge[4][4] =
          {
            { 0, 10, 4, 8},
            { 1, 11, 8, 6},
            { 10, 2, 5, 9},
            { 11, 3, 9, 7}
          };

          for(int i(0); i < 4; ++i)
          {
            MeshPartType* part = new MeshPartType(num_entities, false);

            // manually override target indices
            Index* iv = part->template get_target_set<0>().get_indices();
            Index* ie = part->template get_target_set<1>().get_indices();
            Index* iq = part->template get_target_set<2>().get_indices();
            for(int j(0); j < 4; ++j)
            {
              iv[j] = ivert[i][j];
              ie[j] = iedge[i][j];
            }
            iq[0] = Index(i);

            mesh_node->add_patch(i, part);
          }
        }
      }; // class HierarchUnitCubeDomainControl<...>


      template<typename DomainLevel_>
      class HierarchUnitCubeDomainControl2 :
        public DomainControl<DomainLevel_>
      {
      public:
        typedef DomainControl<DomainLevel_> BaseClass;
        typedef typename BaseClass::LayerType LayerType;
        typedef typename BaseClass::LevelType LevelType;
        typedef typename BaseClass::MeshType MeshType;
        typedef typename BaseClass::AtlasType AtlasType;
        typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
        typedef Geometry::MeshPart<MeshType> MeshPartType;

        static_assert(MeshType::shape_dim == 2, "HierarchUnitCubeDomainControl works only for 2D meshes");

      public:
        explicit HierarchUnitCubeDomainControl2(const Dist::Comm& comm_, const std::deque<int>& lvls, const std::deque<int>& lyrs) :
          BaseClass(comm_)
        {
          _create(lvls, lyrs);
        }

        explicit HierarchUnitCubeDomainControl2(const Dist::Comm& comm_, const std::vector<String>& lvls) :
          BaseClass(comm_)
        {
          std::deque<int> ilvls, ilyrs;
          XASSERT(!lvls.empty());

          ilvls.resize(lvls.size());
          ilyrs.resize(lvls.size(), 1);
          for(std::size_t i(0); i < lvls.size(); ++i)
          {
            std::deque<String> parts;
            lvls.at(i).split_by_string(parts, ":");
            if(!parts.front().parse(ilvls.at(i)))
            {
              comm_.print(std::cerr, "ERROR: failed to parse '" + lvls.at(i) + "' as level");
              FEAT::Runtime::abort();
            }
            if((parts.size() > std::size_t(1)) && !parts.back().parse(ilyrs.at(i)))
            {
              comm_.print(std::cerr, "ERROR: failed to parse '" + lvls.at(i) + "' as layer");
              FEAT::Runtime::abort();
            }
          }
          if(ilvls.size() < std::size_t(2))
          {
            ilvls.push_back(0);
            ilyrs.push_back(0);
          }

          _create(ilvls, ilyrs);
        }

        explicit HierarchUnitCubeDomainControl2(const Dist::Comm& comm_, const std::deque<String>& lvls) :
          BaseClass(comm_)
        {
          std::deque<int> ilvls, ilyrs;
          XASSERT(!lvls.empty());

          ilvls.resize(lvls.size());
          ilyrs.resize(lvls.size(), 1);
          for(std::size_t i(0); i < lvls.size(); ++i)
          {
            std::deque<String> parts;
            lvls.at(i).split_by_string(parts, ":");
            if(!parts.front().parse(ilvls.at(i)))
            {
              comm_.print(std::cerr, "ERROR: failed to parse '" + lvls.at(i) + "' as level");
              FEAT::Runtime::abort();
            }
            if((parts.size() > std::size_t(1)) && !parts.back().parse(ilyrs.at(i)))
            {
              comm_.print(std::cerr, "ERROR: failed to parse '" + lvls.at(i) + "' as layer");
              FEAT::Runtime::abort();
            }
          }
          if(ilvls.size() < std::size_t(2))
          {
            ilvls.push_back(0);
            ilyrs.push_back(0);
          }

          _create(ilvls, ilyrs);
        }

        std::deque<std::pair<int,int>> get_level_indices() const
        {
          std::deque<std::pair<int,int>> lvs;
          for(std::size_t i(0); i < this->_layer_levels.size(); ++i)
            lvs.push_back(
              std::make_pair(
                this->_layer_levels.at(i).front()->get_level_index(),
                this->_layers.at(i)->comm().size()));
          lvs.push_back(
            std::make_pair(
              this->_layer_levels.back().back()->get_level_index(),
              this->_layers.back()->comm().size()));
          return lvs;
        }

      protected:
        static int _ilog4(int x)
        {
          int r = 0;
          for(; x > 1; ++r)
          {
            if((x % 4) > 0)
              return -1;
            x /= 4;
          }
          return r;
        }

        // converts a rank from 2-level to lexi ordering
        static int _2lvl2lexi(int rank, int size)
        {
          const int ls = _ilog4(size);
          int ix(0), iy(0);
          for(int i(0); i < ls; ++i)
          {
            ix |= (rank >> (i+0)) & (1 << i);
            iy |= (rank >> (i+1)) & (1 << i);
          }
          return iy*(1 << ls) + ix;
        }

        // converts a vector of ranks from lexi to 2-level ordering
        static void _lexi22lvl(std::vector<int>& ranks, std::map<int,int>& reranks, int size)
        {
          const int ls = _ilog4(size);
          for(auto& rank : ranks)
          {
            const int ix = rank % (1 << ls);
            const int iy = rank / (1 << ls);
            int rr = 0;
            for(int i(0); i < ls; ++i)
            {
              rr |= (ix & (1 << i)) << (i+0);
              rr |= (iy & (1 << i)) << (i+1);
            }
            reranks.insert(std::make_pair(rank, rr));
            rank = rr;
          }
        }

        void _create(const std::deque<int>& lvls, const std::deque<int>& lyrs)
        {
          const int log4n = _ilog4(this->_comm.size());
          XASSERTM(log4n >= 0, "number of processes must be a power of 4");

          XASSERT(!lvls.empty());
          XASSERT(lvls.size() == lyrs.size());

          // create layers for this process
          _create_layers(lyrs);

          this->_layer_levels.resize(this->_layers.size());

          // loop over all layers
          for(std::size_t i(0); i < this->_layers.size(); ++i)
          {
            auto& layer = *this->_layers.at(i);
            auto& laylevs = this->_layer_levels.at(i);

            const int csize = layer.comm().size();
            const int crank = _2lvl2lexi(layer.comm().rank(), csize);

            std::vector<int> ranks;
            std::map<int,int> reranks; // neighbour rank map: lexi -> 2lvl

            // create root mesh node
            std::shared_ptr<MeshNodeType> mesh_node;
            const int base_lvl = (int)Geometry::UnitCubePatchGenerator<MeshType>::create(crank, csize, mesh_node, ranks);

            // translate neighbour ranks from lexi to 2lvl
            _lexi22lvl(ranks, reranks, csize);

            // rename halos
            mesh_node->rename_halos(reranks);

            // set layer neighbours
            layer.set_neighbour_ranks(ranks);

            // create domain level
            laylevs.push_front(std::make_shared<LevelType>(base_lvl, mesh_node));

            // refine and create child mesh parts
            if(i > std::size_t(0))
            {
              int children = this->_layers.at(i-1)->comm().size() / csize;
              int crs_lvl = base_lvl;
              for(int nel(1); nel < children; nel *= 4)
              {
                std::shared_ptr<MeshNodeType> ref_node = std::shared_ptr<MeshNodeType>(laylevs.front()->get_mesh_node()->refine());
                laylevs.push_front(std::make_shared<LevelType>(++crs_lvl, ref_node));
                //ref_node = std::shared_ptr<MeshNodeType>(ref_node->refine());
              }
              this->_create_child_meshparts(laylevs.front()->get_mesh_node_ptr());
            }

            const int head_lvl = laylevs.front()->get_level_index();
            const int fin_lvl = lvls.at(i);
            const int crs_lvl = lvls.at(i+1);

            // refine up the desired fine level
            for(int l(head_lvl); l < fin_lvl; ++l)
            {
              std::shared_ptr<MeshNodeType> ref_node = std::shared_ptr<MeshNodeType>(laylevs.front()->get_mesh_node()->refine());
              laylevs.push_front(std::make_shared<LevelType>(l+1, ref_node));
            }

            // drop all levels below desired coarse level
            for(int l(base_lvl); l < crs_lvl; ++l)
            {
              laylevs.pop_back();
            }
          }

          // compile the virtual level list
          this->compile_virtual_levels();
        }

        void _create_layers(const std::deque<int>& lyrs)
        {
          // create main layer
          this->_layers.push_back(std::make_shared<LayerType>(this->_comm.comm_dup(), 0));

          // the world comm layer is already pushed
          for(int ilay(1); (ilay+1) < int(lyrs.size()); ++ilay)
          {
            // get child layer
            DomainLayer& child = *this->_layers.back();

            // get the child comm
            const Dist::Comm& comm_c = child.comm();

            // my rank in child comm
            const int rank = comm_c.rank();

            // compute number of siblings
            const int nsib = Math::sqr(1 << lyrs.at(std::size_t(ilay)));
            XASSERT(nsib <= comm_c.size());

            // create sibling comm
            Dist::Comm comm_s = comm_c.comm_create_range_incl(nsib, rank & ~(nsib-1));
            XASSERT(comm_s.size() == nsib);

            // set sibling comm in child layer
            child.set_parent(std::move(comm_s), 0); // rank 0 is parent

            // create parent comm
            Dist::Comm comm_p = comm_c.comm_create_range_incl(comm_c.size()/nsib, 0, nsib);
            if(comm_p.is_null())
              break;

            // finally, add to our layer deque
            this->_layers.push_back(std::make_shared<DomainLayer>(std::move(comm_p), ilay));
          }
        }

        void _create_child_meshparts(std::shared_ptr<MeshNodeType> mesh_node)
        {
          const Index num_elems = mesh_node->get_mesh()->get_num_elements();

          const Index num_entities[] = {4, 4, 1};
          const auto& ivert = mesh_node->get_mesh()->template get_index_set<2,0>();
          const auto& iedge = mesh_node->get_mesh()->template get_index_set<2,1>();

          for(Index i(0); i < num_elems; ++i)
          {
            MeshPartType* part = new MeshPartType(num_entities, false);

            // manually override target indices
            Index* iv = part->template get_target_set<0>().get_indices();
            Index* ie = part->template get_target_set<1>().get_indices();
            Index* iq = part->template get_target_set<2>().get_indices();

            for(int j(0); j < 4; ++j)
            {
              iv[j] = ivert[i][j];
              ie[j] = iedge[i][j];
            }
            iq[0] = i;

            mesh_node->add_patch(int(i), part);
          }
        }
      }; // class HierarchUnitCubeDomainControl2<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_UNIT_CUBE_DOMAIN_CONTROL_HPP

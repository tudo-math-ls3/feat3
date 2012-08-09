#pragma once
#ifndef KERNEL_FOUNDATION_MESH_HH
#define KERNEL_FOUNDATION_MESH_HH 1

#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/mesh_error.hpp>
#include <kernel/foundation/attribute.hpp>
#include <iostream>
#include <cmath>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief Attribute registration operation wrapper class template
     *
     * See specialisations.
     *
     * \tparam MeshType_
     * type of mesh
     *
     *
     * \author Markus Geveler
     */
    class MeshAttributeRegistration
    {
      public:
        /**
         * \brief member function executes a registration of an attribute with a mesh
         *
         * \param[in] mesh
         * target mesh reference
         *
         * \param[in] polytope_level
         * polytope level associated with the attribute
         */
        template<typename MeshType_>
        static unsigned execute(MeshType_ & mesh, const unsigned polytope_level)
        {
          mesh._attribute_polytopelevel_relations.push_back(polytope_level);

          ++mesh._num_attributes;
          return mesh._num_attributes - 1;
        }
    };

    ///required tnumber of topologies for given number of spatial dimensions
    enum RequiredNumTopologies
    {
      rnt_1D = 2,
      rnt_2D = 4,
      rnt_3D = 6
    };

    ///polytope level identifiers
    enum PolytopeLevels
    {
      pl_vertex = 0,
      pl_edge,
      pl_face,
      pl_polyhedron
    };

    ///indices of the stored polytopelevel-to-polytopelevel topologies
    enum InternalPolytopeIndices
    {
      ipi_vertex_edge = 0,
      ipi_edge_vertex,
      ipi_vertex_face,
      ipi_face_vertex,
      ipi_vertex_polyhedron,
      ipi_polyhedron_vertex
    };

    ///indices of primary access when searching for adjacencies
    enum InternalPrimaryAccess
    {
      ipa_vertex_vertex = 0,
      ipa_vertex_edge = 0,
      ipa_vertex_face = 2,
      ipa_vertex_polyhedron = 4,
      ipa_edge_any = 1,
      ipa_face_any = 3,
      ipa_polyhedron_any = 5
    };

    ///indices of secondary access
    enum InternalSecondaryAccess
    {
      isa_vertex_vertex = 1,
      isa_any_edge = 0,
      isa_any_face = 2,
      isa_any_polyhedron = 4
    };

    /**
     * \brief Mesh class template is the base mesh type for all top-level FEAST meshes
     *
     * Mesh relies on a Render Dynamic Mesh technique that is modified to grant explicit (pseudo-constant) time
     * access from every polytope level to any other polytope level without storing all topologies explicitly.
     * Mesh builds a ring of topologies each representing adjacencies from polytope-level k to level k-1 or k+1.
     *
     * \tparam i_
     * required number of topologies (current policy implements the ring but is generally dimension-independent)
     *
     * \tparam TopologyType_
     * type of topology
     *
     * \author Markus Geveler
     */
    template<
      RequiredNumTopologies i_ = rnt_2D,
      typename TopologyType_ = Topology<>,
      template <typename, typename> class OuterStorageType_ = std::vector,
      template <typename, typename> class AttributeStorageType_ = std::vector,
      template <typename, typename> class OuterAttributeStorageType_ = std::vector
      >
    class Mesh
    {
      public:
        friend class MeshAttributeRegistration;

        ///type exports
        typedef TopologyType_ topology_type_;
        typedef typename TopologyType_::index_type_ index_type_;
        typedef typename TopologyType_::storage_type_ storage_type_;

        typedef OuterAttributeStorageType_<
          AttributeBase<AttributeStorageType_>*, std::allocator<AttributeBase<AttributeStorageType_>* > > attr_base_type_;

        ///CTOR
        Mesh(const typename TopologyType_::index_type_ id, attr_base_type_* attrbase = nullptr, const typename TopologyType_::index_type_ pp_rank = 0, const typename TopologyType_::index_type_ mp_rank = 0) :
          _id(id),
          _pp_rank(pp_rank),
          _mp_rank(mp_rank),
          _num_inter_topologies(i_),
          _num_levels((unsigned)(i_/2u) + 1u),
          _topologies(OuterStorageType_<TopologyType_, std::allocator<TopologyType_> >()),
          _attrs(attrbase),
          _num_attributes(0),
          _attribute_polytopelevel_relations(typename TopologyType_::storage_type_())
      {
        for(Index i(0) ; i < i_ ; ++i)
        {
          TopologyType_ t;
          _topologies.push_back(t);
        }
      }

        ///Copy CTOR
        Mesh(const typename TopologyType_::index_type_ new_id, Mesh & other, attr_base_type_* attrbase = nullptr) :
          _id(new_id),
          _pp_rank(other._pp_rank),
          _mp_rank(other._mp_rank),
          _num_inter_topologies(other._num_inter_topologies),
          _num_levels(other._num_levels),
          _attrs(attrbase),
          _num_attributes(0)
      {
        for(Index i(0) ; i < _num_inter_topologies ; ++i)
        {
          _topologies.push_back(other._topologies.at(i));
        }

        _attribute_polytopelevel_relations = other._attribute_polytopelevel_relations;
      }

        ///DTOR
        ~Mesh()
        {
        }

        void add_polytope(const PolytopeLevels level)
        {
          CONTEXT("When adding polytope");
#ifdef FOUNDATION_DEBUG
          if(level > _num_levels - 1)
            throw MeshInternalIndexOutOfBounds(level, _num_levels - 1);
#endif
          switch(level)
          {
            case pl_vertex:
              {
                if(_num_levels > pl_edge)
                  _topologies.at(ipi_vertex_edge).push_back();
                if(_num_levels > pl_face)
                  _topologies.at(ipi_vertex_face).push_back();
                if(_num_levels > pl_polyhedron)
                  _topologies.at(ipi_vertex_polyhedron).push_back();
              }
              break;

            case pl_edge:
              {
                _topologies.at(ipi_edge_vertex).push_back();
              }
              break;

            case pl_face:
              {
                _topologies.at(ipi_face_vertex).push_back();
              }
              break;

            case pl_polyhedron:
              {
                _topologies.at(ipi_polyhedron_vertex).push_back();
              }
              break;
          }

        }

        ///Add an adjacency to the associated topology and automatically add it to its transpose
        void add_adjacency(
            const PolytopeLevels from_level,
            const PolytopeLevels to_level,
            const index_type_ polytope_index,
            const index_type_ value)
        {
          CONTEXT("When adding adjacency");
#ifdef FOUNDATION_DEBUG
          if(from_level > _num_levels - 1)
            throw MeshInternalIndexOutOfBounds(from_level, _num_levels - 1);

          if(to_level > _num_levels - 1)
            throw MeshInternalIndexOutOfBounds(to_level, _num_levels - 1);

          if(from_level == to_level)
            throw MeshError("Explicit neighbourhood of polytopes on same level currently unsupported. Use implict storage instead.");//TODO
#endif

          InternalPolytopeIndices ipi(_get_internal_index(from_level, to_level));
          InternalPolytopeIndices ipit(_transpose_internal_index(ipi));

          _topologies.at(ipi).at(polytope_index).push_back(value);
          _topologies.at(ipit).at(value).push_back(polytope_index);

        }

        ///get all adjacent polytopes at level to_level for polytope i on level from_level
        typename TopologyType_::storage_type_ get_adjacent_polytopes(
                                                                     const PolytopeLevels from_level,
                                                                     const PolytopeLevels to_level,
                                                                     index_type_ i
                                                                    )
        {
          CONTEXT("When calculating adjacent polytopes");

#ifdef FOUNDATION_DEBUG
          if(from_level > _num_levels - 1)
            throw MeshInternalIndexOutOfBounds(from_level, _num_levels - 1);
          if(to_level > _num_levels - 1)
            throw MeshInternalIndexOutOfBounds(to_level, _num_levels - 1);
#endif
          InternalPrimaryAccess ipa(_get_primary_index(from_level, to_level));

          std::cout << "from: " << from_level << std::endl;
          std::cout << "to: " << to_level << std::endl;
          std::cout << "IPA: " << ipa << std::endl;

          if(!_secondary_access_needed(from_level, to_level))
            return _topologies.at(ipa).at(i);
          else
          {
            typename TopologyType_::storage_type_ result;

            InternalSecondaryAccess isa(_get_secondary_index(from_level, to_level));
            std::cout << "ISA: " << isa << std::endl;

            std::cout << "i: " << i << std::endl;
            std::cout << "_topologies.at(ipa).at(i).size() = " << _topologies.at(ipa).at(i).size() << std::endl;
            for(index_type_ j(0) ; j < _topologies.at(ipa).at(i).size() ; ++j)
            {
              std::cout << "j: " << j << std::endl;
              std::cout << "_topologies.at(isa).at(j).size() = " << _topologies.at(isa).at(j).size() << std::endl;
              for(index_type_ k(0) ; k < _topologies.at(isa).at(_topologies.at(ipa).at(i).at(j)).size() ; ++k)
              {
                index_type_ to_insert(_topologies.at(isa).at(_topologies.at(ipa).at(i).at(j)).at(k));
                std::cout << "Candidate: " << to_insert << std::endl;

                bool insert(true);
                for(index_type_ search_index(0) ; search_index < result.size(); ++search_index)
                {
                  if(result.at(search_index) == to_insert)
                  {
                    insert = false;
                    break;
                  }
                }
                if(insert)
                {
                  result.push_back(to_insert);
                }
              }
            }
            return result;
          }

          /*
          const int sweepdir(_sweep_direction(from_level, to_level));
          int current_sweepdir(sweepdir == 0 ? -1 : sweepdir);
          current_sweepdir = from_level == 0 && sweepdir == 0 ? 1 : current_sweepdir;
          const unsigned level_diff(_level_difference(from_level, to_level));
          const unsigned path_length(level_diff + 1);
          unsigned current_level(from_level);
          typename TopologyType_::compound_storage_type_ search_data;


          if(level_diff == 0u) //we can directly give the neighbours
          {
            if(sweepdir == -1) //sweep down
              return _topologies.at(_downward_index(from_level)).at(i);
            else if(sweepdir == 1) //sweep up
              return _topologies.at(_upward_index(from_level)).at(i);
#ifdef FOUNDATION_DEBUG
            else
              throw MeshError("Sweep direction / sweep length mismatch.");
#endif
          }
          else
          {

            for(unsigned j(0) ; j < path_length ; ++j)
            {
              //create current search datastructure
              typename TopologyType_::storage_type_ temp;
              search_data.push_back(temp);

              //fill current datastructure:
              if(j == 0)
              {
                search_data.at(0) = current_sweepdir == -1 ? _topologies.at(_downward_index(current_level)).at(i)
                  : _topologies.at(_upward_index(current_level)).at(i);

                current_level = current_sweepdir == 1 ? current_level + 1
                  : (current_sweepdir == 0 ? current_level + 1
                      : current_level - 1);


                current_sweepdir = sweepdir == 0 && current_sweepdir == -1 ? 1
                  : (sweepdir == 0 && current_sweepdir == 1 ? -1
                      : current_sweepdir);
              }
              else
              {
                //for all entries in previous polytope list get all polytope lists and store sequentially
                for(IndexType_ k(0) ; k < (IndexType_)search_data.at(j - 1).size(); ++k)
                {
                  IndexType_ l_upper(
                      current_sweepdir == -1 ? (IndexType_)_topologies.at(_downward_index(current_level)).at(search_data.at(j - 1).at(k)).size()
                      : (IndexType_)_topologies.at(_upward_index(current_level)).at(search_data.at(j - 1).at(k)).size()
                      );

                  for(IndexType_ l(0) ; l < l_upper ; ++l)
                  {
                    //TODO optimise search
                    IndexType_ to_insert(
                        current_sweepdir == -1 ? _topologies.at(_downward_index(current_level)).at(search_data.at(j - 1).at(k)).at(l)
                        : _topologies.at(_upward_index(current_level)).at(search_data.at(j - 1).at(k)).at(l)
                        );

                    bool insert(true);
                    for(IndexType_ search_index(0) ; search_index < (IndexType_)search_data.at(j).size() ; ++search_index)
                    {
                      if((IndexType_)search_data.at(j).at(search_index) == (IndexType_)to_insert)
                      {
                        insert = false;
                        break;
                      }
                    }
                    if(insert)
                    {
                      search_data.at(j).push_back(to_insert);
                    }
                  }
                }
                current_level = current_sweepdir == 1 ? current_level + 1
                  : current_level - 1;
              }
            }
            return search_data.at(search_data.size() - 1);
          }
          return search_data.at(search_data.size() - 1);
          */
        }

        ///needed public access functions
        typename TopologyType_::index_type_ get_num_levels()
        {
          return _num_levels;
        }

        int get_downward_index(const unsigned pl)
        {
          return _downward_index(pl);
        }

        int get_upward_index(const unsigned pl)
        {
          return _upward_index(pl);
        }

        ///get primary communication neighbours for element polytope i (that is only non-diagonal)
        template<typename IndexType_>
          typename TopologyType_::storage_type_ get_primary_comm_neighbours(IndexType_ i)
          {
            typename TopologyType_::storage_type_ pre_result(this->get_adjacent_polytopes(_num_levels - 1, _num_levels - 1, i));
            typename TopologyType_::storage_type_ result;

            //remove self
            for(IndexType_ j(0) ; j < (IndexType_)pre_result.size() ; ++j)
              if((IndexType_)pre_result.at(j) != i)
                result.push_back(pre_result.at(j));

            return result;
          }

        ///get all communication neighbours for element polytope i
        template<typename IndexType_>
          typename TopologyType_::storage_type_ get_all_comm_neighbours(IndexType_ i)
          {
            const int sweepdir(-1);
            int current_sweepdir(sweepdir);
            const unsigned path_length(_num_inter_topologies);
            unsigned current_level(_num_levels - 1);
            typename TopologyType_::compound_storage_type_ search_data;

            for(unsigned j(0) ; j < path_length ; ++j)
            {
              //create current search datastructure
              typename TopologyType_::storage_type_ temp;
              search_data.push_back(temp);

              //fill current datastructure:
              if(j == 0)
              {
                search_data.at(0) = _topologies.at(_downward_index(current_level)).at(i);

                --current_level;

                current_sweepdir = current_level == 0 ? 1 : current_sweepdir;
              }
              else
              {
                //for all entries in previous polytope list get all polytope lists and store sequentially
                for(IndexType_ k(0) ; k < (IndexType_)search_data.at(j - 1).size(); ++k)
                {
                  IndexType_ l_upper(
                      current_sweepdir == -1 ? (IndexType_)_topologies.at(_downward_index(current_level)).at(search_data.at(j - 1).at(k)).size()
                      : (IndexType_)_topologies.at(_upward_index(current_level)).at(search_data.at(j - 1).at(k)).size()
                      );

                  for(IndexType_ l(0) ; l < l_upper ; ++l)
                  {
                    //TODO optimise search
                    IndexType_ to_insert(
                        current_sweepdir == -1 ? _topologies.at(_downward_index(current_level)).at(search_data.at(j - 1).at(k)).at(l)
                        : _topologies.at(_upward_index(current_level)).at(search_data.at(j - 1).at(k)).at(l)
                        );

                    bool insert(true);
                    for(IndexType_ search_index(0) ; search_index < (IndexType_)search_data.at(j).size() ; ++search_index)
                    {
                      if((IndexType_)search_data.at(j).at(search_index) == (IndexType_)to_insert)
                      {
                        insert = false;
                        break;
                      }
                    }
                    if(insert)
                    {
                      search_data.at(j).push_back(to_insert);
                    }
                  }
                }
                current_level = current_sweepdir == 1 ? current_level + 1
                  : current_level - 1;

                current_sweepdir = current_level == 0 ? 1 : current_sweepdir;
              }
            }

            typename TopologyType_::storage_type_ result;

            //remove self
            for(IndexType_ j(0) ; j < (IndexType_)search_data.at(search_data.size() - 1).size() ; ++j)
              if((IndexType_)search_data.at(search_data.size() - 1).at(j) != i)
                result.push_back(search_data.at(search_data.size() - 1).at(j));

            return result;
          }

        ///Further needed access functions
        unsigned get_num_attributes()
        {
          return _num_attributes;
        }


        typename TopologyType_::storage_type_ & get_attribute_polytopelevel_relations()
        {
          return _attribute_polytopelevel_relations;
        }

        typename TopologyType_::index_type_ get_id()
        {
          return _id;
        }

        typename TopologyType_::index_type_ get_pp_rank()
        {
          return _pp_rank;
        }

        typename TopologyType_::index_type_ get_mp_rank()
        {
          return _mp_rank;
        }

        void set_pp_rank(typename TopologyType_::index_type_ rank)
        {
          _pp_rank = rank;
        }

        void set_mp_rank(typename TopologyType_::index_type_ rank)
        {
          _mp_rank = rank;
        }

        attr_base_type_* get_attributes()
        {
          return _attrs;
        }

      private:
        const typename TopologyType_::index_type_ _id;
        const typename TopologyType_::index_type_ _pp_rank;
        const typename TopologyType_::index_type_ _mp_rank;
        const unsigned _num_inter_topologies;
        const unsigned _num_levels;
        OuterStorageType_<TopologyType_, std::allocator<TopologyType_> > _topologies;

        attr_base_type_* _attrs;

        unsigned _num_attributes;

        typename TopologyType_::storage_type_ _attribute_polytopelevel_relations;

        inline unsigned _level_difference(const unsigned from, const unsigned to)
        {
          return from == to ? 1u : ( from > to ? (unsigned)std::abs((double)(from - to)) - 1u : (unsigned)std::abs((double)(to - from)) - 1u);
        }

        inline int _downward_index(const unsigned pl)
        {
          return (pl == 0u || pl >= _num_levels) ? - 1 : (pl == 3u ? 5 : (1 << pl) - 1);
        }

        inline int _upward_index(const unsigned pl)
        {
          return pl >= _num_levels - 1? -1 : ( pl > 0 ? (1 << pl) : 0);
        }

        inline int _sweep_direction(const unsigned from_level, const unsigned to_level)
        {
          return from_level == to_level ? 0 : (from_level > to_level ? -1 : 1);
        }

        InternalPrimaryAccess _get_primary_index(PolytopeLevels from, PolytopeLevels to)
        {
          switch(from)
          {
            case pl_vertex:
              {
                switch(to)
                {
                  case pl_vertex:
                    {
                      return ipa_vertex_vertex;
                    }
                    break;

                  case pl_edge:
                    {
                      return ipa_vertex_edge;
                    }
                    break;

                  case pl_face:
                    {
                      return ipa_vertex_face;
                    }
                    break;

                  case pl_polyhedron:
                    {
                      return ipa_vertex_polyhedron;
                    }
                    break;
                }
              }
              break;

            case pl_edge:
              {
                return ipa_edge_any;
              }
              break;

            case pl_face:
              {
                return ipa_face_any;
              }
              break;

            case pl_polyhedron:
              {
                return ipa_polyhedron_any;
              }
          }
        }

        InternalSecondaryAccess _get_secondary_index(PolytopeLevels from, PolytopeLevels to)
        {
          return from == pl_vertex && to == pl_vertex ? InternalSecondaryAccess(1) : InternalSecondaryAccess((to - 1) * 2);
        }

        InternalPolytopeIndices _get_internal_index(PolytopeLevels from, PolytopeLevels to)
        {
          switch(from)
          {
            case pl_vertex:
              {
                switch(to)
                {
                  case pl_edge:
                    {
                      return ipi_vertex_edge;
                    }
                    break;

                  case pl_face:
                    {
                      return ipi_vertex_face;
                    }
                    break;

                  case pl_polyhedron:
                    {
                      return ipi_vertex_polyhedron;
                    }
                    break;

                  default:
                    throw MeshError("Invalid impairment of polytope levels! Polytope list does not exist. Note - stored is: (vertex,.), (.,vertex) ..." );
                }
              }
              break;

            case pl_edge:
              {
                switch(to)
                {
                  case pl_vertex:
                    {
                      return ipi_edge_vertex;
                    }
                    break;

                  default:
                    throw MeshError("Invalid impairment of polytope levels! Polytope list does not exist. Note - stored is: (vertex,.), (.,vertex) ..." );
                }
              }
              break;

            case pl_face:
              {
                switch(to)
                {
                  case pl_vertex:
                    {
                      return ipi_face_vertex;
                    }
                    break;

                  default:
                    throw MeshError("Invalid impairment of polytope levels! Polytope list does not exist. Note - stored is: (vertex,.), (.,vertex) ..." );
                }
              }
              break;

            case pl_polyhedron:
              {
                switch(to)
                {
                  case pl_vertex:
                    {
                      return ipi_polyhedron_vertex;
                    }
                    break;

                  default:
                    throw MeshError("Invalid impairment of polytope levels! Polytope list does not exist. Note - stored is: (vertex,.), (.,vertex) ..." );
                }
              }
              break;

            default:
              throw MeshError("Invalid impairment of polytope levels! Polytope list does not exist. Note - stored is: (vertex,.), (.,vertex) ..." );
          }
        }

        bool _secondary_access_needed(PolytopeLevels from, PolytopeLevels to)
        {
          return (from != pl_vertex && to != pl_vertex) || (from == pl_vertex && to == from);
        }

        InternalPolytopeIndices _transpose_internal_index(InternalPolytopeIndices ipi)
        {
          return (InternalPolytopeIndices)(ipi == ipi_vertex_edge ? ipi_edge_vertex : (ipi % 2 == 0 ? ipi + 1 : ipi - 1));
        }
    };
  }
}
#endif

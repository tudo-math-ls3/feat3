#pragma once
#ifndef KERNEL_FOUNDATION_MESH_HPP
#define KERNEL_FOUNDATION_MESH_HPP 1

#include <kernel/foundation/base.hpp>
#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/mesh_error.hpp>
#include <kernel/foundation/buffer.hpp>
#include <kernel/foundation/communication.hpp>
#include <iostream>
#include <cmath>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief Indices of the stored polytopelevel-to-polytopelevel topologies
     *
     * \author Markus Geveler
     */
    enum InternalPolytopeIndices
    {
      ipi_vertex_edge = 0,
      ipi_edge_vertex,
      ipi_vertex_face,
      ipi_face_vertex,
      ipi_vertex_polyhedron,
      ipi_polyhedron_vertex
    };

    /**
     * \brief Indices of primary access when searching for adjacencies
     *
     * \author Markus Geveler
     */
    enum InternalPrimaryAccess
    {
      ipa_vertex_vertex = 0,
      ipa_vertex_edge = 0,
      ipa_vertex_face = 2,
      ipa_vertex_polyhedron = 4,
      ipa_edge_any = 1,
      ipa_face_any = 3,
      ipa_polyhedron_any = 5,
      ipa_none = 6
    };

    /**
     * \brief Indices of secondary access
     *
     * \author Markus Geveler
     */
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
     *
     * \tparam i_
     * required number of topologies
     *
     * \tparam TopologyType_
     * type of topology
     *
     * \tparam OuterStorageType_
     * STL conformal container template type storing the topologies
     *
     * \author Markus Geveler
     */
    template<
      typename Dim_ = Dim2D,
      typename TopologyType_ = Topology<>,
      template <typename, typename> class OuterStorageType_ = std::vector
      >
    class Mesh :
      public CommunicateableByAggregates<TopologyType_, com_send_receive>
    {
      public:
        ///type exports
        typedef TopologyType_ topology_type_;
        typedef typename TopologyType_::index_type_ index_type_;
        typedef typename TopologyType_::storage_type_ storage_type_;
        typedef typename TopologyType_::buffer_type_ buffer_type_;

        ///CTOR
        Mesh(const typename TopologyType_::index_type_ id = 0, const typename TopologyType_::index_type_ pp_rank = 0, const typename TopologyType_::index_type_ mp_rank = 0) :
          _id(id),
          _pp_rank(pp_rank),
          _mp_rank(mp_rank),
          _num_inter_topologies(Dim_::required_num_topologies),
          _num_levels((unsigned)(Dim_::required_num_topologies/2u) + 1u),
          _topologies(Dim_::required_num_topologies)
        {
        }

        ///Copy CTOR
        Mesh(const typename TopologyType_::index_type_ new_id, Mesh& other) :
          _id(new_id),
          _pp_rank(other._pp_rank),
          _mp_rank(other._mp_rank),
          _num_inter_topologies(other._num_inter_topologies),
          _num_levels(other._num_levels),
          _topologies(other._topologies)
        {
        }

        ///assignment operator overload needed for copy-assignability
        Mesh& operator=(const Mesh& rhs)
        {
          if(this == &rhs)
            return *this;

          this-> _id = rhs._id;
          this-> _pp_rank = rhs._pp_rank;
          this-> _mp_rank = rhs._mp_rank;
          this-> _mp_rank = rhs._mp_rank;
          this-> _num_inter_topologies = rhs._num_inter_topologies;
          this-> _num_levels = rhs._num_levels;
          this-> _topologies = rhs._topologies;

          return *this;
        }

        ///DTOR
        ~Mesh()
        {
        }

        /**
         * \brief adds a polytope to the mesh
         *
         * \param[in] level
         * level to add a polytope to anonymously
         *
         */
        void add_polytope(const PolytopeLevels level)
        {
          CONTEXT("When adding polytope to Mesh<>");
#ifdef FOUNDATION_DEBUG
          if(level + 1 > _num_levels)
            throw MeshInternalIndexOutOfBounds(level, _num_levels - 1);
#endif
          switch(level)
          {
            case pl_vertex:
              {
                if(_num_levels > pl_edge)
                {
                  _topologies.at(ipi_vertex_edge).push_back();

                }
                if(_num_levels > pl_face)
                {
                  _topologies.at(ipi_vertex_face).push_back();

                }
                if(_num_levels > pl_polyhedron)
                {
                  _topologies.at(ipi_vertex_polyhedron).push_back();

                }
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

            default:
              {
              }
          }
        }

        /**
         * \brief removes a polytope from the mesh
         *
         * \param[in] level
         * polytope level from which a polytope is to be removed
         *
         * \param[in] index
         * index of the polytope to be removed
         */
        void remove_polytope(const PolytopeLevels /*level*/, index_type_ /*i*/)
        {
          //TODO
        }

        /**
         * \brief adds an adjacency to the associated topology and automatically adds it to its transpose
         *
         * \param[in] from_level
         * polytope level (source of the relation)
         *
         * \param[in] to_level
         * polytope level (target of the relation)
         *
         * \param[in] polytope_index
         * index of the polytope (source)
         *
         * \param[in] value
         * index of the polytope (target)
         */
        void add_adjacency(
            const PolytopeLevels from_level,
            const PolytopeLevels to_level,
            const index_type_ polytope_index,
            const index_type_ value)
        {
          CONTEXT("When adding adjacency in Mesh<>");
#ifdef FOUNDATION_DEBUG
          if(from_level + 1 > _num_levels)
            throw MeshInternalIndexOutOfBounds(from_level, _num_levels - 1);

          if(to_level + 1 > _num_levels)
            throw MeshInternalIndexOutOfBounds(to_level, _num_levels - 1);

          if(from_level == to_level)
            throw MeshError("Explicit neighbourhood of polytopes on same level currently unsupported. Use implict storage instead.");//TODO
#endif

          InternalPolytopeIndices ipi(_get_internal_index(from_level, to_level));
          InternalPolytopeIndices ipit(_transpose_internal_index(ipi));

          _topologies.at(ipi).insert(polytope_index, value);
          _topologies.at(ipit).insert(value, polytope_index);

        }

        /**
         * \brief get all adjacent polytopes at level to_level for polytope i on level from_level
         *
         * \param[in] from_level
         * polytope level of polytope i
         *
         * \param[in] to_level
         * polytope level of the target adjacencies
         *
         * \param[in] i
         * index of the polytope for which the adjacencies are to be computed
         *
         */
        const typename TopologyType_::storage_type_ get_adjacent_polytopes(
                                                                     const PolytopeLevels from_level,
                                                                     const PolytopeLevels to_level,
                                                                     index_type_ i) const
        {
          CONTEXT("When calculating adjacent polytopes in Mesh<>");

          if(from_level == pl_none || to_level == pl_none)
            return typename TopologyType_::storage_type_();

#ifdef FOUNDATION_DEBUG
          if(from_level + 1 > _num_levels)
            throw MeshInternalIndexOutOfBounds(from_level, _num_levels - 1);
          if(to_level + 1 > _num_levels)
            throw MeshInternalIndexOutOfBounds(to_level, _num_levels - 1);
#endif
          InternalPrimaryAccess ipa(_get_primary_index(from_level, to_level));

          if(!_secondary_access_needed(from_level, to_level))
            return _topologies.at(ipa).at(i);
          else
          {
            typename TopologyType_::storage_type_ result;

            InternalSecondaryAccess isa(_get_secondary_index(from_level, to_level));
            for(index_type_ j(0) ; j < _topologies.at(ipa).at(i).size() ; ++j)
            {
              for(index_type_ k(0) ; k < _topologies.at(isa).at(_topologies.at(ipa).at(i).at(j)).size() ; ++k)
              {
                index_type_ to_insert(_topologies.at(isa).at(_topologies.at(ipa).at(i).at(j)).at(k));

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
                  bool match_super(true);
                  if(to_level < from_level)
                  {
                    //check if all vertices of to_insert are in super-polytope's vertex list
                    typename TopologyType_::storage_type_ vertices(get_adjacent_polytopes(to_level, pl_vertex, to_insert));
                    for(index_type_ l(0) ; l < vertices.size() ; ++l)
                    {
                      bool found(false);
                      for(index_type_ m(0) ; m < _topologies.at(ipa).at(i).size() ; ++m)
                      {
                        if(_topologies.at(ipa).at(i).at(m) == vertices.at(l))
                        {
                          found = true;
                          break;
                        }
                      }
                      if(!found)
                      {
                        match_super = false;
                      }
                    }
                  }
                  if(match_super)
                    result.push_back(to_insert);
                }
              }
            }
            return result;
          }
        }

        /**
         * \brief for two given polyopes at the same level, compute intersections
         *
         * \param[in] intersecting_elements_level
         * polytope level of the elements to check intersection for
         *
         * \param[in] intersection_level
         * polytope level of the intersection polytope
         *
         * \param[in] a
         * index of the first element
         *
         * \param[in] b
         * index of the second element
         */
        const typename TopologyType_::storage_type_ get_comm_intersection(
                                                                     const PolytopeLevels intersecting_elements_level,
                                                                     const PolytopeLevels intersection_level,
                                                                     index_type_ a,
                                                                     index_type_ b) const
        {
          CONTEXT("When calculating communication intersection in Mesh<>");

          typename TopologyType_::storage_type_ result;

          if(a == b)
            return result;

          typename TopologyType_::storage_type_ adj_a(get_adjacent_polytopes(intersecting_elements_level, intersection_level, a));
          typename TopologyType_::storage_type_ adj_b(get_adjacent_polytopes(intersecting_elements_level, intersection_level, b));

          for(index_type_ i(0) ; i < adj_a.size() ; ++i)
            for(index_type_ j(0) ; j < adj_b.size() ; ++j)
              if(adj_a.at(i) == adj_b.at(j))
                result.push_back(adj_a.at(i));

          return result;
        }

        ///needed public access functions
        typename TopologyType_::index_type_ get_num_levels()
        {
          return _num_levels;
        }

        ///Further needed access functions
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

        ///further needed getter fcts
        OuterStorageType_<TopologyType_, std::allocator<TopologyType_> >& get_topologies()
        {
          return _topologies;
        }

        const OuterStorageType_<TopologyType_, std::allocator<TopologyType_> >& get_topologies() const
        {
          return _topologies;
        }

        const index_type_ num_polytopes(PolytopeLevels level) const
        {
          switch(level)
          {
            case pl_vertex:
              return _topologies.at(ipi_vertex_edge).size();
            case pl_edge:
              return _topologies.at(ipi_edge_vertex).size();
            case pl_face:
              return _topologies.at(ipi_face_vertex).size();
            case pl_polyhedron:
              return _topologies.at(ipi_polyhedron_vertex).size();
            default:
              return 0;
          }
        }


      private:
        typename TopologyType_::index_type_ _id;
        typename TopologyType_::index_type_ _pp_rank;
        typename TopologyType_::index_type_ _mp_rank;
        unsigned _num_inter_topologies;
        unsigned _num_levels;

        OuterStorageType_<TopologyType_, std::allocator<TopologyType_> > _topologies;

        inline unsigned _level_difference(const unsigned from, const unsigned to)
        {
          return from == to ? 1u : ( from > to ? (unsigned)std::abs((double)(from - to)) - 1u : (unsigned)std::abs((double)(to - from)) - 1u);
        }

        InternalPrimaryAccess _get_primary_index(const PolytopeLevels from, const PolytopeLevels to) const
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

                  default:
                    {
                      return ipa_none;
                    }
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

            default:
              {
                return ipa_none;
              }
          }
        }

        InternalSecondaryAccess _get_secondary_index(PolytopeLevels from, PolytopeLevels to) const
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

        bool _secondary_access_needed(PolytopeLevels from, PolytopeLevels to) const
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

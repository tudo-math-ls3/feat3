#pragma once
#ifndef KERNEL_FOUNDATION_MESH_HH
#define KERNEL_FOUNDATION_MESH_HH 1

#include <kernel/foundation/base.hpp>
#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/mesh_error.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/foundation/buffer.hpp>
#include <kernel/foundation/communication.hpp>
#include <iostream>
#include <cmath>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief Attribute registration operation wrapper class template
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

    /**
     * \brief Required number of topologies for given number of spatial dimensions
     *
     * \author Markus Geveler
     */
    enum RequiredNumTopologies
    {
      rnt_1D = 2,
      rnt_2D = 4,
      rnt_3D = 6
    };

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
      ipa_polyhedron_any = 5
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
     * \tparam AttributeStorageType_
     * STL conformal container template type storing the attributes' data
     *
     * \tparam OuterAttributeStorageType_
     * STL conformal container template type storing the attributes
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
    class Mesh :
      public CommunicateableByAggregates<TopologyType_, com_send_receive>
    {
      public:
        friend class MeshAttributeRegistration;

        ///type exports
        typedef TopologyType_ topology_type_;
        typedef typename TopologyType_::index_type_ index_type_;
        typedef typename TopologyType_::storage_type_ storage_type_;
        typedef typename TopologyType_::buffer_type_ buffer_type_;

        typedef OuterAttributeStorageType_<
          std::shared_ptr<AttributeBase<AttributeStorageType_> >, std::allocator<std::shared_ptr<AttributeBase<AttributeStorageType_> > > > attr_base_type_;

        ///CTOR
        Mesh(const typename TopologyType_::index_type_ id, attr_base_type_* attrbase = nullptr, const typename TopologyType_::index_type_ pp_rank = 0, const typename TopologyType_::index_type_ mp_rank = 0) :
          _id(id),
          _pp_rank(pp_rank),
          _mp_rank(mp_rank),
          _num_inter_topologies(i_),
          _num_levels((unsigned)(i_/2u) + 1u),
          _topologies(i_),
          _history(true),
          _attrs(attrbase),
          _num_attributes(0),
          _attribute_polytopelevel_relations()
        {
        }

        ///Copy CTOR
        Mesh(const typename TopologyType_::index_type_ new_id, Mesh & other, attr_base_type_* attrbase = nullptr) :
          _id(new_id),
          _pp_rank(other._pp_rank),
          _mp_rank(other._mp_rank),
          _num_inter_topologies(other._num_inter_topologies),
          _num_levels(other._num_levels),
          _topologies(other._topologies),
          _history(),
          _attrs(attrbase),
          _num_attributes(0),
          _attribute_polytopelevel_relations(other._attribute_polytopelevel_relations)
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
          this-> _history = rhs._history;
          this-> _attrs = rhs._attrs;
          this-> _num_attributes = rhs._num_attributes;
          this-> _attribute_polytopelevel_relations = rhs._attribute_polytopelevel_relations;

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
                _history.add_functor(new CompoundFunctor<OuterStorageType_>(true)); //do not execute since this is done by low-level functor

                if(_num_levels > pl_edge)
                {
                  _topologies.at(ipi_vertex_edge).push_back();

                  ((CompoundFunctor<OuterStorageType_>*)(_history.get_functors().at(_history.size() - 1).get()))->add_functor(_topologies.at(ipi_vertex_edge).get_history().get_functors().at(_topologies.at(ipi_vertex_edge).get_history().size() - 1));
                }
                if(_num_levels > pl_face)
                {
                  _topologies.at(ipi_vertex_face).push_back();

                  ((CompoundFunctor<OuterStorageType_>*)(_history.get_functors().at(_history.size() - 1).get()))->add_functor(_topologies.at(ipi_vertex_face).get_history().get_functors().at(_topologies.at(ipi_vertex_face).get_history().size() - 1));
                }
                if(_num_levels > pl_polyhedron)
                {
                  _topologies.at(ipi_vertex_polyhedron).push_back();

                  ((CompoundFunctor<OuterStorageType_>*)(_history.get_functors().at(_history.size() - 1).get()))->add_functor(_topologies.at(ipi_vertex_polyhedron).get_history().get_functors().at(_topologies.at(ipi_vertex_polyhedron).get_history().size() - 1));
                }
              }
              break;

            case pl_edge:
              {
                  _topologies.at(ipi_edge_vertex).push_back();

                  _history.add_functor(_topologies.at(ipi_edge_vertex).get_history().get_functors().at(_topologies.at(ipi_edge_vertex).get_history().size() - 1));
              }
              break;

            case pl_face:
              {
                  _topologies.at(ipi_face_vertex).push_back();

                  _history.add_functor(_topologies.at(ipi_face_vertex).get_history().get_functors().at(_topologies.at(ipi_face_vertex).get_history().size() - 1));
              }
              break;

            case pl_polyhedron:
              {
                  _topologies.at(ipi_polyhedron_vertex).push_back();

                  _history.add_functor(_topologies.at(ipi_polyhedron_vertex).get_history().get_functors().at(_topologies.at(ipi_polyhedron_vertex).get_history().size() - 1));
              }
              break;
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

          _history.add_functor(new CompoundFunctor<OuterStorageType_>(true));
          ((CompoundFunctor<OuterStorageType_>*)(_history.get_functors().at(_history.size() - 1).get()))->add_functor(_topologies.at(ipi).get_history().get_functors().at(_topologies.at(ipi).get_history().size() - 1));

          ((CompoundFunctor<OuterStorageType_>*)(_history.get_functors().at(_history.size() - 1).get()))->add_functor(_topologies.at(ipit).get_history().get_functors().at(_topologies.at(ipit).get_history().size() - 1));
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

        /**
         * \brief undo last insertion or removal of polytopes or their adjacencies
         *
         */
        std::shared_ptr<FunctorBase> undo()
        {
          if(_history.size() == 0)
            throw MeshError("Already cleared!");

          _history.get_functors().at(_history.size() - 1).get()->undo();
          std::shared_ptr<FunctorBase> func(_history.get_functors().at(_history.size() - 1));
          _history.get_functors().pop_back();
          return func;
        }

        /**
         * \brief reset all data (undo the complete history)
         *
         */
        void reset()
        {
          if(_history.size() == 0)
            throw MeshError("Already cleared!");

          _history.undo();
        }

        /**
         * \brief clear all data (undo the complete history) but return it
         *
         */
        CompoundFunctor<OuterStorageType_> clear()
        {
          if(_history.size() == 0)
            throw MeshError("Already cleared!");

          _history.get_functors().clear();

          for(index_type_ i(0) ; i < i_ ; ++i)
            _topologies.at(i).get_topology().clear();

          return _history;
        }


        /**
         * \brief redo last undone action
         *
         */
        void redo(const std::shared_ptr<FunctorBase> func)
        {
          func.get()->execute();
          _history.add_functor(func);
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

        CompoundFunctor<OuterStorageType_>& get_history()
        {
          return _history;
        }

        const CompoundFunctor<OuterStorageType_>& get_history() const
        {
          return _history;
        }

      private:
        typename TopologyType_::index_type_ _id;
        typename TopologyType_::index_type_ _pp_rank;
        typename TopologyType_::index_type_ _mp_rank;
        unsigned _num_inter_topologies;
        unsigned _num_levels;

        OuterStorageType_<TopologyType_, std::allocator<TopologyType_> > _topologies;
        CompoundFunctor<OuterStorageType_> _history;

        attr_base_type_* _attrs;

        unsigned _num_attributes;

        typename TopologyType_::storage_type_ _attribute_polytopelevel_relations;

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

          throw InternalError("Invalid polytope level combination");
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

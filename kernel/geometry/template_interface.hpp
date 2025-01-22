// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include "kernel/geometry/intern/congruency_sampler.hpp"
#include "kernel/geometry/intern/face_index_mapping.hpp"
#include <kernel/geometry/intern/adaptive_mesh_storage.hpp>
#include <kernel/shape.hpp>

namespace FEAT::Geometry
{
  /**
   * \brief Mesh Query. Allows access to local mesh entities for template queries.
   *
   * \author Markus Muegge
   */
  template<typename MeshStorage_, typename QueryShape_>
  class MeshQuery
  {
    // This class is in a strange place where it both requires access to
    // low-level details of the adaptive mesh and is part of the public API. To
    // avoid leaking implementation details directly most types are typedefed
    // below. The constructor is marked private to avoid making the mesh storage
    // part of the public API. AdaptiveMesh is marked as a friend to allow it to
    // construct instances of MeshQuery.
    template<
      typename TemplateSet_,
      typename Shape_,
      int num_coords_,
      typename Coord_>
    friend class AdaptiveMesh;

    /// Dimension of element this query is about
    static constexpr int dimension = QueryShape_::dimension;

    /// Mesh storage type
    using MeshStorage = MeshStorage_;

    /// Topology of this query
    using Topology = typename MeshStorage::template ElementTopologyByDim<dimension>;

    /// Children of this query
    using Children = typename MeshStorage::template ElementChildrenByDim<dimension>;

    /// Template set being used
    using TemplateSet = typename MeshStorage::TemplateSet;

  public:
    /// Shape of entity this query is about
    using QueryShape = QueryShape_;

    /// Vertex type of mesh this query is about
    using VertexType = typename MeshStorage_::VertexType;

    /// Type of array of vertex coordinates
    using VertexCoordinateArray = std::array<VertexType, Shape::FaceTraits<QueryShape_, 0>::count>;

    /// Accessor for topology types
    template<int dim_>
    using TopologyByDim = typename MeshStorage::template ElementTopologyByDim<dim_>;

    /// Accessor for reference types
    template<int dim_>
    using RefByDim = typename MeshStorage::template ElementRefByDim<dim_>;

    /// Accessor for oriented reference types
    template<int dim_>
    using OrientedRefByDim = typename MeshStorage::template OrientedElementRefByDim<dim_>;

    /// Oriented vertex reference type
    using OrientedVertex = OrientedRefByDim<0>;
    /// Oriented edge reference type
    using OrientedEdge = OrientedRefByDim<1>;
    /// Oriented face reference type
    using OrientedFace = OrientedRefByDim<2>;

    //// Vertex reference type
    using VertexRef = RefByDim<0>;
    /// Edge reference type
    using EdgeRef = RefByDim<1>;
    /// Face reference type
    using FaceRef = RefByDim<2>;
    /// Cell reference type
    using CellRef = RefByDim<3>;

  private:
    /// Storage of mesh
    MeshStorage& _storage;
    /// Topology of entity being queried
    Topology& _topology;
    /// Children of entity being queried
    Children& _children;

    /**
     * \brief Constructor
     *
     * \param[in] storage Mesh storage
     * \param[in] topology Topology of the entity being queried
     * \param[in] children Children of the entity being queried
     */
    MeshQuery(MeshStorage& storage, Topology& topology, Children& children) :
      _storage(storage),
      _topology(topology),
      _children(children)
    {
    }

  public:
    MeshQuery() = delete;

    /**
     * \brief Accessor for vertex coordinates
     *
     * \returns An array of vertex coordinates
     *
     * \warning A new array is created on each call. Cache the result if you
     * require repeated access to the array.
     */
    VertexCoordinateArray vertex_coordinates() const
    {
      VertexCoordinateArray result;

      auto& vertices = _topology.template by_dim<0>();
      for(int i{0}; i < _topology.template size<0>(); i++)
      {
        result[i] = _storage[vertices[i]].vertex;
      }
      return result;
    }

    /**
     * \brief Retrieve refineent type of an entity
     *
     * \tparam dim_ The dimension of the entity
     * \param[in] ref Reference to the entity
     *
     * \returns The refinement type of the entity
     */
    template<int dim_>
    typename TemplateSet::template RefinementTypeByDim<dim_> type(RefByDim<dim_> ref) const
    {
      return _storage[ref].type;
    }

    /**
     * \brief Determines to orientat
     */
    template<int topo_dim_, int entity_dim_>
    int orientation(const TopologyByDim<topo_dim_>& topology, int elem) const
    {
      static_assert(entity_dim_ < topo_dim_);

      using TopoShape = typename Shape::FaceTraits<QueryShape, topo_dim_>::ShapeType;
      using EntityShape = typename Shape::FaceTraits<QueryShape, entity_dim_>::ShapeType;

      using LocalMapping = Intern::FaceIndexMapping<TopoShape, entity_dim_, 0>;
      using CongruencySampler = Intern::CongruencySampler<EntityShape>;

      static constexpr int num_vertices = Shape::FaceTraits<EntityShape, 0>::count;

      std::array<RefByDim<0>, num_vertices> local;
      std::array<RefByDim<0>, num_vertices> foreign;

      const auto& entity = _storage[topology.template key_by_dim<entity_dim_>(elem)];

      for(int i(0); i < num_vertices; i++)
      {
        local[i] = topology.template key_by_dim<0>(LocalMapping::map(elem, i)).key;
        foreign[i] = entity.topology.template key_by_dim<0>(i).key;
      }

      int result = CongruencySampler::compare(local.data(), foreign.data());
      if(result == -1)
      {
        for(int i(0); i < num_vertices; ++i)
        {
          for(int j(0); j < num_vertices; ++j)
          {
            std::cout << "Comparing:\n";
            std::cout << "  Local: " << local[i] << "\n";
            std::cout << "  Foreign: " << foreign[j] << "\n";
          }
        }

        XABORTM("Invalid orientation found.");
      }
      return result;
    }

    /**
     * \brief Determine topology size of current entity
     *
     * \tparam dim_ The dimension to query
     *
     * \returns Number of entities of dimension \c dim in the topology
     */
    template<int dim_>
    Index topology_size() const
    {
      return Shape::FaceTraits<QueryShape, dim_>::count;
    }

    /**
     * \brief Topology accessor. Retrieves a reference from the topology.
     *
     * The returned reference points to an entity one layer above an entity
     * currently being constructed by a template query.
     *
     * \tparam dim_ Dimension of entity asked for
     * \param[in] idx Index (in topology) of entity asked for
     *
     * \return A reference to the entity
     */
    template<int dim_>
    OrientedRefByDim<dim_> topology(Index idx) const
    {
      static_assert(dim_ >= 0, "Invalid shape dimension in MeshQuery::get_from_topology()!");
      static_assert(dim_ <= 3, "Invalid shape dimension in MeshQuery::get_from_topology()!");

      return _topology.template key_by_dim<dim_>(idx);
    }

    /**
     * \brief Boundary accessor. Retrieves a reference from a child of an entity of the topology.
     *
     * The returned reference points to an entity on the same mesh layer as a
     * entity currently being constructed by a template query.
     *
     * \tparam parent_dim_ Dimension of topology entity
     * \tparam child_dim_ Dimension of child of topology entity
     * \param[in] parent_idx Index (in topology) of topology entity
     * \param[in] child_idx Index (in children) of topology child
     *
     * \return A reference to the entity
     */
    template<int parent_dim_, int child_dim_>
    RefByDim<child_dim_> boundary(Index parent_idx, Index child_idx) const
    {
      static_assert(child_dim_ <= parent_dim_);
      static_assert(child_dim_ >= 0);

      auto& parent = _storage[topology<parent_dim_>(parent_idx).key];
      return parent.children.template by_dim<child_dim_>()[child_idx];
    }

    /**
     * \brief Sibling accessor. Retrieves a reference from the children of the current entity.
     *
     * The returned reference points to an entity on the same mesh layer as a
     * entity currently being constructed by a template query.
     *
     * \tparam dim_ Dimension of sibling entity
     * \param[in] idx Index (in children) of entity asked for
     *
     * \return A reference to the entity
     */
    template<int dim_>
    RefByDim<dim_> sibling(Index sibling_idx) const
    {
      static_assert(dim_ >= 0, "Invalid shape dimension in MeshQuery::get_from_topology()!");
      static_assert(dim_ <= 3, "Invalid shape dimension in MeshQuery::get_from_topology()!");

      return _children.template by_dim<dim_>()[sibling_idx];
    }
  };

#ifdef DOXYGEN
  /**
   * \brief Template query for vertices.
   *
   * This interface is expected by the AdaptiveMesh to construct all vertices
   * required to refine a mesh entity.
   *
   * \author Markus Muegge
   */
  class VertexQuery
  {
    /**
     * \brief Returns the number of vertices to construct with this query
     *
     * \param[in] query Mesh query
     *
     * \return Number of vertices to construct
     */
    Index num_entities(const MeshQuery& query);

    /**
     * \brief Constructs a vertex
     *
     * Must return a vertex for all indices idx with 0 <= idx < num_entities(query)
     *
     * \param[in] query Mesh query
     * \param[in] vertex_idx Index of vertex to construct
     *
     * \returns A vertex
     */
    typename MeshQuery_::VertexType_ vertex(const MeshQuery_& query, Index vertex_idx)
  };

  /*
   * \brief Template query for vertices.
   *
   * This interface is expected by the AdaptiveMesh for constructing edges,
   * faces, and cells.
   *
   * \author Markus Muegge
   */
  template<int dim_>
  class EntityQuery
  {
    using Topology = typename MeshQuery_::template TopologyByDim<dim_>;

    /**
     * \brief Returns the number of entities to construct with this query
     *
     * \param[in] query Mesh query
     *
     * \return Number of vertices to construct
     */
    Index num_entities(const MeshQuery& query);

    /**
     * \brief Constructs a topology
     *
     * Must return a valid topology for all indices idx with 0 <= idx <= num_entities(query)
     *
     * \param[in] query Mesh query
     * \param[in] idx Index of entity to construct
     */
    Topology topology(const MeshQuery& query, Index idx)

      /**
       * \brief Determine subdivision levels for child entity
       *
       * \param[in] query Mesh query
       * \param[in] idx Index of child entity to determine levels for
       *
       * \returns Subdivision levels for the entity
       */
      SubdivisionLevelTuple levels(const MeshQuery_& /*query*/, Index edge_idx)
    {
      return {
        _level_spread[_template.template vertex_of<1>(edge_idx, 0)],
        _level_spread[_template.template vertex_of<1>(edge_idx, 1)],
      };
    }
  };

  /**
   * \brief Template set interface
   *
   * This is the interface for a template set expected by the AdaptiveMesh.
   *
   * \author Markus Muegge
   */
  class TemplateSet
  {
  public:
    /// Refinement type of an element requiring no further refinement
    static constexpr std::uint8_t zero_refinement_type = 0;
    /// Refinement type of an element requiring full refinement
    static constexpr std::uint8_t full_refinement_type = 0b11111111;

    /**
     * \brief Shape compatability test
     *
     * Used by the AdaptiveMesh to ensure this template set is compatible with the meshes shape type.
     *
     * \tparam Shape_ Mesh shape this template set is to be applied to
     *
     * \returns True, if the template set is compatible with the mesh shape, false otherwise.
     */
    template<typename Shape_>
    static constexpr bool is_shape_compatible();

    /**
     * \brief Returns maximum number of children a template produces
     *
     * \tparam template_dim_ Dimension of the template shape
     * \tparam child_dim_ Dimension of the children
     *
     * Used by the AdaptiveMesh to determine the size of child data structures at compile time.
     */
    template<int template_dim_, int child_dim_>
    static constexpr int max_children();

    /**
     * \brief Computes refinement type of a SubdivisionLevelTuple
     *
     * The refinement type is meant as a summary of the vertex markings. It is used to compare different refinements.
     *
     * \tparam n Length of tuple
     *
     * \param[in] levels Subdivision level tuple
     */
    template<int n>
    static std::uint8_t refinement_type(const Geometry::SubdivisionLevelTuple<n>& levels);

    /**
     * \brief Adjusts SubdivisionLevels to match TemplateSet requirements
     *
     * A templat set might have special requirements for the vertex markings.
     * This method will be called by the AdaptiveMesh before any refinement
     * takes place. It allows the template set to modify the vertex markings to
     * ensure all requiremens are met.
     *
     * This method should only add subdivision levels, to ensure the user get
     * at least the refinement he asked for.
     *
     * \param[inout] sdls SubdivisionLevels to adjust
     */
    template<typename Mesh_>
    static void convert(const Mesh_& mesh, Geometry::SubdivisionLevels& sdls);

    /**
     * \brief Returns number of children a specific template produces
     *
     * This is intended as a utility function for iterating over children without creating full TemplateQuerys.
     *
     * \tparam template_dim_ Template dimension
     * \tparam child_dim_ Child dimension
     *
     * \param[in] type Type of template
     *
     * \returns Number of children of dimension \c child_dim_ of template with dimension \c template_dim_ and type \c
     * type
     */
    template<int template_dim_, int child_dim_>
    static Index num_children(std::uint8_t type);

    /**
     * \brief Create Vertex-Query for instantiating a template
     *
     * This method is the start of the query chain. The AdaptiveMesh will use
     * it to create a vertex query for any template it wants to apply.
     * The AdaptiveMesh expects the object returned by this method to conform
     * to the VertexQuery interface.
     *
     * \tparam MeshQuery_ Mesh query type
     *
     * \param[in] query Mesh query
     * \param[in] levels Subdivision levels
     */
    template<typename MeshQuery_>
    static TemplateQuery<MeshQuery_, 0> get_query(
      const MeshQuery_& query,
      Geometry::SubdivisionLevelTuple<Shape::FaceTraits<typename MeshQuery_::QueryShape, 0>::count> levels);

    /**
     * \brief Advances a template query
     *
     * This method is meant to take, for example, a vertex query and advance it
     * into an edge query. This gives the template set a chance to do any setup
     * it needs to do before the next dimension of entities is to be created.
     * The AdaptiveMesh make no assumptions on the types returned by this
     * method. The template set can, if it so wishes, return a completely
     * independent type for each kind of query.
     *
     * The AdaptiveMesh expects the objects returned by this method to conform
     * to the EntityQuery interface.
     *
     * The old query can be moved from, for example to keep cached data.
     *
     * \tparam MeshQuery_ Mesh query type
     * \tparam n_ Current query stage
     *
     * \param[in] query Mesh query
     * \param[in] template_query Query at currentstage
     *
     * \returns Template query for next higher dimension of entities.
     */
    static NextQuery advance_query(const MeshQuery_& query, PreviousQuery&& template_query)

#endif
  } // namespace FEAT::Geometry

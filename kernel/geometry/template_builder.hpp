// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/geometry/intern/congruency_mapping.hpp>
#include <kernel/geometry/intern/congruency_sampler.hpp>
#include <kernel/geometry/intern/congruency_trafo.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/geometry/raw_refinement_templates.hpp>
#include <kernel/geometry/refinement_types.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <iostream>

namespace FEAT::Geometry
{
  namespace Intern
  {
    /**
     * \brief Constexpr function for number of possible orientations of a shape
     */
    template<typename Shape_>
    constexpr int orientations()
    {
      if constexpr(std::is_same_v<Shape_, Shape::Hypercube<1>>)
      {
        return 2;
      }

      if constexpr(std::is_same_v<Shape_, Shape::Hypercube<2>>)
      {
        return 8;
      }

      return 0;
    }

    /**
     * \brief Comparison predicate for vertices.
     *
     * \param[in] a First vertex
     * \param[in] b Second vertex
     *
     * \returns Returns true if not coordinate of the vertices differs by more than an epsilon.
     */
    template<int n_>
    bool is_same_vertex(const Tiny::Vector<Real, n_>& a, const Tiny::Vector<Real, n_>& b)
    {
      bool is_same = true;

      const Real tol = Math::pow(Math::eps<Real>(), Real(0.7));
      for(int i(0); i < n_; ++i)
      {
        if(Math::abs(a[i] - b[i]) > tol)
        {
          is_same = false;
        }
      }

      return is_same;
    }

    /**
     * \brief Comparison predicate for raw entities
     *
     * Two raw entities are the same of they contains the same vertices in any order.
     * This is enough for our purposes, as vertex order will be corrected by the
     * OrientationMappings later.
     *
     * \param[in] a First entity
     * \param[in] b Second Entity
     *
     * \returns True, if the entities coordinates are permutations of each other,
     * according to is_same_vertex
     */
    template<typename Shape_, int num_coords_ = Shape_::dimension>
    bool is_same_raw_entity(const RawEntity<Shape_, num_coords_>& a, const RawEntity<Shape_, num_coords_>& b)
    {
      return std::is_permutation(a.coords.begin(), a.coords.end(), b.coords.begin(), Intern::is_same_vertex<num_coords_>);
    }

    /**
     * \brief Returns a vector of all vertices used in the given entities.
     *
     * \param[in] entities List of raw entities
     *
     * \note The vertices are not deduplicated. If a vertex is used by multiple
     * entities, it will be included multiple times.
     */
    template<typename Shape_, int num_coords_ = Shape_::dimension>
    std::vector<Tiny::Vector<Real, num_coords_>> all_vertices(const std::vector<RawEntity<Shape_, num_coords_>>& entities)
    {
      std::vector<Tiny::Vector<Real, num_coords_>> vertices;

      for(const auto& entity : entities)
      {
        for(const auto& vertex : entity.coords)
        {
          vertices.push_back(vertex);
        }
      }

      return vertices;
    }

    /**
     * \brief Internal vertex predicate
     *
     * \param[in] vertex Vertex to check
     *
     * An internal vertex is any vertex that does not lie on the boundary of its
     * shape. This implementation is for Hypercubes, so it ensures all
     * coordinates lie in (0, 1).
     *
     * \returns True, if the given vertex is internal, false otherwise.
     */
    template<int n_>
    bool is_internal_vertex(Tiny::Vector<Real, n_> vertex)
    {
      const Real tol = Math::pow(Math::eps<Real>(), Real(0.7));

      // A vertex is internal if all its coordinates are in (0, 1)
      bool is_internal = true;
      for(int i(0); i < n_; ++i)
      {
        if(vertex[i] < tol || (1 - tol) < vertex[i])
        {
          is_internal = false;
        }
      }
      return is_internal;
    }

    /**
     * \brief Internal entity predicate
     *
     * \param[in] entity Entity to check
     *
     * Calculates the entities centroid and checks whether it is internal. This
     * covers edge cases such as an edge running from one boundary to another,
     * crossing the interior. Such an edge is internal, but all its vertices are
     * not.
     *
     * \returns True, if the given entity is internal, false otherwise.
     */
    template<typename Shape_, int num_coords_ = Shape_::dimension>
    bool is_internal_entity(const RawEntity<Shape_, num_coords_>& entity)
    {
      Tiny::Vector<Real, num_coords_> centroid(0);
      for(const auto& vertex : entity.coords)
      {
        centroid += vertex;
      }
      return is_internal_vertex((1.0 / entity.coords.size()) * centroid);
    }

    /// Deduplicates the given vector, using Comp for comparisons.
    template<typename T_, typename Comp>
    std::vector<T_> deduplicate(Comp comparison, const std::vector<T_>& list)
    {
      std::vector<T_> deduplicated;
      for(const T_& current : list)
      {
        auto pred = [&](T_& t) { return comparison(current, t); };
        bool is_new = !std::any_of(deduplicated.begin(), deduplicated.end(), pred);
        if(is_new)
        {
          deduplicated.push_back(current);
        }
      }
      return deduplicated;
    }

    /// Filters the given vector
    template<typename T_, typename Predicate_>
    std::vector<T_> filter(Predicate_ predicate, const std::vector<T_>& list)
    {
      std::vector<T_> filtered;
      for(const T_& current : list)
      {
        if(predicate(current))
        {
          filtered.push_back(current);
        }
      }
      return filtered;
    }

    /// Returns a deduplicated list of internal vertices of the template
    template<typename Shape_>
    std::vector<Tiny::Vector<Real, Shape_::dimension>> internal_vertices(const RawTemplate<Shape_>& tmplt)
    {
      std::vector<Tiny::Vector<Real, Shape_::dimension>> vertices;

      for(const auto& entity : tmplt.entities)
      {
        for(const auto& vertex : entity.coords)
        {
          vertices.push_back(vertex);
        }
      }

      vertices = filter([](auto vertex) { return is_internal_vertex(vertex); }, vertices);
      vertices = deduplicate([](auto a, auto b) { return Intern::is_same_vertex(a, b); }, vertices);

      return vertices;
    }

    /// Returns a deduplicated list of internal entities of the template
    template<typename Shape_, int face_dim_, int num_coords_ = Shape_::dimension>
    std::vector<RawEntity<typename Shape::FaceTraits<Shape_, face_dim_>::ShapeType, num_coords_>>
    internal_entities(const std::vector<RawEntity<Shape_, num_coords_>>& entities)
    {
      if constexpr(Shape_::dimension == face_dim_)
      {
        return entities;
      }
      else
      {
        static constexpr const int faces_per_entity = Shape::FaceTraits<Shape_, face_dim_>::count;

        std::vector<RawEntity<typename Shape::FaceTraits<Shape_, face_dim_>::ShapeType, num_coords_>> result;

        for(const auto& entity : entities)
        {
          for(int i(0); i < faces_per_entity; ++i)
          {
            result.push_back(entity.template face<face_dim_>(i));
          }
        }

        result = filter([](auto& entity) { return is_internal_entity(entity); }, result);
        result = deduplicate([](auto& a, auto& b) { return Intern::is_same_raw_entity(a, b); }, result);

        return result;
      }
    }

    /// Transforms a given vertex to match the given orientation
    template<typename Shape_>
    static Tiny::Vector<Real, Shape_::dimension>
    orient_vertex(const Tiny::Vector<Real, Shape_::dimension>& vertex, int orientation)
    {
      static constexpr int num_coords = Shape_::dimension;
      using Trafo = Intern::CongruencyTrafo<Shape_>;

      Tiny::Vector<Real, num_coords> offset(0.5);

      Tiny::Matrix<Real, num_coords, num_coords> transform(0.0);
      Tiny::Vector<Real, num_coords> unused(0.0);
      Trafo::compute(transform, unused, orientation);

      Tiny::Vector<Real, num_coords> result;
      result.set_mat_vec_mult(transform, (vertex - offset));

      return result + offset;
    }

    /// Transforms a given RawEntity to match the given orientation
    template<typename Shape_, int num_coords_ = Shape_::dimension>
    static RawEntity<Shape_, num_coords_> orient_raw_entity(const RawEntity<Shape_, num_coords_>& entity, int orientation)
    {
      RawEntity<Shape_, num_coords_> result;
      for(int i(0); i < entity.coords.size(); i++)
      {
        result.coords[i] = Intern::orient_vertex<Shape_>(entity.coords[i], orientation);
      }
    }


    /**
     * \brief Embeds a vertex of a lower dimensional shape in a higher dimensional shape.
     *
     * Embedds for example an edge vertex (out of [0, 1]) into a face ([0, 1]^2).
     *
     * \tparam Shape_ Shape to embed into
     * \tparam souce_dim_ Dimension to embed from
     *
     * See specializations for details.
     */
    template<typename Shape_, int source_dim_>
    struct VertexEmbedder
    {
    };

    /// VertexEmbedder specialization for embedding edges in faces
    template<>
    struct VertexEmbedder<Shape::Quadrilateral, 1>
    {
      using SourceVertex = Tiny::Vector<Real, 1>;
      using TargetVertex = Tiny::Vector<Real, 2>;

      /**
       * \brief Embed an edge vertex into a Quadrilateral
       *
       * \parma[in] elem Index of the edge in the face, determines where the
       * vertex is embedded.
       * \param[in] vertex Vertex to embedded
       */
      static TargetVertex embed(int elem, const SourceVertex& vertex)
      {
        switch(elem)
        {
        case 0:
          return Tiny::Vector<Real, 2>{vertex[0], 0.0};
        case 1:
          return Tiny::Vector<Real, 2>{vertex[0], 1.0};
        case 2:
          return Tiny::Vector<Real, 2>{0.0, vertex[0]};
        case 3:
          return Tiny::Vector<Real, 2>{1.0, vertex[0]};
        default:
          XABORTM("Invalid face-edge index.");
        }
      }
    };

    /// VertexEmbedder specialization for embedding edges in cells
    template<>
    struct VertexEmbedder<Shape::Hexahedron, 1>
    {
      using SourceVertex = Tiny::Vector<Real, 1>;
      using TargetVertex = Tiny::Vector<Real, 3>;

      /**
       * \brief Embed an edge vertex into a Hexahedron
       *
       * \parma[in] elem Index of the edge in the face, determines where the
       * vertex is embedded.
       * \param[in] vertex Vertex to embedded
       */
      static TargetVertex embed(int elem, const SourceVertex& vertex)
      {
        switch(elem)
        {
        case 0:
          return Tiny::Vector<Real, 3>{vertex[0], 0.0, 0.0};
        case 1:
          return Tiny::Vector<Real, 3>{vertex[0], 1.0, 0.0};
        case 2:
          return Tiny::Vector<Real, 3>{vertex[0], 0.0, 1.0};
        case 3:
          return Tiny::Vector<Real, 3>{vertex[0], 1.0, 1.0};
        case 4:
          return Tiny::Vector<Real, 3>{0.0, vertex[0], 0.0};
        case 5:
          return Tiny::Vector<Real, 3>{1.0, vertex[0], 0.0};
        case 6:
          return Tiny::Vector<Real, 3>{0.0, vertex[0], 1.0};
        case 7:
          return Tiny::Vector<Real, 3>{1.0, vertex[0], 1.0};
        case 8:
          return Tiny::Vector<Real, 3>{0.0, 0.0, vertex[0]};
        case 9:
          return Tiny::Vector<Real, 3>{1.0, 0.0, vertex[0]};
        case 10:
          return Tiny::Vector<Real, 3>{0.0, 1.0, vertex[0]};
        case 11:
          return Tiny::Vector<Real, 3>{1.0, 1.0, vertex[0]};
        default:
          XABORTM("Invalid cell-edge index.");
        }
      }
    };

    /// VertexEmbedder specialization for embedding faces in cells
    template<>
    struct VertexEmbedder<Shape::Hexahedron, 2>
    {
      using SourceVertex = Tiny::Vector<Real, 2>;
      using TargetVertex = Tiny::Vector<Real, 3>;

      /**
       * \brief Embed a face vertex into a Hexahedron
       *
       * \parma[in] elem Index of the edge in the face, determines where the
       * vertex is embedded.
       * \param[in] vertex Vertex to embedded
       */
      static TargetVertex embed(int elem, const SourceVertex& vertex)
      {
        switch(elem)
        {
        case 0:
          return Tiny::Vector<Real, 3>{vertex[0], vertex[1], 0.0};
        case 1:
          return Tiny::Vector<Real, 3>{vertex[0], vertex[1], 1.0};
        case 2:
          return Tiny::Vector<Real, 3>{vertex[0], 0.0, vertex[1]};
        case 3:
          return Tiny::Vector<Real, 3>{vertex[0], 1.0, vertex[1]};
        case 4:
          return Tiny::Vector<Real, 3>{0.0, vertex[0], vertex[1]};
        case 5:
          return Tiny::Vector<Real, 3>{1.0, vertex[0], vertex[1]};
        default:
          XABORTM("Invalid cell-edge index.");
        }
      }
    };

    /// Utility function to construct all vertices of a Hypercube topology
    template<typename Shape_>
    std::array<Tiny::Vector<Real, Shape_::dimension>, Shape::FaceTraits<Shape_, 0>::count> topology_vertices()
    {
      static constexpr int num_coords = Shape_::dimension;
      static constexpr int num_vertices = Shape::FaceTraits<Shape_, 0>::count;

      std::array<Tiny::Vector<Real, num_coords>, num_vertices> vertices;

      for(int i(0); i < num_vertices; i++)
      {
        // Build topology vertex from reference cell
        Tiny::Vector<Real, num_coords> vertex;
        for(int j(0); j < num_coords; ++j)
        {
          vertex[j] = 0.5 * (Shape::ReferenceCell<Shape_>::template vertex<Real>(i, j) + 1.0);
        }
        vertices[i] = vertex;
      }

      return vertices;
    }

  } // namespace Intern

  ////////////////////////////////////
  // Refinement template data types
  ////////////////////////////////////

  /**
   * \brief Enumeration of possible entity sources
   *
   * During adaptive refinement, entity topologies have to be found. This enum
   * indicates where an entity can be found.
   */
  enum class EntitySource : std::uint8_t
  {
    ParentTopology, ///< Entity is part of the parents topology
    Sibling,        ///< Entity is a sibling in the refinement tree
    BoundaryEdge,   ///< Entity is part of a boundary edge, i.e. a child of one of the parents edges
    BoundaryFace,   ///< Entity is part of a boundary face, i.e. a child of one of the parents faces
  };

  /// Output operator for EntitySource
  inline std::ostream& operator<<(std::ostream& stream, EntitySource source)
  {
    static std::unordered_map<EntitySource, String> mapping{
      {EntitySource::ParentTopology, "EntitySource::ParentTopology"},
      {EntitySource::Sibling, "EntitySource::Sibling"},
      {EntitySource::BoundaryEdge, "EntitySource::BoundaryEdge"},
      {EntitySource::BoundaryFace, "EntitySource::BoundaryFace"},
    };
    stream << mapping[source];
    return stream;
  }

  /**
   * \brief Reference to another local mesh entity.
   *
   * This is used to define how to retrieve the topologies of entities created during adpative refinement.
   */
  struct EntityReference
  {
    /// Where to retrieve the entity from
    EntitySource source;
    /// Index of the entity in the target location. Might need to be mapped to
    /// account for template types and orientations.
    Index index;
    /// Orientation of the retrieved entity relative to the retrieving entity. Only set for sibling entities.
    int orientation;
    /// Index of boundary entity to retrieve entity from. Only set for boundary entities.
    int entity;

    EntityReference() = default;
    EntityReference(EntitySource s, Index i, int o, int e) : source(s), index(i), orientation(o), entity(e)
    {
    }

    // Equals operator for EntityReferences
    bool operator==(const EntityReference& rhs) const
    {
      return (source == rhs.source) && (index == rhs.index) && (orientation == rhs.orientation) &&
            (entity == rhs.entity);
    }

    // Not-equals operator for EntityReferences
    bool operator!=(const EntityReference& rhs) const
    {
      return !(*this == rhs);
    }
  };

  /// Output operator for EntityReferences
  inline std::ostream& operator<<(std::ostream& stream, const EntityReference& reference)
  {
    stream << "EntityReference {"
           << "source: " << reference.source << ", index: " << stringify(reference.index)
           << ", orientation: " << stringify(reference.orientation) << ", entity: " << stringify(reference.entity)
           << "}";
    return stream;
  }

  /**
   * \brief Construction rule for a single entities topology.
   *
   * \tparam Shape_ Shape type of the entity this topology describes.
   * \tparam reference_dim_ Dimension of entities that references refer to
   *
   * This is a dimension-recursive data structure that stores references to all
   * elements of the topology of an entity with shape \c Shape_.
   */
  template<typename Shape_, int reference_dim_ = Shape_::dimension - 1>
  struct TopologyTemplate : TopologyTemplate<Shape_, reference_dim_ - 1>
  {
    /// Number of references of dimension \c reference_dim_
    static constexpr int num_entities = Shape::FaceTraits<Shape_, reference_dim_>::count;

    /// Array of references
    std::array<EntityReference, num_entities> references;

    /**
     * \brief Reference-array accessor
     *
     * \tparam dim_ Dimension of entities to retrieve references for
     */
    template<int dim_>
    std::array<EntityReference, Shape::FaceTraits<Shape_, dim_>::count>& get_references()
    {
      return TopologyTemplate<Shape_, dim_>::references;
    }

    /**
     * \brief Const reference-array accessor
     *
     * \tparam dim_ Dimension of entities to retrieve references for
     */
    template<int dim_>
    const std::array<EntityReference, Shape::FaceTraits<Shape_, dim_>::count>& get_references() const
    {
      return TopologyTemplate<Shape_, dim_>::references;
    }

    template<int dim_>
    EntityReference get_reference(int idx) const
    {
      return TopologyTemplate<Shape_, dim_>::references[idx];
    }
  };

  /**
   * \brief Construction rule for a single entities topology.
   *
   * \tparam Shape_ Shape type of the entity this topology describes.
   * \tparam reference_dim_ Dimension of entities that references refer to
   *
   * This is a dimension-recursive data structure that stores references to all
   * elements of the topology of an entity with shape \c Shape_.
   *
   * This is the base case of the recursion.
   */
  template<typename Shape_>
  struct TopologyTemplate<Shape_, 0>
  {
    /// Number of vertex references
    static constexpr int num_entities = Shape::FaceTraits<Shape_, 0>::count;

    /// Array of vertex references
    std::array<EntityReference, num_entities> references;

    /**
     * \brief Reference-array accessor
     *
     * \tparam dim_ Dimension of entities to retrieve references for
     */
    template<int dim_>
    std::array<EntityReference, num_entities>& get_references()
    {
      static_assert(dim_ == 0, "invalid reference dimension");
      return references;
    }

    /**
     * \brief Const reference-array accessor
     *
     * \tparam dim_ Dimension of entities to retrieve references for
     */
    template<int dim_>
    const std::array<EntityReference, num_entities>& get_references() const
    {
      static_assert(dim_ == 0, "invalid reference dimension");
      return references;
    }

    template<int dim_>
    EntityReference get_reference(Index idx) const
    {
      return TopologyTemplate<Shape_, dim_>::references[idx];
    }
  };

  /**
   * \brief Template for refining a mesh entity.
   *
   * This is a collection of TopologyTemplates, which describe how to construct
   * all children of all dimension of the template.
   *
   * \tparam Shape_ Shape of this template
   * \tparam CoShape_ Current shape of this data structure definition.
   *
   * This is a dimension-recursive data structure via recursion over the \c CoShape_.
   */
  template<typename Shape_, typename CoShape_ = Shape_>
  struct RefinementTemplate
    : RefinementTemplate<Shape_, typename Shape::FaceTraits<CoShape_, CoShape_::dimension - 1>::ShapeType>
  {
    /// Accessor for topology types
    template<int dim_>
    using TopologyTemplateByDim = TopologyTemplate<typename Shape::FaceTraits<Shape_, dim_>::ShapeType>;

    /**
     * \brief Entities of shape CoShape_ added by this template.
     *
     * Presented as topology templates, which allow collecting all entities
     * making up the templates children.
     */
    std::vector<TopologyTemplate<CoShape_>> topologies;

    /**
     * \brief TopologyTemplate accessor
     *
     * \tparam Dimension of topologies to access
     */
    template<int dim_>
    std::vector<TopologyTemplateByDim<dim_>>& get_topologies()
    {
      using QueryShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      return RefinementTemplate<Shape_, QueryShape>::topologies;
    }

    /**
     * \brief const TopologyTemplate accessor
     *
     * \tparam Dimension of topologies to access
     */
    template<int dim_>
    const std::vector<TopologyTemplateByDim<dim_>>& get_topologies() const
    {
      using QueryShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      return RefinementTemplate<Shape_, QueryShape>::topologies;
    }

    /**
     * \brief Returns the number of entities of dimension \c dim_ in this template
     *
     * \tparam dim_ Dimension to query
     */
    template<int dim_>
    Index num_entities() const
    {
      if constexpr(dim_ > Shape_::dimension)
      {
        return 0;
      }
      else if constexpr(dim_ > 0 && dim_ <= Shape_::dimension)
      {
        using QueryShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
        return RefinementTemplate<Shape_, QueryShape>::topologies.size();
      }
      else
      {
        using QueryShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
        return RefinementTemplate<Shape_, QueryShape>::vertices.size();
      }
    }
  };

  /**
   * \brief Template for refining a mesh entity.
   *
   * This is a collection of TopologyTemplates, which describe how to construct
   * all children of all dimension of the template.
   *
   * \tparam Shape_ Shape of this template
   * \tparam CoShape_ Current shape of this data structure definition.
   *
   * This is the base case for the RefinementTemplate definition.
   */
  template<typename Shape_>
  struct RefinementTemplate<Shape_, Shape::Vertex>
  {
    /// Number of coefficients for creating a vertex of this template via linear interpolation
    static constexpr int num_coefficients = Shape::FaceTraits<Shape_, 0>::count;

    /// Vertices added by this template.
    std::vector<Tiny::Vector<Real, Shape_::dimension>> vertices;
    /// Vertices added by this template. Precomputed as coefficients for linear interpolation.
    std::vector<std::array<Real, num_coefficients>> vertex_coefficients;

    /// Vertex-accessor
    std::vector<Tiny::Vector<Real, Shape_::dimension>>& get_vertices()
    {
      return vertices;
    }

    /// const Vertex-accessor
    const std::vector<Tiny::Vector<Real, Shape_::dimension>>& get_vertices() const
    {
      return vertices;
    }

    /// Vertex coefficient accessor
    std::vector<std::array<Real, num_coefficients>>& get_vertex_coefficients()
    {
      return vertex_coefficients;
    }

    /// const Vertex coefficient accessor
    const std::vector<std::array<Real, num_coefficients>>& get_vertex_coefficients() const
    {
      return vertex_coefficients;
    }

    /**
     * \brief Returns the number of vertices in this template
     *
     * \note Calling with any \dim_ other than 0 is a compile-error.
     */
    template<int dim_>
    Index num_entities() const
    {
      static_assert(dim_ == 0);

      return vertices.size();
    }
  };

  /**
   * \brief Mapping between different orientations of a refinement template
   *
   * During adaptive mesh refinements, certain entities are shared between
   * multiple entities of the mesh. To deal with these shared elements in an
   * orderly fashion, the adaptive mesh refines them standalone, starting with
   * edges and working upward. Adjacent entities can then retrieve the shared
   * entities during their own construction. Because the orientation of
   * entities can differ between the local orientation seen by a
   * hihger-dimensional entity and the entities actual orientation, we need to
   * correct for orientation when retrieving such a shared entity.
   *
   * Consider for example an edge shared between two faces
   *
   * 2---3       0
   * |   |       |
   * 0---1       1
   *
   * Face 1     Edge
   *
   * By the usual FEAT convention the right edge of face 1 runs from bottom to
   * top, that is aginst the actual direction of the edge. When retrieving a
   * child of the edge, the order of the children as seen by the face and as
   * seen by the edge are thus reversed. The OrientationMappings correct for
   * this discrepancy.
   *
   * A OrientationMapping represents a function between the local child-indices
   * of some shared entity and the actual child entities. The
   * OrientationMappingHolder and OrienationMappingWrapper use the
   * OrientationMapping to build the full Type x Orientation x Index -> Index
   * mapping required for orientation-correction.
   */
  class OrientationMapping
  {
    /// Index -> Index mapping
    std::vector<std::pair<Index, int>> _mapping;

  public:
    /// Apply operator
    std::pair<Index, int> operator()(Index idx) const
    {
      return _mapping.at(idx);
    }

    /// Accessor for full orientation mapping
    std::vector<std::pair<Index, int>>& get_mapping()
    {
      return _mapping;
    }

    /// Const accessor for full orientation mapping
    const std::vector<std::pair<Index, int>>& get_mapping() const
    {
      return _mapping;
    }
  };

  /**
   * \brief Dimension-recursive container for OrientationMappings
   *
   * See OrientationMapping for details.
   */
  template<typename RawData_, typename Shape_, int dim_ = Shape_::dimension>
  struct OrientationMappingWrapper : public OrientationMappingWrapper<RawData_, Shape_, dim_ - 1>
  {
    static constexpr int num_orientations = Intern::orientations<Shape_>();

    template<typename S>
    using RefinementTypeByShape = typename RawData_::template RefinementTypeByShape<S>;

    std::unordered_map<RefinementTypeByShape<Shape_>, std::array<OrientationMapping, num_orientations>> mappings = {};

    template<int n_>
    OrientationMapping& get_mapping(RefinementTypeByShape<Shape_> type, int orientation)
    {
      static_assert(n_ <= dim_);

      ASSERT(orientation < num_orientations);

      auto& mappings_n = OrientationMappingWrapper<RawData_, Shape_, n_>::mappings;
      if(mappings_n.find(type) == mappings_n.end())
      {
        mappings_n[type] = {};
      }

      return mappings_n.at(type)[orientation];
    }

    template<int n_>
    const OrientationMapping& get_mapping(RefinementTypeByShape<Shape_> type, int orientation) const
    {
      static_assert(n_ <= dim_);

      ASSERT(mappings.find(type) != mappings.end());
      ASSERT(orientation < num_orientations);

      auto& mappings_n = OrientationMappingWrapper<RawData_, Shape_, n_>::mappings;
      return mappings_n.at(type)[orientation];
    }
  };

  /**
   * \brief Dimension-recursive container for OrientationMappings
   *
   * See OrientationMapping for details.
   * Base case.
   */
  template<typename RawData_, typename Shape_>
  struct OrientationMappingWrapper<RawData_, Shape_, 0>
  {
    static constexpr int num_orientations = Intern::orientations<Shape_>();

    template<typename S>
    using RefinementTypeByShape = typename RawData_::template RefinementTypeByShape<S>;

    std::unordered_map<RefinementTypeByShape<Shape_>, std::array<OrientationMapping, num_orientations>> mappings = {};

    template<int n_>
    OrientationMapping& get_mapping(RefinementTypeByShape<Shape_> type, int orientation)
    {
      static_assert(n_ == 0);

      ASSERT(mappings.find(type) != mappings.end());
      ASSERT(orientation < num_orientations);

      return OrientationMappingWrapper<RawData_, Shape_, n_>::mappings.at(type)[orientation];
    }

    template<int n_>
    const OrientationMapping& get_mapping(RefinementTypeByShape<Shape_> type, int orientation) const
    {
      static_assert(n_ == 0);

      ASSERT(mappings.find(type) != mappings.end());
      ASSERT(orientation < num_orientations);

      return OrientationMappingWrapper<RawData_, Shape_, n_>::mappings.at(type)[orientation];
    }
  };

  /**
   * \brief Dimension-recursive container for OrientationMappingWrappers
   *
   * See OrientationMapping for details.
   */
  template<typename RawData_, typename Shape_>
  struct OrientationMappingHolder
    : OrientationMappingHolder<RawData_, typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType>
  {
    OrientationMappingWrapper<RawData_, Shape_> wrapper = {};

    template<int dim_>
    using RefinementTypeByDim = typename RawData_::template RefinementTypeByDim<dim_>;

    template<int dim_, int codim_>
    OrientationMapping&
    get_mapping(RefinementTypeByDim<dim_> type, int orientation)
    {
      using DimShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      return OrientationMappingHolder<RawData_, DimShape>::wrapper.template get_mapping<codim_>(
        type,
        orientation);
    }

    template<int dim_, int codim_>
    const OrientationMapping&
    get_mapping(RefinementTypeByDim<dim_> type, int orientation) const
    {
      using DimShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      return OrientationMappingHolder<RawData_, DimShape>::wrapper.template get_mapping<codim_>(
        type,
        orientation);
    }
  };

  /**
   * \brief Dimension-recursive container for OrientationMappingWrappers
   *
   * See OrientationMapping for details.
   * Base case.
   */
  template<typename RawData_>
  struct OrientationMappingHolder<RawData_, Shape::Vertex>
  {
  };

  /**
   * \brief Tuple of raw entity and corresponding reference
   *
   * \tparam Shape_ Shape of the raw entity
   * \tparam num_coords_ Number of coordinates for each vertex of the raw entity
   */
  template<typename Shape_, int num_coords_>
  struct EntitySearchEntry
  {
    /// Raw entity
    RawEntity<Shape_, num_coords_> raw_entity;

    /// EntityReference refering to the raw entity.
    EntityReference reference;

    /// Constructor
    EntitySearchEntry(RawEntity<Shape_, num_coords_> ent, EntityReference ref) : raw_entity(ent), reference(ref)
    {
    }
  };

  /// Output-operator for EntitySearchEntry
  template<typename Shape_, int num_coords_>
  std::ostream& operator<<(std::ostream& stream, const EntitySearchEntry<Shape_, num_coords_>& entry)
  {
    stream << "EntitySearchEntry<" << Shape_::name() << ", " << num_coords_ << "> { raw_entity: " << entry.raw_entity
           << ", reference: " << entry.reference << " }";
    return stream;
  }

  /**
   * \brief Mapping between RawEntities and EntityReferences
   *
   * \tparam Shape_ Template shape
   * \tparam num_coords_ Number of coordinates for each vertex
   *
   * The TemplateBuilder constructs RefinementTemplates from lists of
   * RawEntities. During this process RawEntities must be mapped to
   * EntityReferences, so that the AdaptiveMesh can retrieve the correct
   * entities during refinement.
   *
   * This class manages a list of EntitySearchEntries to allow this mapping.
   * It is dimension-recursive to allows storing entities of multiple dimensions.
   */
  template<typename Shape_, int num_coords_ = Shape_::dimension>
  class TemplateSearchSpace
    : public TemplateSearchSpace<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType, num_coords_>
  {
  public:
    /// Current shape
    using ShapeType = Shape_;

    /// Type of search entry for current shape
    using Entry = EntitySearchEntry<Shape_, num_coords_>;
    static const constexpr int num_coords = num_coords_;

    /// Our base class
    using BaseClass =
      TemplateSearchSpace<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType, num_coords_>;

    /// Accessor for EntitySearchEntry types
    template<int dim_>
    using EntryByDim = EntitySearchEntry<typename Shape::FaceTraits<Shape_, dim_>::ShapeType, num_coords_>;

    /// Accessor for RawTemplate types
    template<int dim_>
    using RawTemplateByDimes = RawTemplate<typename Shape::FaceTraits<ShapeType, dim_>::ShapeType>;

    /// Accessor for search space entries
    template<int dim_>
    std::vector<EntryByDim<dim_>>& entries()
    {
      static_assert(dim_ <= Shape_::dimension);
      static_assert(dim_ >= 0);

      using DimShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      return TemplateSearchSpace<DimShape, num_coords_>::_entries;
    }

    /// Const accessor for search space entries
    template<int dim_>
    const std::vector<EntryByDim<dim_>>& entries() const
    {
      static_assert(dim_ <= Shape_::dimension);
      static_assert(dim_ >= 0);

      using DimShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
      return TemplateSearchSpace<DimShape, num_coords_>::_entries;
    }

    /**
     * \brief Search for a RawEntity in this search space
     *
     * \param[in] entity A raw entity
     *
     * \note O(n) runtime, n number of search space entries.
     *
     * \returns A EntityReference that refers to the given raw entity if one
     * exists. Aborts otherwise.
     */
    template<typename Shape__>
    EntityReference search(const RawEntity<Shape__, num_coords_>& entity) const
    {
      return TemplateSearchSpace<Shape__, num_coords>::_search(entity);
    }

    /**
     * \brief Add a sibling entity to the search space
     *
     * \param[in] entity The entity to add
     * \param[in] idx The index of the entity
     *
     * These sibling entities get created by the template builder during
     * template construction. We hence present a way to easily add these to the
     * search space.
     */
    template<typename Shape__>
    void add_sibling(const RawEntity<Shape__, num_coords_>& entity, Index idx)
    {
      return TemplateSearchSpace<Shape__, num_coords>::_add_sibling(entity, idx);
    }

    /**
     * \brief Adds boundary entities to the search space
     *
     * This function adds the entities of a template to the search space as
     * boundary elements.
     *
     * \param[in] parent_index Which boundary entity we are adding
     * \param[in] tmplt Entities to add
     */
    template<typename EntityShape_, int dim_ = EntityShape_::dimension>
    void add_boundary_entity(int parent_index, const RawTemplate<EntityShape_>& tmplt)
    {
      static_assert(EntityShape_::dimension <= ShapeType::dimension);

      using DimShape = typename Shape::FaceTraits<EntityShape_, dim_>::ShapeType;
      using Embedder = Intern::VertexEmbedder<ShapeType, EntityShape_::dimension>;

      static constexpr int num_vertices = Shape::FaceTraits<DimShape, 0>::count;

      EntitySource source = EntityShape_::dimension == 1 ? EntitySource::BoundaryEdge : EntitySource::BoundaryFace;

      if constexpr(dim_ == 0)
      {
        auto& vertices = this->template entries<0>();
        Index vertex_index = 0;
        for(const auto& vertex : Intern::internal_vertices(tmplt))
        {
          RawEntity<Shape::Vertex, num_coords> entity(Embedder::embed(parent_index, vertex));
          vertices.emplace_back(entity, EntityReference{source, vertex_index++, 0, parent_index});
        }
      }
      else
      {
        add_boundary_entity<EntityShape_, dim_ - 1>(parent_index, tmplt);

        auto& entries = this->template entries<dim_>();

        // Add internal elements of dimension dim_ to search space
        Index idx = 0;
        for(const auto& entity : Intern::internal_entities<EntityShape_, dim_>(tmplt.entities))
        {
          RawEntity<DimShape, num_coords> embedded;

          for(int i(0); i < num_vertices; i++)
          {
            embedded.coords[i] = Embedder::embed(parent_index, entity.coords[i]);
          }
          entries.emplace_back(embedded, EntityReference{source, idx++, 0, parent_index});
        }
      }
    }

    /**
     * \brief Adds topology entities to a TemplateSearchSpace
     *
     * This function adds the entities of a shapes topology to a search space.
     *
     * \param[in] search_space The search space to add to
     */
    template<int dim_ = ShapeType::dimension - 1>
    void add_topology()
    {
      static constexpr int num_vertices = Shape::FaceTraits<ShapeType, 0>::count;

      std::array<Tiny::Vector<Real, num_coords>, num_vertices> vertices = Intern::topology_vertices<ShapeType>();
      _add_topology_rec(vertices);
    }


  protected:
    /// Helper-function for add_topology_to_search_space
    template<int dim_ = ShapeType::dimension - 1>
    void _add_topology_rec(
      const std::array<Tiny::Vector<Real, num_coords>, Shape::FaceTraits<ShapeType, 0>::count>& vertices)
    {
      if constexpr(dim_ == 0)
      {
        static constexpr int num_vertices = Shape::FaceTraits<ShapeType, 0>::count;

        auto& entries = this->template entries<0>();

        for(int vertex(0); vertex < num_vertices; vertex++)
        {
          RawEntity<Shape::Vertex, num_coords> entity;
          entity.coords[0] = vertices[vertex];
          entries.emplace_back(entity, EntityReference{EntitySource::ParentTopology, static_cast<Index>(vertex), 0, 0});
        }
      }
      else
      {
        using FaceShape = typename Shape::FaceTraits<ShapeType, dim_>::ShapeType;
        using FaceMapping = Intern::FaceIndexMapping<ShapeType, dim_, 0>;

        static constexpr int num_faces = Shape::FaceTraits<ShapeType, dim_>::count;
        static constexpr int num_vertices = Shape::FaceTraits<FaceShape, 0>::count;

        auto& entries = this->template entries<dim_>();
        for(int face(0); face < num_faces; face++)
        {
          RawEntity<FaceShape, num_coords> entity;

          for(int vertex(0); vertex < num_vertices; vertex++)
          {
            entity.coords[vertex] = vertices[FaceMapping::map(face, vertex)];
          }
          entries.emplace_back(entity, EntityReference{EntitySource::ParentTopology, static_cast<Index>(face), 0, 0});
        }

        _add_topology_rec<dim_ - 1>(vertices);
      }
    }

    /// Inner search logic for the current shape
    EntityReference _search(const RawEntity<Shape_, num_coords_>& entity) const
    {
      auto comp = [&](auto entry) { return Intern::is_same_raw_entity(entry.raw_entity, entity); };
      auto it = std::find_if(_entries.begin(), _entries.end(), comp);

      if(it == _entries.end())
      {
        std::cout << "Searching for entity " << entity << "\n";
        XABORTM("Entity not in search space");
      }

      // We found a valid reference. Next we need to determine the orientation of entity,
      // as compared to the found reference.
      // Because the CongurencySampler can not directly compare vertices, we produce index-lists here.
      using CongruencySampler = Intern::CongruencySampler<Shape_>;

      int orientation = CongruencySampler::compare_with(entity.coords.data(), it->raw_entity.coords.data(), Intern::is_same_vertex<num_coords_>);
      EntityReference ref = it->reference;

      ref.orientation = orientation;

      if constexpr(num_coords_ < 3)
      {
        // 3D or less templates should not produce boundary face references
        ASSERT(ref.source != EntitySource::BoundaryFace);
      }

      return ref;
    }

    /// Inner add_sibling logic for the current shape.
    void _add_sibling(const RawEntity<Shape_, num_coords_>& entity, Index idx)
    {
      _entries.push_back({entity, EntityReference{EntitySource::Sibling, idx, 0, 0}});
    }

    /// Vector of search space entries.
    std::vector<Entry> _entries;
  };

  /**
   * \brief Mapping between RawEntities and EntityReferences
   *
   * \tparam num_coords_ Number of coordinates for each vertex
   *
   * Base case. See the general case for details.
   */
  template<int num_coords_>
  class TemplateSearchSpace<Shape::Vertex, num_coords_>
  {
  public:
    /// Current shape
    using ShapeType = Shape::Vertex;

    /// Search space entry type
    using Entry = EntitySearchEntry<Shape::Vertex, num_coords_>;

    /// Number of coordinates for vertices
    static const constexpr int num_coords = num_coords_;

    /// Accessor for search space entries
    template<int dim_>
    std::vector<Entry>& entries()
    {
      static_assert(dim_ == 0);
      return _entries;
    }

    /// Const accessor for search space entries
    template<int dim_>
    const std::vector<Entry>& entries() const
    {
      static_assert(dim_ == 0);

      return _entries;
    }

    /**
     * \brief Search for a vertex in this search space
     *
     * \param[in] vertx Vertex to search for
     *
     * \note O(n) runtime, n number of search space entries.
     *
     * \returns A EntityReference that refers to the given vertex if one
     * exists. Aborts otherwise.
     */
    EntityReference search_vertex(const Tiny::Vector<Real, num_coords>& vertex) const
    {
      for(const auto& entry : _entries)
      {
        if(Intern::is_same_vertex(vertex, entry.raw_entity.coords[0]))
        {
          return entry.reference;
        }
      }
      std::cout << "Vertex: " << vertex << " not found\n";
      XABORTM("Vertex not in search space");
    }

    /**
     * \brief Add a sibling vertex to the search space
     *
     * \param[in] vertex The vertex to add
     * \param[in] idx The index of the entity
     *
     * These sibling vertices get created by the template builder during
     * template construction. We hence present a way to easily add these to the
     * search space.
     */
    void add_sibling_vertex(const Tiny::Vector<Real, num_coords>& vertex, Index idx)
    {
      _entries.push_back(
        {RawEntity<Shape::Vertex, num_coords>(vertex), EntityReference{EntitySource::Sibling, idx, 0, 0}});
    }

    /**
     * \brief Search for a vertex given a EntityReference
     *
     * \param[in] reference Entity reference to find vertex for
     *
     * \returns The vertex refered to by \c reference. Aborts otherwise.
     */
    Tiny::Vector<Real, num_coords_> vertex_by_reference(const EntityReference& reference) const
    {
      for(const auto& entry : _entries)
      {
        if(entry.reference == reference)
        {
          return entry.raw_entity.coords[0];
        }
      }
      XABORTM("EntifyReference not in search space!");
    }

    /// Vector of search space entries
    std::vector<Entry> _entries;
  };


  /**
   * \brief Constructs RefinementTemplates from RawTemplates
   *
   * \tparam RawData_ Raw template data. See Geometry::SchneidersData for an example.
   *
   * The TemplateBuilder is lazy. Templates will be constructed the first time
   * they are accessed. This especially means hexahedral templates will not be
   * constructed for 2D meshes and defined but unused template sets inccur no
   * overhead.
   */
  template<class RawData_>
  class TemplateBuilder
  {
    /// Maximum shape supported by the raw data
    using MaxShape = typename RawData_::MaxShape;

    /// Accessor for raw template maps
    template<int dim_>
    using RawTemplateMapByDim = typename RawData_::template TemplateMapByDim<dim_>;

    /// Map type from refinement types to raw edge templates
    using RawEdgeMap = typename RawData_::EdgeMap;
    /// Map type from refinement types to raw face templates
    using RawFaceMap = typename RawData_::FaceMap;
    /// Map type from refinement types to raw cell templates
    using RawCellMap = typename RawData_::CellMap;

    /// Accessor for refinement types
    template<int dim_>
    using RefinementTypeByDim = typename RawData_::template RefinementTypeByDim<dim_>;

    /// Accessor for refinement template types
    template<int dim_>
    using RefinementTemplateByDim = RefinementTemplate<typename Shape::FaceTraits<MaxShape, dim_>::ShapeType>;

    /// Accessor for template map types
    template<int dim_>
    using TemplateMapByDim = std::unordered_map<RefinementTypeByDim<dim_>, RefinementTemplate<Shape::Hypercube<dim_>>>;

    /// Map type from refinement types to edge templates
    using EdgeMap = TemplateMapByDim<1>;
    /// Map type from refinement types to face templates
    using FaceMap = TemplateMapByDim<2>;
    /// Map type from refinement types to cell templates
    using CellMap = TemplateMapByDim<3>;

    /// Edge-templates constructed by this template builder
    EdgeMap _edge_templates = {};
    /// Face-templates constructed by this template builder
    FaceMap _face_templates = {};
    /// Cell-templates constructed by this template builder
    CellMap _cell_templates = {};

    /// Container for orientation mappings
    OrientationMappingHolder<RawData_, Shape::Hexahedron> _orientation_mappings = {};

    /// Type adjustments for edges
    std::array<std::array<std::uint64_t, 2>, 4> _edge_type_adjustments = {};
    /// Type adjustments for faces
    std::array<std::array<std::uint64_t, 4>, 16> _face_type_adjustments = {};
    /// Type adjustments for cells
    std::array<std::array<std::uint64_t, 8>, 256> _cell_type_adjustments = {};

    /// Bookkeeping for which type adjustments have been created already.
    std::array<bool, 4> _type_adjustment_initialized = {false, false, false, false};

    /// Edge-template accessor
    const EdgeMap& edges()
    {
      if(_edge_templates.empty())
      {
        _edge_templates = _build<1>(RawData_::raw_edges());
      }
      return _edge_templates;
    }

    /// Face-template accessor
    const FaceMap& faces()
    {
      if(_face_templates.empty())
      {
        _face_templates = _build<2>(RawData_::raw_faces());
      }
      return _face_templates;
    }

    /// Cell-template accessor
    const CellMap& cells()
    {
      if(_cell_templates.empty())
      {
        _cell_templates = _build<3>(RawData_::raw_cells());
      }
      return _cell_templates;
    }

  public:
    template<int dim_>
    const RefinementTemplateByDim<dim_>& get_template(RefinementTypeByDim<dim_> type)
    {
      if constexpr(dim_ == 1)
      {
        return edges().at(type);
      }
      else if constexpr(dim_ == 2)
      {
        return faces().at(type);
      }
      else if constexpr(dim_ == 3)
      {
        return cells().at(type);
      }
      XABORTM("Unsupported template dimension!");
    }

    template<int dim_>
    bool has_template(RefinementTypeByDim<dim_> type)
    {
      if constexpr(dim_ == 1)
      {
        return edges().find(type) != edges().end();
      }
      else if constexpr(dim_ == 2)
      {
        return faces().find(type) != faces().end();
      }
      else if constexpr(dim_ == 3)
      {
        return cells().find(type) != cells().end();
      }
      XABORTM("Unsupported template dimension!");
    }

    /**
     * \brief Correct an index for orientation
     */
    template<int dim_, int codim_>
    std::pair<Index, int> correct_for_orientation(RefinementTypeByDim<dim_> type, int orientation, Index idx)
    {
      return _orientation_mappings.template get_mapping<dim_, codim_>(type, orientation)(idx);
    }

    /*
     * \brief Accessor for type adjustments.
     *
     * \tparam dim_ Dimension to retrieve adjustments for
     *
     * \param[in] type Current markings of a mesh element
     *
     * \returns An array containing a delta of subdivision levels.
     *
     * A template set does not necessarily contain templates for all possible
     * markings of a mesh element. In that case the AdaptiveMesh corrects the
     * markings during the creation of the final refinement field. This method
     * supplies the necessary corrections.
     *
     * \note The adjustments will only ever add subdivision levels.
     */
    template<int dim_>
    const std::array<std::uint64_t, Shape::FaceTraits<Shape::Hypercube<dim_>, 0>::count>&
    type_adjustment(Index type)
    {
      if constexpr(dim_ == 1)
      {
        if(!_type_adjustment_initialized[1])
        {
          _create_type_adjustments<Shape::Hypercube<dim_>>(edges());
        }
        return _edge_type_adjustments.at(type);
      }
      else if constexpr(dim_ == 2)
      {
        if(!_type_adjustment_initialized[2])
        {
          _create_type_adjustments<Shape::Hypercube<dim_>>(faces());
        }
        return _face_type_adjustments.at(type);
      }
      else if constexpr(dim_ == 3)
      {
        if(!_type_adjustment_initialized[3])
        {
          _create_type_adjustments<Shape::Hypercube<dim_>>(cells());
        }
        return _cell_type_adjustments.at(type);
      }
      XABORTM("Unsupported template dimension!");
    }

    /// print some stats about maximum children in the raw data
    void stats()
    {
      std::array<Index, 4> maxima = {0, 0, 0, 0};

      for(const auto& entry : edges())
      {
        maxima[0] = Math::max(maxima[0], entry.second.template num_entities<0>());
        maxima[1] = Math::max(maxima[1], entry.second.template num_entities<1>());
      }
      std::cout << "Edges:\n";
      std::cout << "  Max Vertices: " << stringify(unsigned(maxima[0])) << "\n";
      std::cout << "  Max Edges: " << stringify(unsigned(maxima[1])) << "\n";

      maxima = {0, 0, 0, 0};
      for(const auto& entry : faces())
      {
        maxima[0] = Math::max(maxima[0], entry.second.template num_entities<0>());
        maxima[1] = Math::max(maxima[1], entry.second.template num_entities<1>());
        maxima[2] = Math::max(maxima[2], entry.second.template num_entities<2>());
      }
      std::cout << "Faces:\n";
      std::cout << "  Max Vertices: " << stringify(unsigned(maxima[0])) << "\n";
      std::cout << "  Max Edges: " << stringify(unsigned(maxima[1])) << "\n";
      std::cout << "  Max Faces: " << stringify(unsigned(maxima[2])) << "\n";

      maxima = {0, 0, 0, 0};
      for(const auto& entry : cells())
      {
        maxima[0] = Math::max(maxima[0], entry.second.template num_entities<0>());
        maxima[1] = Math::max(maxima[1], entry.second.template num_entities<1>());
        maxima[2] = Math::max(maxima[2], entry.second.template num_entities<2>());
        maxima[3] = Math::max(maxima[3], entry.second.template num_entities<3>());
      }
      std::cout << "Faces:\n";
      std::cout << "  Max Vertices: " << stringify(unsigned(maxima[0])) << "\n";
      std::cout << "  Max Edges: " << stringify(unsigned(maxima[1])) << "\n";
      std::cout << "  Max Faces: " << stringify(unsigned(maxima[2])) << "\n";
      std::cout << "  Max Cells: " << stringify(unsigned(maxima[3])) << "\n";
    }

  private:
    /// Determines the index of a vertex in a list of vertices
    template<int n_>
    Index _find_vertex(const std::vector<Tiny::Vector<Real, n_>>& vertices, const Tiny::Vector<Real, n_>& vertex)
    {
      auto pred = [&](auto v) { return Intern::is_same_vertex(vertex, v); };
      return std::distance(vertices.begin(), std::find_if(vertices.begin(), vertices.end(), pred));
    }

    /// dim-based accessor for type adjustment arrays
    template<int dim_>
    auto& type_adjustments()
    {
      if constexpr(dim_ == 1)
      {
        return _edge_type_adjustments;
      }
      else if constexpr(dim_ == 2)
      {
        return _face_type_adjustments;
      }
      else if constexpr(dim_ == 3)
      {
        return _cell_type_adjustments;
      }
      XABORTM("Unsupported template dimension!");
    }

    /**
     * \brief Constructs a TopologyTemplate for the given RawEntity
     *
     * \tparam TemplateShape_ Shape of the template the topology is part of
     * \tparam TopologyShape_ Shape of the entity the topology is for
     * \tparam dim_ Dimension the algorithm is currently working on
     *
     * \param[in] raw_entity The raw entity to construct a topology template for
     * \param[inout] topo The topology being constructed
     * \param[in] search_space The TemplateSearchSpace for the current template.
     */
    template<typename TemplateShape_, typename TopologyShape_, int dim_ = TopologyShape_::dimension - 1>
    void _build_topology(
      const RawEntity<TopologyShape_, TemplateShape_::dimension>& raw_entity,
      const TemplateSearchSpace<TemplateShape_>& search_space,
      TopologyTemplate<TopologyShape_>& topo)
    {
      static constexpr int num_entities = Shape::FaceTraits<TopologyShape_, dim_>::count;
      if constexpr(dim_ == 0)
      {
        auto& vertex_refs = topo.template get_references<0>();
        for(int i(0); i < num_entities; i++)
        {
          vertex_refs[i] = search_space.search_vertex(raw_entity.coords[i]);
        }
      }
      else
      {
        auto& refs = topo.template get_references<dim_>();
        for(int i(0); i < num_entities; i++)
        {
          refs[i] = search_space.search(raw_entity.template face<dim_>(i));
        }

        _build_topology<TemplateShape_, TopologyShape_, dim_ - 1>(raw_entity, search_space, topo);
      }
    }

    /**
     * \brief Constructs a RefinementTemplate from the given RawTemplate
     *
     * \tparam TemplateShape_ Shape of the template
     * \tparam dim_ Dimension the algorithm is currently working on
     *
     * \param[in] raw_template The raw template to construct a refinement template for
     * \param[inout] tmplt The refinement template being constructed
     * \param[inout] search_space The TemplateSearchSpace for the current template.
     */
    template<typename TemplateShape_, int dim_ = TemplateShape_::dimension>
    void _build_template(
      const RawTemplate<TemplateShape_>& raw_template,
      RefinementTemplate<TemplateShape_>& tmplt,
      TemplateSearchSpace<TemplateShape_>& search_space)
    {
      if constexpr(dim_ == 0)
      {
        for(const auto& vertex : Intern::internal_vertices(raw_template))
        {
          // Add vertex and its coefficients to template
          tmplt.get_vertices().push_back(vertex);
          tmplt.get_vertex_coefficients().push_back(Intern::vertex_coefficients<TemplateShape_>(vertex));

          // Add entry to search space for building topologies later
          search_space.add_sibling_vertex(vertex, tmplt.get_vertices().size() - 1);
        }
      }
      else
      {
        _build_template<TemplateShape_, dim_ - 1>(raw_template, tmplt, search_space);

        using TopologyShape = typename Shape::FaceTraits<TemplateShape_, dim_>::ShapeType;

        auto& topologies = tmplt.template get_topologies<dim_>();
        for(const auto& entity : Intern::internal_entities<TemplateShape_, dim_>(raw_template.entities))
        {
          TopologyTemplate<TopologyShape> topology{};
          _build_topology<TemplateShape_, TopologyShape>(entity, search_space, topology);
          topologies.push_back(topology);

          search_space.add_sibling(entity, topologies.size() - 1);
        }
      }
    }

    /**
     * \brief Fill orientation mappings for a given template
     *
     * \tparam TemplateShape_ Shape of the template
     * \tparam dim_ Current working dim
     *
     * \param[in] oriented_type The refinement type of the given oriented template
     * \param[in] orientation The orientation to fill the mappings of
     * \param[in] local_template The refinement template, as seen by the current mesh entity
     * \parma[in] local_search_space Search space containing all elements of local_template
     * \param[in] oriented_template Refinement template as seen by sub-entity with orientation \c orientation
     * \param[in] oriented_search_space Search space containing all elements of oriented_template
     *
     * Fills orientation mapping for (oriented_type, orientation) as a side effect.
     *
     * To compute these mappings, we use the following idea.
     *
     * Let e be a mesh entity with dimension >= 2, i.e. at least a face, and let
     * s be a sub-entity of e, i.e. an edge or face of e's topology.
     *
     * The sub-entity s then has two refinement types. A local type, as
     * determined by the vertex-ordering of e, and a foreign tye, as determined
     * by the vertex ordering of s itself. We can predict the foreign type
     * using the local type and the relative orientation of s as stored in the
     * entities topology.
     *
     * This allows us to fill the orientation mappings by finding the foreign
     * type for each pair of local type and possible orientation. Given a local
     * template and the corresponding oriented template, as are passed to this
     * function, we can retrieve the vertices of all entities from their
     * respective search spaces, and determine the orientation mapping by
     * finding identical entities in both templates.
     *
     * Because the templates determine how mesh entities will be constructed in
     * the adaptive mesh, we can also pre-compute these entities relative
     * orientations to each other at this point. This means no orientations
     * have to be computed by the adaptive mesh beyond the intial collection of
     * root topologies.
     */
    template<typename TemplateShape_, int dim_ = TemplateShape_::dimension>
    void _fill_orientation_mappings(
      const RefinementTypeByDim<TemplateShape_::dimension>& oriented_type,
      const int orientation,
      const RefinementTemplate<TemplateShape_>& local_template,
      const TemplateSearchSpace<TemplateShape_>& local_search_space,
      const RefinementTemplate<TemplateShape_>& oriented_template,
      const TemplateSearchSpace<TemplateShape_>& oriented_search_space)
    {

      using VertexType = Tiny::Vector<Real, TemplateShape_::dimension>;

      if constexpr(dim_ == 0)
      {
        OrientationMapping& vertex_mapping =
          _orientation_mappings.template get_mapping<TemplateShape_::dimension, 0>(oriented_type, orientation);

        // For each vertex of the local template
        // * map it into the oriented frame of reference
        // * search for the corresponding vertex of the oriented template
        // * save the found index in the orientation mapping
        for(const VertexType& local_vertex : local_template.get_vertices())
        {
          VertexType oriented_vertex = Intern::orient_vertex<TemplateShape_>(local_vertex, orientation);
          Index mapped = _find_vertex(oriented_template.get_vertices(), oriented_vertex);
          vertex_mapping.get_mapping().emplace_back(mapped, 0);
        }
      }
      else
      {
        _fill_orientation_mappings<TemplateShape_, dim_ - 1>(
          oriented_type,
          orientation,
          local_template,
          local_search_space,
          oriented_template,
          oriented_search_space);

        using CurrentShape = typename Shape::FaceTraits<TemplateShape_, dim_>::ShapeType;
        static constexpr int num_vertices = Shape::FaceTraits<CurrentShape, 0>::count;

        OrientationMapping& mapping =
          _orientation_mappings.template get_mapping<TemplateShape_::dimension, dim_>(oriented_type, orientation);

        for(const auto& local_topo : local_template.template get_topologies<dim_>())
        {
          // Collect local vertices from the search space
          std::array<VertexType, num_vertices> local = {};

          const auto& local_refs = local_topo.template get_references<0>();
          for(int i(0); i < num_vertices; i++)
          {
            local[i] =
              Intern::orient_vertex<TemplateShape_>(local_search_space.vertex_by_reference(local_refs[i]), orientation);
          }

          Index entity_idx = 0;
          for(const auto& oriented_topo : oriented_template.template get_topologies<dim_>())
          {
            // Collect oriented vertices from the oriented search space
            std::array<VertexType, num_vertices> oriented = {};

            const auto& oriented_refs = oriented_topo.template get_references<0>();
            for(int i(0); i < num_vertices; i++)
            {
              oriented[i] = oriented_search_space.vertex_by_reference(oriented_refs[i]);
            }

            // Compare local and oriented entities
            RawEntity<CurrentShape, TemplateShape_::dimension> ent_local;
            ent_local.coords = local;

            RawEntity<CurrentShape, TemplateShape_::dimension> ent_oriented;
            ent_oriented.coords = oriented;

            if(Intern::is_same_raw_entity(ent_local, ent_oriented))
            {
              // The entities are identical. This means local_topo and
              // oriented_topo are the same mesh element and we can save this
              // mapping in the orientation mapping. Additionally we can
              // pre-compute their relative orientation at this point, because
              // we know how both of these entities will be constructed in the
              // final mesh and how their parents are oriented relative to each
              // other.
              using CongruencySampler = Intern::CongruencySampler<CurrentShape>;
              int o = CongruencySampler::compare_with(
                ent_local.coords.data(),
                ent_oriented.coords.data(),
                Intern::is_same_vertex<TemplateShape_::dimension>);

              mapping.get_mapping().emplace_back(entity_idx, o);
              break;
            }
            entity_idx++;
          }
        }
      }
    }

    /**
     * \brief Build all templates of a certain dimension
     *
     * \tparam dim_ Dimension to build templates for
     *
     * \param[in] raw_templates Map of raw templates
     *
     * \returns A map of refinement templates
     */
    template<int dim_>
    TemplateMapByDim<dim_> _build(const RawTemplateMapByDim<dim_>& raw_templates)
    {
      using TemplateShape = typename Shape::FaceTraits<Shape::Hypercube<dim_>, dim_>::ShapeType;

      TemplateMapByDim<dim_> result;

      // We require the full template search spaces to create orientation mappings with.
      std::unordered_map<RefinementTypeByDim<dim_>, TemplateSearchSpace<TemplateShape>> search_spaces;

      for(const auto& entry : raw_templates)
      {
        RefinementTemplate<Shape::Hypercube<dim_>> tmplt;
        TemplateSearchSpace<TemplateShape> search_space = _setup_search_space<Shape::Hypercube<dim_>>(entry.first);

        _build_template(entry.second, tmplt, search_space);

        result.insert({entry.first, tmplt});
        search_spaces.insert({entry.first, search_space});
      }

      // We have now created all templates for this dimension.
      // Next, we create the corresponding orientation mappings.

      // NOTE(mmuegge): _fill_orientation_mappings is fully dimension-generic,
      // apart from CongruencySamplers for higher dimensions. We thus disable
      // orientation mapping creation for all dimensions above 2
      if constexpr(dim_ <= 2)
      {
        for(const auto& entry : raw_templates)
        {
          RefinementTypeByDim<dim_> local_type = entry.first;
          auto& search_space = search_spaces[local_type];

          for(int orientation = 0; orientation < Intern::orientations<Shape::Hypercube<dim_>>(); orientation++)
          {
            const RefinementTypeByDim<dim_> oriented_type = local_type.orient(orientation);
            auto& oriented_search_space = search_spaces[oriented_type];

            const auto& local_tmplt = result.at(local_type);
            const auto& oriented_tmplt = result.at(oriented_type);

            _fill_orientation_mappings<TemplateShape>(
              oriented_type,
              orientation,
              local_tmplt,
              search_space,
              oriented_tmplt,
              oriented_search_space);
          }
        }
      }
      return result;
    }

    /**
     * \brief Helper function for creating a search space with boundary and topology entities
     *
     * \tparam Shape_ Shape of the search space
     *
     * \param type Refinement type for determining boundary entities
     *
     * \returns A TemplateSearchSpace with all necessary boundary and topology elements for creating a refinement
     * template for the given type.
     */
    template<typename Shape_>
    TemplateSearchSpace<Shape_> _setup_search_space(RefinementTypeByDim<Shape_::dimension> type) const
    {
      // Setup result.
      TemplateSearchSpace<Shape_> result = {};

      // For faces and cells we also have to add the boundary edges and internal edges
      if constexpr(Shape_::dimension > 1)
      {
        static const constexpr int num_edges = Shape::FaceTraits<Shape_, 1>::count;
        for(int i(0); i < num_edges; ++i)
        {
          auto edge_type = type.template face_type<1>(i);
          result.add_boundary_entity(i, RawData_::raw_edges().at(edge_type));
        }
      }

      // For cells, we have to add boundary faces and internal faces
      if constexpr(Shape_::dimension > 2)
      {
        static const constexpr int num_faces = Shape::FaceTraits<Shape_, 2>::count;
        for(int i(0); i < num_faces; ++i)
        {
          auto face_type = type.template face_type<2>(i);
          result.add_boundary_entity(i, RawData_::raw_faces().at(face_type));
        }
      }

      result.add_topology();

      return result;
    }

    /**
     * \brief Calculate adjustments of the subdivision levels
     *
     * \tparam Shape_ Shape to calculate adjustments for
     *
     * A template set does not necessarily contain templates for all possible
     * markings of a mesh element. In that case the AdaptiveMesh corrects the
     * markings during the creation of the final refinement field. This method
     * supplies the necessary corrections.
     *
     * It does so by determining the "closest" template, i.e. the template with
     * the least required changes to the subidivision levels, to all possible
     * markings and saving the neccessary changes to apply that template to the
     * element.
     */
    template<typename Shape_, typename TemplateMap_>
    void _create_type_adjustments(const TemplateMap_& available_templates)
    {
      auto& adjustments = type_adjustments<Shape_::dimension>();

      static const constexpr int num_verts = Shape::FaceTraits<Shape_, 0>::count;
      for(Index type(0); type < (1 << num_verts); ++type)
      {
        auto type_marking = VertexMarking<Shape_>(type);

        std::size_t closest_count = num_verts;
        auto closest_marking = VertexMarking<Shape_>::all_marked();

        for(auto& entry : available_templates)
        {
          VertexMarking<Shape_> entry_marking(entry.first.to_vertex_marking());

          // Is the current type a subtype of this type?
          if(!entry_marking.covers(type_marking))
          {
            continue;
          }
          std::size_t distance = type_marking.distance(entry_marking);

          if(distance < closest_count)
          {
            closest_count = distance;
            closest_marking = entry_marking;
          }
        }

        // Found the type that requires the least changes
        // Next, determine exact required changes

        auto& adjustment = adjustments.at(type);
        for(int i(0); i < num_verts; ++i)
        {
          if(closest_marking.is_vertex_marked(i) && !type_marking.is_vertex_marked(i))
          {
            adjustment[i] = 1;
          }
        }
      }
    }
  };
} // namespace FEAT::Geometry

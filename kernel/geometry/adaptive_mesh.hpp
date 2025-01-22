// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include "kernel/geometry/intern/adaptive_refinement_utils.hpp"
#include <kernel/geometry/templates/two_refinement_data.hpp>
#include <kernel/geometry/intern/congruency_mapping.hpp>
#include <kernel/geometry/intern/congruency_sampler.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/geometry/templates/schneiders_data.hpp>
#include <kernel/geometry/templates/sun_zhao_ma_data.hpp>
#include <kernel/geometry/templates/sun_zhao_ma_expansion_data.hpp>
#include <kernel/shape.hpp>
#include <kernel/base_header.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/intern/adaptive_mesh_algorithms.hpp>
#include <kernel/geometry/intern/adaptive_mesh_storage.hpp>
#include <kernel/geometry/intern/refinement_field.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/geometry/template_interface.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/geometry/template_sets.hpp>
#include <kernel/geometry/template_builder.hpp>

// includes, system
#include <array>
#include <limits>
#include <set>
#include <unordered_map>

namespace FEAT::Geometry
{
  /// Re-export of layer type
  using Layer = Intern::Layer;

  // Possible template sets
  using SchneidersTemplates = StandardTemplateSet<SchneidersData>;
  using SunZhaoMaTemplates = IsolatedPointTemplateSet<SunZhaoMaData>;
  using SunZhaoMaExpansionTemplates = IsolatedPointTemplateSet<SunZhaoMaExpansionData>;
  using TwoRefinementTemplates = StandardTemplateSet<TwoRefinementData>;

  namespace Intern
  {
    /**
     * \brief Mapping from ConformalMesh Indices to MeshStorage keys
     */
    template<int dim_>
    struct MeshRoots : public MeshRoots<dim_ - 1>
    {
      std::unordered_map<Index, Intern::ElementKey<dim_>> roots;

      template<int query_dim_>
      std::unordered_map<Index, Intern::ElementKey<query_dim_>>& by_dim()
      {
        static_assert(query_dim_ <= dim_);
        return MeshRoots<query_dim_>::roots;
      }

      template<int query_dim_>
      const std::unordered_map<Index, Intern::ElementKey<query_dim_>>& by_dim() const
      {
        static_assert(query_dim_ <= dim_);
        return MeshRoots<query_dim_>::roots;
      }
    };

    /**
     * \brief Mapping from ConformalMesh Indices to MeshStorage keys
     */
    template<>
    struct MeshRoots<0>
    {
      std::unordered_map<Index, Intern::ElementKey<0>> roots;

      template<int query_dim_>
      std::unordered_map<Index, Intern::ElementKey<0>>& by_dim()
      {
        static_assert(query_dim_ == 0);
        return roots;
      }

      template<int query_dim_>
      const std::unordered_map<Index, Intern::ElementKey<0>>& by_dim() const
      {
        static_assert(query_dim_ == 0);
        return roots;
      }
    };
  }

  /**
   * \brief Statistics about entities added, removed, and kept during a mesh adaption
   *
   * \author Markus Muegge
   */
  struct AdaptionStats
  {
    /**
     * \brief Statistics about added entities
     *
     * \c added[dim] contains the number of entities of dimension \c dim added to the mesh during adaption.
     */
    std::array<Index, 4> added;

    /**
     * \brief Statistics about removed entities
     *
     * \c added[dim] contains the number of entities of dimension \c dim removed from the mesh during adaption.
     */
    std::array<Index, 4> removed;

    /**
     * \brief Statistics about kept entities
     *
     * \c added[dim] contains the number of entities of dimension \c dim reused from the mesh during adaption.
     */
    std::array<Index, 4> kept;

    /**
     * \brief Reset statistics to zero
     */
    void reset()
    {
      added.fill(0);
      removed.fill(0);
      kept.fill(0);
    }

    /**
     * \brief Indicate a entity was added
     *
     * \tparam dim_ Dimension of added entities
     * \param[in] count Number of entities added
     */
    template<int dim_>
    void added_element(Index count = 1)
    {
      added[dim_] += count;
    }

    /**
     * \brief Indicate a entity was removed
     *
     * \tparam dim_ Dimension of removed entities
     * \param[in] count Number of entities removed
     */
    template<int dim_>
    void removed_element(Index count = 1)
    {
      removed[dim_] += count;
    }

    /**
     * \brief Indicate entities were removed
     *
     * \param[in] counts Number of entities of each dimension removed
     */
    void removed_elements(const std::array<Index, 4>& counts)
    {
      for(Index i = 0; i < 4; i++)
      {
        removed.at(i) += counts.at(i);
      }
    }

    /**
     * \brief Indicate a entity was kept
     *
     * \tparam dim_ Dimension of kep entities
     * \param[in] count Number of entities kept
     */
    template<int dim_>
    void kept_element(Index count = 1)
    {
      kept[dim_] += count;
    }
  };

  /**
   * \brief Dynamic mesh data structure
   *
   * Mesh data structure that allows adaptively refining a ConformalMesh guided
   * by subdivision levels assigned to the vertices of the ConformalMesh. The
   * refinement is template based. That means mesh entities are replaced with
   * refined versions chosen from a fixed set of possible refinements.
   *
   * See tutorial_09_adaptivemesh for details on how to use this class and the
   * "Adaptive Meshing" wiki entry for documentation of the ideas and concepts
   * used in this implementation.
   *
   * The AdaptiveMesh supports:
   * - quadrilateral and hexahedral meshes,
   * - adaptation to different subdivision levels with minimal work,
   * - exchangable template sets,
   * - a ConformalMesh-like interface for integrating with the rest of FEAT3,
   * - efficient retrieval of mesh entities, allowing assembly directly on the
   *   AdaptiveMesh.
   *
   * \author Markus Muegge
   */
  template<
    typename TemplateSet_,
    typename Shape_,
    int num_coords_ = Shape_::dimension,
    typename Coord_ = Real>
  class AdaptiveMesh
  {
    static_assert(TemplateSet_::template is_shape_compatible<Shape_>());

  public:
    /// Shape of this mesh
    using ShapeType = Shape_;

    /// Shape dimension
    static constexpr int shape_dim = ShapeType::dimension;

    /// World dimension
    static constexpr int world_dim = ShapeType::dimension;

    /// Coordinate type
    using CoordType = Coord_;

    /// Type of underlying mesh
    using FoundationMeshType = ConformalMesh<Shape_, num_coords_, Coord_>;

    /// Vertex type
    using VertexType = typename FoundationMeshType::VertexType;

    /// Template set used by this mesh
    using TemplateSet = TemplateSet_;

    /**
     * \brief ImportBehaviour for this mesh
     *
     * The AdaptiveMesh imports some mesh entities into its own data format.
     * This enum indicate what entities to import.
     */
    enum class ImportBehaviour : std::uint8_t
    {
      RequiredOnly, ///< Import only entites required for the current subdivision levels.
      All,          ///< Import all entities. Useful for exporting the complete mesh later.
    };

  private:
    /// MeshStorage type
    using MeshStorage = Intern::AdaptiveMeshStorage<Shape_, TemplateSet, VertexType>;

    // The following definitions are helpers for the code below

    using VertexKey = typename MeshStorage::template ElementRefByDim<0>;
    using EdgeKey = typename MeshStorage::template ElementRefByDim<1>;
    using FaceKey = typename MeshStorage::template ElementRefByDim<2>;
    using CellKey = typename MeshStorage::template ElementRefByDim<3>;

    using OrientedEdge = typename MeshStorage::template OrientedElementRefByDim<1>;
    using OrientedFace = typename MeshStorage::template OrientedElementRefByDim<2>;

    using AdaptiveVertex = typename MeshStorage::AdaptiveVertexType;

    template<int dim_>
    using LevelTuple =
      Intern::RefinementFieldTuple<
        typename TemplateSet::VertexMarkerType,
        Shape::FaceTraits<typename Shape::FaceTraits<Shape_, dim_>::ShapeType, 0>::count>;

    template<int dim_>
    using ElementRef = typename MeshStorage::template ElementRefByDim<dim_>;

    template<int dim_>
    using OrientedElementRef = typename MeshStorage::template OrientedElementRefByDim<dim_>;

    template<int dim_>
    using ElementTopology = typename MeshStorage::template ElementTopologyByDim<dim_>;

    template<int dim_>
    using ElementChildren = typename MeshStorage::template ElementChildrenByDim<dim_>;

    template<int dim_>
    using AdaptiveElement = typename MeshStorage::template AdaptiveElementByDim<dim_>;

    using MeshRoots = Intern::MeshRoots<shape_dim>;

    using VertexMarkerType = typename TemplateSet::VertexMarkerType;

    // Algorithms that interact with the foundation mesh are written to recurse
    // over the dimension, similar to other algorithms in kernel/geometry. Some
    // of these algorithms require access to internals of the AdaptiveMesh.

    template<typename AdaptiveMeshType_, typename TargetMeshType_, typename AlgShape_>
    friend struct Intern::MeshPartProjector;

    template<typename AdaptiveMeshType_, int topology_dim_, int collection_dim_>
    friend struct Intern::FoundationTopologyCollector;

    /// Reference to foundation mesh
    const FoundationMeshType& _foundation_mesh;

    /// Mesh roots
    Intern::MeshRoots<shape_dim> _roots = {};

    /// Mesh storage
    MeshStorage _storage = {};

    /// Stats of latest adaption
    AdaptionStats _stats = {};

  public:
    /**
     * \brief Constructor
     *
     *  \param[in] mesh Reference to foundation mesh. The reference is stored and must be valid for the lifetime of the
     * AdaptiveMesh \param[in] behaviour The import behaviour to use during adaptions
     */
    explicit AdaptiveMesh(const FoundationMeshType& mesh) : _foundation_mesh(mesh) {};

    AdaptiveMesh(AdaptiveMesh& other) = delete;
    AdaptiveMesh(AdaptiveMesh&& other) = delete;

    AdaptiveMesh& operator=(AdaptiveMesh& other) = delete;
    AdaptiveMesh& operator=(AdaptiveMesh&& other) = delete;

    /**
     * \brief Adapt the mesh to new SubdivisionLevels
     *
     * \param[in] foundation_levels SubdivisionLevels on the foundation mesh
     * \param[in] import_behaviour Enum indicating what entities should be
     * added to layer 0. If RequiredOnly, only elements required for the
     * (partial) refinmement are added, if All, all entities of the foundation
     * mesh are added.
     *
     * \returns Stats about added, removed, and kept elements
     */
    AdaptionStats
    adapt(SubdivisionLevels& foundation_levels, ImportBehaviour import_behaviour = ImportBehaviour::RequiredOnly)
    {
      // We expect every vertex to be marked with a subdivision level
      ASSERTM(foundation_levels.size() == _foundation_mesh.get_num_vertices(), "Invalid SDL size");

      // Algorithm types
      using EntityCollector = Intern::EntityCollector<ShapeType, TemplateSet, FoundationMeshType>;

      // Reset the adaption stats to zero
      _stats.reset();

      // Change foundation_levels to match template set requirements.
      Intern::RefinementField<typename TemplateSet::VertexMarkerType> r_field =
        TemplateSet::make_refinement_field(_foundation_mesh, foundation_levels);

      // We collect sets of all entities that are required to be added to layer
      // 0 of the adaptive mesh. This includes in all cases all elements with
      // non-zero refinement type and all entities making up those elements. Note
      // that sub-entities are added irrespective of their own refinement type,
      // as they are always required for producing a complete refinement tree.
      Intern::MeshIndexSet set;

      // NOTE: Building the sets of all required entities is nice for
      // determining which entities can be safely removed from the mesh, as we
      // can just remove all elements not in the set. Otherwise we would need to
      // determine which entities have both a zero refinement type and are not a
      // sub-entity of an element with non-zero refinement type. We can change
      // this if the hashing becomes a performance bottleneck or the memory-usage
      // causes problems.

      // Actually collect those sets of entities
      bool import_all = import_behaviour == ImportBehaviour::All;
      EntityCollector::collect(set, _foundation_mesh, r_field, import_all);

      // Create, adapt, or replace refinement trees to match new subdivision levels.
      _adapt_roots(set, r_field);

      // Delete all refinement trees of entiites that were previously refined, but are no longer refined.
      _garbage_collect(_roots, set);

      // Re-index the storage to allow for index-based access
      _storage.reindex();

      // Return a copy of the stats for this adapatation.
      return _stats;
    }

    /**
     * \brief Use the AdaptiveMesh as a generator for other meshes
     *
     * \tparam VertexMarker_ Marker for creating subdivision levels.
     * Callable object, which takes an Index as the only argument and returns a std::uint64_t.
     * Used to assign a subdivision level to each vertex.
     *
     * \param[in] foundation_mesh Foundation mesh which gets adaptively refined to produce the new mesh
     * \param[in] marker Vertex marker which determines subdivision levels
     *
     * \returns A new, independent mesh of type FoundationMeshType
     *
     * \note Use this for one off refinements. Create and keep an AdaptiveMesh
     * if you intend to change the mesh's refinement more than once.
     */
    template<typename VertexMarker_>
    static FoundationMeshType create_refined_mesh(const FoundationMeshType& foundation_mesh, const VertexMarker_& marker)
    {
      // Create subdivision levels
      const Index num_vertices = foundation_mesh.get_num_vertices();
      SubdivisionLevels sdls(num_vertices);

      for(Index i(0); i < num_vertices; ++i)
      {
        sdls[i] = marker(i);
      }

      // Create mesh
      AdaptiveMesh adaptive_mesh(foundation_mesh);
      adaptive_mesh.adapt(sdls, ImportBehaviour::All);
      return adaptive_mesh.to_conformal_mesh(Layer{adaptive_mesh.num_layers() - 1});
    }

    /**
     * \brief Create a ConformalMesh from a layer
     *
     * \param[in] layer Index of layer to create a ConformalMesh of
     * \returns Layer \l as a ConformalMesh
     */
    FoundationMeshType to_conformal_mesh(Layer layer) const
    {
      // Algorithm type
      using MeshWriter = Intern::ConformalMeshWriter<AdaptiveMesh, FoundationMeshType, ShapeType>;

      // Determine result mesh size
      std::array<Index, 4> entities = {};
      entities[0] = get_num_entities(layer, 0);
      entities[1] = get_num_entities(layer, 1);
      entities[2] = get_num_entities(layer, 2);
      entities[3] = get_num_entities(layer, 3);
      FoundationMeshType result(entities.data());

      // Write mesh
      MeshWriter::write_to_mesh(result, *this, layer);
      result.fill_neighbors();

      return result;
    }

    /**
     * Convert a mesh part on the regular mesh to a mesh part on the adaptive mesh.
     *
     * All components of the adaptive mesh, that are a (grand)-child of any of
     * the regular components are included. In other words, this method
     * computes the intersection of a regular mesh part with an adaptive mesh,
     * accounting for refinement of elements.
     *
     * Maintains relative order of entities. This means, for example, if the
     * original MeshPart contains edges 1, 4 and 6, then the projected mesh
     * part will first contain all (grand)-children of edge 1, then all
     * (grand-)children of edge 4, and then finally all (grand-)children of
     * edge 6. This, by extension, means that mesh part projection is suitable
     * for mirror assembly.
     *
     * \tparam[in]
     * \param[in] part The MeshPart to project
     * \param[in] layer The layer to project to
     *
     * \returns A mesh part corresponding to the regular mesh part
     */
    template<typename TargetMeshType_>
    MeshPart<TargetMeshType_> project_meshpart(Layer layer, const MeshPart<FoundationMeshType>& part) const
    {
      // Algorithm types
      using Projector = Intern::MeshPartProjector<AdaptiveMesh, TargetMeshType_, ShapeType>;
      using Visitor = Intern::MeshPartProjectorVisitor<MeshStorage>;

      // For projecting a mesh part, we need to collect all child keys on the
      // appropriate layers. This means both the keys of type zero entities
      // that are (grand)-children of any entity in the original mesh part, as
      // well as non-zero type elements on the target layer. The visitor will
      // do so in a PreOrder iteration through the refinement trees,
      // maintaining the relative order of entities.
      Visitor visitor(_storage, layer);

      // Collect keys of mesh part, performs a PreOrder iteration on all
      // refinement trees rooted at an entity of the original mesh part
      Projector::collect_keys(visitor, part, layer, _roots, _storage);

      // Create result mesh part
      std::array<Index, 4> num_elements =
        {visitor.vertices.size(), visitor.edges.size(), visitor.faces.size(), visitor.cells.size()};
      MeshPart<TargetMeshType_> result(num_elements.data(), false);

      // Convert collected keys to indices and write into result mesh part
      Projector::write_meshpart(visitor, result, _storage);

      return result;
    }

    /**
     * \brief Returns the total number of entities of dimension \c dim_ across all mesh layers
     */
    template<int dim_>
    Index num_total_entities() const
    {
      return _storage.template num_total_entities<dim_>();
    }

    /**
     * \brief Retrieve vertex at index \c v from layer \c layer
     */
    VertexType& vertex(Layer layer, Index vertex_index)
    {
      return _storage.get_vertex_by_index(layer, vertex_index).vertex;
    }

    /**
     * \copydoc AdaptiveMesh::vertex()
     */
    const VertexType& vertex(Layer layer, Index vertex_index) const
    {
      return _storage.get_vertex_by_index(layer, vertex_index).vertex;
    }

    /**
     * \brief Returns an index of a sub-entity, i.e. the index of a face of a cell.
     *
     * \tparam dim_ Dimension of the entity
     * \tparam codim_ Dimension of the sub-entity
     * \param[in] layer Layer to retrieve index on
     * \param[in] entity_idx Index of the entity of \c layer
     * \param[in] face_idx Index of the sub-entity
     */
    template<int dim_, int codim_>
    Index get_face_index(Layer layer, Index entity_idx, Index face_idx) const
    {
      auto& entity = _storage.template get_by_index<dim_>(layer, entity_idx);

      return _storage.get_index(entity.topology.template key_by_dim<codim_>(face_idx));
    }

    /**
     * \brief Returns number of elements of dimension \c dim in layer \c layer
     *
     * \param[in] layer
     * The mesh layer to retrieve numbentr of elements from
     * \param[in] dim
     * The dimension of the elements
     */
    Index get_num_entities(Layer layer, int dim) const
    {
      return _storage.num_entities(layer, dim);
    }

    /**
     * \brief Retrieve index of a child element
     *
     * \tparam dim_ Dimension of the child
     * \param[in] layer Layer of parent element
     * \param[in] parent_idx Index of parent
     * \param[in] child Child to retrieve index of. Should fulfill 0 <= child <= get_num_children(layer, parent_idx).
     *
     * \returns The index of the parent on the next layer, if the parent has no
     * children, the requested child's index, if the parent has at least \c
     * child children, std::nullopt else.
     */
    template<int dim_>
    std::optional<Index> get_child(Layer layer, Index parent_idx, Index child)
    {
      auto& element = _storage.template get_by_index<dim_>(layer, parent_idx);

      if(element.type.is_zero_refinement())
      {
        // If the parent has no children, it will be part of all further layers with the same Index.
        // For the purposes of transfer operators and similar, it is thus its own child
        return parent_idx;
      }

      if(child < TemplateSet::template num_children<dim_, dim_>(element.type))
      {
        // A child with the wanted index exists. Retrieve the key and convert to index.
        return _storage.get_index(element.children.template by_dim<dim_>()[child]);
      }

      // No child with the wanted index exists
      return std::nullopt;
    }

    /**
     * \brief Returns the number of children of an element in the mesh
     *
     * \param[in] layer Layer of parent element
     * \param[in] elem Index of parent
     *
     * \returns Number of children of parent
     */
    template<int dim>
    Index get_num_children(Layer layer, Index elem)
    {
      auto type = _storage.template get_by_index<dim>(layer, elem).type;
      if(!type.is_zero_refinement())
      {
        return TemplateSet::template num_children<dim, dim>(type);
      }

      return 1;
    }

    /**
     * \brief Indicates whether any mesh element adjacent to the given vertex has changed on the given layer
     *
     * This is intended to be used with the BPX solver to determine which vertices participate  on the given layer.
     */
    bool has_vertex_changed(Layer layer, Index vertex_idx) const
    {
      auto& vertex = _storage.get_vertex_by_index(layer, vertex_idx);
      return vertex.layer <= layer && layer <= vertex.last_changed;
    }

    /**
     * \brief Computes mapping from elements of the foundation mesh to elements of the AdaptiveMesh
     */
    std::vector<std::pair<Index, Index>> bridge_pairs()
    {
      std::vector<std::pair<Index, Index>> pairs;

      for(auto& iter : _roots.vertices)
      {
        pairs.emplace_back(std::make_pair(iter.first, _storage.get_index(iter.second)));
      }

      return pairs;
    }

    /**
     * \brief Accessor for foundation mesh
     */
    FoundationMeshType& foundation_mesh()
    {
      return _foundation_mesh;
    }

    /**
     * \brief Mapping of elements of the foundation mesh to elements of the AdaptiveMesh
     *
     * \returns Index of element corresponding to given foundation mesh
     * element, if one exists. Index belongs to layer 0.
     */
    std::optional<Index> get_overlap_cell(Index foundation_index)
    {
      std::unordered_map<Index, ElementRef<shape_dim>>& root_map = _roots.template by_dim<ShapeType::dimension>();

      auto iter = root_map.find(foundation_index);
      if(iter != root_map.end())
      {
        return _storage.get_index(iter->second);
      }

      return std::nullopt;
    }

    /**
     * \brief Computes mapping from AdaptiveMesh elements to foundation mesh elements
     *
     * \returns A vector v, where v[i] = j indicates that element i of layer 0
     * of the AdaptiveMesh corresponds to element j of the foundation mesh.
     */
    std::vector<Index> overlap_mapping()
    {
      std::vector<Index> result(_storage.num_entities(0, 0));

      for(auto& entry : _roots.vertices)
      {
        result[_storage.get_index(entry.second)] = entry.first;
      }

      return result;
    }

    /**
     * \brief Returns number of layers of the AdaptiveMesh
     *
     * Note that this contains layer 0, which is a copy of (part of) the
     * foundation mesh. There are num_layers() - 1 new mesh layers.
     */
    Index num_layers() const
    {
      return _storage.num_layers();
    }

    /**
     * \brief Exclusive MeshPart factory
     *
     * Creates a MeshPart (without topology) of all elements of the foundation
     * mesh that are not part of the AdaptiveMesh, i.e. all elements that are
     * not overlapped by the AdaptiveMesh. That is, all elements which are not
     * further refined (if ImportBehaviour::RequiredOnly is used) or no
     * elements (if ImportBehaviour::All is used).
     *
     * The created MeshPart is the inverse of the MeshPart created by overlap_meshpart().
     *
     * \returns A  MeshPart describing all cells of the regular mesh not overlapped by this adaptive mesh
     */
    MeshPart<FoundationMeshType> exclusive_meshpart() const
    {
      // Get element roots
      const std::unordered_map<Index, ElementRef<shape_dim>>& root_map = _roots.template by_dim<ShapeType::dimension>();

      // Determine mesh part size
      const Index total_cells = _foundation_mesh.get_num_entities(shape_dim);
      const Index overlap_cells = root_map.size();
      const Index part_cells = total_cells - overlap_cells;
      std::array<Index, shape_dim + 1> part_size = {};
      part_size[shape_dim] = part_cells;

      // Create mesh part without topology
      MeshPart<FoundationMeshType> result(part_size.data(), false);
      auto& target_set = result.template get_target_set<shape_dim>();

      // Write foundation mesh indices into mesh part
      Index tset_index = 0;
      for(Index i(0); i < total_cells; ++i)
      {
        if(root_map.find(i) == root_map.end())
        {
          // Cell i is not a root of the adaptive mesh. Add it.
          target_set[tset_index] = i;
          tset_index++;
        }
      }

      return result;
    }

    /**
     * \brief Overlap MeshPart factory
     *
     * Creates a MeshPart (without topology) of all elements of the foundation
     * mesh that are also part of the AdaptiveMesh, i.e. all elements that are
     * overlapped by the AdaptiveMesh. That is, all elements which are further
     * refined (if ImportBehaviour::RequiredOnly is used) or all elements (if
     * ImportBehaviour::All is used).
     *
     * The created MeshPart is the inverse of the MeshPart created by exclusive_meshpart().
     *
     * \returns A  MeshPart describing all cells of the regular mesh overlapped by this adaptive mesh
     */
    MeshPart<FoundationMeshType> overlap_meshpart() const
    {
      // Get element roots
      const std::unordered_map<Index, ElementRef<shape_dim>>& root_map = _roots.template by_dim<ShapeType::dimension>();

      // Determine result target set size
      const Index overlap_cells = root_map.size();
      std::array<Index, shape_dim + 1> part_size = {};
      part_size[shape_dim] = overlap_cells;

      // Create mesh part without topology
      MeshPart<FoundationMeshType> result(part_size.data(), false);
      auto& target_set = result.template get_target_set<shape_dim>();

      // Write foundation mesh indices into mesh part
      Index tset_index = 0;
      for(auto& entry : root_map)
      {
        target_set[tset_index] = entry.first;
        tset_index++;
      }

      return result;
    }

    /**
     * \brief Transfer subdivision levels from the adaptive mesh to the foundation mesh
     *
     * For each vertex of the finest layer of the AdaptiveMesh the closest
     * foundation mesh vertex is found and its subdivision level is set, such
     * that each founation vertex is assigned the maximum subdivision level of
     * any of its closests vertices.
     *
     * \param[in] lvl_in SubdivisionLevels on the AdaptiveMesh. Assigns a
     * subdivision level to each vertex of the finest adaptive mesh level
     *
     * \returns SubdivisionLevels on the foundation mesh. Computed such that subdivision levels assigned to
     */
    SubdivisionLevels transfer_sdls(SubdivisionLevels& lvls_in) const
    {
      const Index num_layers = _storage.num_layers();
      const Layer fine_layer = Layer{num_layers - 1};
      const Index num_base_vertices = _storage.num_entities(Layer{0}, 0);
      const Index num_fine_vertices = _storage.num_entities(fine_layer, 0);

      SubdivisionLevels result(num_base_vertices);

      for(Index v_fine_idx = 0; v_fine_idx < num_fine_vertices; v_fine_idx++)
      {
        // Find vertex in foundation mesh that is closest to the current fine vertex
        Real closest_distance = std::numeric_limits<Real>::max();
        Index closest_vertex = 0;
        for(Index v_base_idx = 0; v_base_idx < num_base_vertices; v_base_idx++)
        {
          const AdaptiveVertex& v_fine = _storage.get_vertex_by_index(fine_layer, v_fine_idx);
          const AdaptiveVertex& v_base = _storage.get_vertex_by_index(Layer{0}, v_base_idx);
          const Real distance = (v_fine.vertex - v_base.vertex).norm_euclid_sqr();
          if(distance < closest_distance)
          {
            closest_distance = distance;
            closest_vertex = v_base_idx;
          }
        }

        // Assign the subdivision levels of the current vertex to the closest vertex in the foundation mesh
        result[closest_vertex] = std::max(result[closest_vertex], lvls_in[v_fine_idx]);
      }
      return result;
    }

    /**
     * \brief Computes neighboring elements for all mesh elements
     */
    void fill_neighbors()
    {
      _storage.fill_neighbors();
    }

    /**
     * \brief Returns the index of the n-th neighbor of the given element
     *
     * \returns The index of the neighbor, if one exists, ~Index(0) otherwise.
     */
    Index get_neighbor(Layer layer, Index element_idx, Index neighbor_idx) const
    {
      return _storage.get_neighbor(layer, element_idx, neighbor_idx);
    }

    /**
     * \brief Applies a "proper rigid" transformation onto the mesh vertices.
     *
     * Let \e v denote the \p origin world point, \e w the \p offset world point and \e R
     * the rotation matrix corresponding to the \p angles, then this function applies the
     * following transformation for any vertex \e x of the vertex set:
     *
     *   \f[ x \mapsto w + R\cdot (x - v) \f]
     *
     * \param[in] origin
     * The origin of the transformation. This is subtracted from any vertex before applying the
     * rotation.
     *
     * \param[in] angles
     * The angles of the rotation matrix.
     * - 2D: the rotation angle in radians is stored as:
     *   - angles(0): rotation angle
     *   - angles(1): \e ignored
     * - 3D: the rotation angles in radians stored as:
     *   - angles(0): yaw angle
     *   - angles(1): pitch angle
     *   - angles(2): roll angle
     *
     * \param[in] offset
     * The offset of the transformation. This is added to any vertex after applying the rotation.
     *
     * \warning Affects all layers of the mesh
     */
    void transform(const VertexType& origin, const VertexType& angles, const VertexType& offset)
    {
      _storage.transform(origin, angles, offset);
    }

  private:
    /**
     * \brief Construct a new subtree
     *
     * The root element will be refined using templates from the TemplateSet template parameter.
     *
     * \param[in] depth
     * The starting depth for the new subtree
     * \param[in] levels
     * The refinement levels for the element, dictates further refinement
     * \param[in] topology
     * The surrounding elements of the new element
     *
     * \returns The key of the root element of the new subtree
     */
    template<int dim_>
    Intern::ElementKey<dim_> _build(Index depth, const LevelTuple<dim_>& levels, const ElementTopology<dim_>& topology)
    {
      ASSERT(_dbg_check_topology_unique<dim_>(topology));
      ASSERT(_dbg_check_topology_orientations<dim_>(topology));
      ASSERT(_dbg_is_topology_consistent<dim_>(topology));
      ASSERT(_dbg_is_topology_layering_consistent<dim_>(Intern::Layer{depth}, topology));

      _stats.added_element<dim_>();

      auto type = typename TemplateSet::template RefinementTypeByDim<dim_>(levels);

      //std::cout << "Building element of type " << type << "\n";
      // Create tree root
      auto self = AdaptiveElement<dim_>(type, Layer{depth}, topology);

      if(type.is_zero_refinement())
      {
        // Zero type elements have no children. We are done here.
        auto key = _storage.insert(self);
        //std::cout << "Built element " << key << "\n";
        return key;
      }

      using TemplateShape = typename Shape::FaceTraits<ShapeType, dim_>::ShapeType;

      static constexpr int num_vertices = Shape::FaceTraits<TemplateShape, 0>::count;

      const RefinementTemplate<TemplateShape>& tmplt = TemplateSet::get_template(type);

      // Build child vertices
      _stats.added_element<0>(tmplt.template num_entities<0>());

      // Retrieve vertex coordinates of topology for following vertex interpolation
      std::array<VertexType, num_vertices> vertex_coordinates;

      for(int i{0}; i < num_vertices; i++)
      {
        vertex_coordinates[i] = _storage[topology.template key_by_dim<0>(i)].vertex;
      }

      // Produce child vertices by linear interpolation
      auto& child_vertices = self.children.template by_dim<0>();
      for(Index i = 0; i < tmplt.template num_entities<0>(); i++)
      {
        const VertexType vertex(Intern::interpolate(vertex_coordinates, tmplt.get_vertex_coefficients()[i]));
        child_vertices[i] = _storage.insert(AdaptiveVertex(vertex, Layer{depth + 1}));
      }

      // We are recording the last mesh layer on which any of the entities surrounding a vertex has changed.
      // This allows disregarding these vertices in a BPX solver later.
      // NOTE: This is quite specific. Can we move this out of the build routine somehow?
      if(tmplt.template num_entities<0>() > 0)
      {
        // Parent vertices created new vertices at this depth.
        // Increase referenced depth for BPX star-sum
        auto& topo_vertices = self.topology.template by_dim<0>();
        for(Index v = 0; v < topo_vertices.size(); v++)
        {
          auto& vertex = _storage[topo_vertices[v]];
          if(depth + 1 > vertex.last_changed.idx)
          {
            vertex.last_changed = Layer{depth + 1};
          }
        }
      }

      // Build child edges
      if constexpr(dim_ >= 1)
      {
        auto& children = self.children.template by_dim<1>();
        const auto& topo_templates = tmplt.template get_topologies<1>();
        for(Index i = 0; i < tmplt.template num_entities<1>(); i++)
        {
          ElementTopology<1> child_topo;
          _build_topology<dim_, 1>(topo_templates[i], child_topo, topology, self.children);

          children[i] = _build<1>(depth + 1, _levels<dim_, 1>(topo_templates[i], levels), child_topo);
        }
      }

      // Build faces
      if constexpr(dim_ >= 2)
      {
        auto& children = self.children.template by_dim<2>();
        const auto& topo_templates = tmplt.template get_topologies<2>();
        for(Index i = 0; i < tmplt.template num_entities<2>(); i++)
        {
          ElementTopology<2> child_topo;
          _build_topology<dim_, 2>(topo_templates[i], child_topo, topology, self.children);

          children[i] = _build<2>(depth + 1, _levels<dim_, 2>(topo_templates[i], levels), child_topo);
        }
      }

      // Build cells
      if constexpr(dim_ >= 3)
      {
        auto& children = self.children.template by_dim<3>();
        const auto& topo_templates = tmplt.template get_topologies<3>();
        for(Index i = 0; i < tmplt.template num_entities<3>(); i++)
        {
          ElementTopology<3> child_topo;
          _build_topology<dim_, 3>(topo_templates[i], child_topo, topology, self.children);

          children[i] = _build<3>(depth + 1, _levels<dim_, 3>(topo_templates[i], levels), child_topo);
        }
      }

      return _storage.insert(self);
    }

    /**
     * \brief Adapts roots of refinement trees
     *
     * \tparam dim_ Dimension of entities
     * \param[in] set Set of indices that should be part of the mesh after adaptation
     * \param[in] sdls SubdivisionLevels to guide adaptation
     *
     * This method iterates over all entities of dimension dim_ in the MeshIndexSet and either
     * * adapt its children, if the entity already exists in the adaptive mesh and its type is unchanged
     * * replaces it, if the entity already exists in the adaptive mesh and its type is changed
     * * builds it, if the entity does not yet exist
     *
     * This method does not erase entities that are part of the AdaptiveMesh but not in the MeshIndexSet.
     */
    template<int dim_>
    void _adapt_roots_of_dim(const Intern::MeshIndexSet& set, Intern::RefinementField<typename TemplateSet::VertexMarkerType>& sdls)
    {
      // Algorithm for collecting topologies from the foundation mesh. Entry
      // dimension is dim_ - 1, because topologies only contain entities of
      // smaller dimension than the entity the topology belongs to.
      using TopologyCollector = Intern::FoundationTopologyCollector<AdaptiveMesh, dim_, dim_ - 1>;

      // Vertices of foundation entities, for collecting subdivision levels
      auto& foundation_entities = _foundation_mesh.template get_index_set<dim_, 0>();
      auto& roots = _roots.template by_dim<dim_>();

      for(Index idx : set.template by_dim<dim_>())
      {
        // Determine refinement type of current entity
        auto markings = sdls.get_tuple(foundation_entities[idx]);
        auto type = typename TemplateSet::template RefinementTypeByDim<dim_>(markings);

        // Determine entities entry in mesh roots, if it exists
        auto iter = roots.find(idx);
        bool exists = iter != roots.end();

        if(exists && _storage.type(iter->second) == type)
        {
          // Element already exists and its type has not changed.
          // Adapt its children further.
          _adapt_children<dim_>(iter->second, markings);
        }
        else
        {
          //std::cout << "Rebuilding element with new type " << type << "\n";
          //std::cout << "Markings: ";
          for(Index i(0); i < markings.size; i++)
          {
            //std::cout << markings[i] << ", ";
          }
          //std::cout << "\n";
          // Element either does not exist or its type has changed.
          // (Re)build it.

          // We need the entity's topology to rebuild it.
          // NOTE: Even if the entity already exists in the mesh and is only
          // going to be replaced, we can not reuse its topology. If we entered
          // this branch the entity's refinement type has changed, which means
          // the refinement type of at least one sub-entity has also changed.
          // That sub-entity has itself been replaced at this point, and the
          // reference stored in this entity's topology is invalid.

          // OPTIMIZATION: While it is true that we can't reuse the previous
          // topology, if it exists, we are currently doing some needless work
          // in that case. The topology collection also determines the
          // orientation of any references in the topology. Those orientations
          // don't change, so it would be enough to only update the references,
          // if the topology has been collected before.
          ElementTopology<dim_> topology;
          TopologyCollector::collect(topology, idx, *this);

          if(exists)
          {
            // Entity already existed. Replace it.
            _erase(iter->second);
            iter->second = _build<dim_>(0, markings, topology);
          }
          else
          {
            // Entity is entirely new. Build it.
            roots[idx] = _build<dim_>(0, markings, topology);
          }
        }
      }
    }

    /**
     * \brief Adapts roots of refinement trees
     *
     * \param[in] set Set of indices that should be part of the mesh after adaptation
     * \param[in] sdls SubdivisionLevels to guide adaptation
     *
     * This method iterates over all entities in the MeshIndexSet and either
     * * adapt its children, if the entity already exists in the adaptive mesh and its type is unchanged
     * * replaces it, if the entity already exists in the adaptive mesh and its type is changed
     * * builds it, if the entity does not yet exist
     *
     * This method does not erase entities that are part of the AdaptiveMesh but not in the MeshIndexSet.
     */
    void _adapt_roots(const Intern::MeshIndexSet& set, Intern::RefinementField<typename TemplateSet::VertexMarkerType>& sdls)
    {
      if constexpr(shape_dim >= 0)
      {
        auto& vertex_set = _foundation_mesh.get_vertex_set();
        auto& vertex_roots = _roots.template by_dim<0>();
        for(const auto& vert : set.vertices)
        {
          if(vertex_roots.find(vert) == vertex_roots.end())
          {
            _stats.template added_element<0>();
            vertex_roots[vert] = _storage.insert(vertex_set[vert], Layer{0});
          }
          else
          {
            _stats.template kept_element<0>();
          }
        }
      }
      if constexpr(shape_dim >= 1)
      {
        _adapt_roots_of_dim<1>(set, sdls);
      }
      if constexpr(shape_dim >= 2)
      {
        _adapt_roots_of_dim<2>(set, sdls);
      }
      if constexpr(shape_dim >= 3)
      {
        _adapt_roots_of_dim<3>(set, sdls);
      }
    }

    /**
     * \brief Adjusts the children of the given element to match new subdivision levels
     *
     * Children will be
     * * further adapted, if their types are unchanged
     * * replaced, if their type has changed.
     *
     * \param[in] element
     * The element whose children will be adapted
     * \param[in] levels
     * New subdivision levels of the element
     */
    template<int dim_>
    void _adapt_children(ElementRef<dim_> key, const LevelTuple<dim_>& levels)
    {
      // Calculate element type
      const auto type = _storage.type(key);
      const Index depth = _storage[key].layer.idx;

      //std::cout << "Adapting element of type " << type << "\n";
      // Adapting an entity's children requires a valid topology for the parent
      // element.  If this entity's type changed, then it should have been
      // rebuilt with a valid topology and calling _adapt_children is an error.
      XASSERTM(type == typename TemplateSet::template RefinementTypeByDim<dim_>(levels), "Called adapt_children on element with changed type.");

      // We are keeping the current entity as part of the mesh
      _stats.kept_element<dim_>();

      if(type.is_zero_refinement())
      {
        // Zero-type elements have no children that need further adaptation.
        return;
      }

      // Retrieve reference to current entity
      auto& self = _storage[key];

      using TemplateShape = typename Shape::FaceTraits<ShapeType, dim_>::ShapeType;
      const RefinementTemplate<TemplateShape>& tmplt = TemplateSet::get_template(type);

      // Retrieve vertex query from template set. Vertices do not need
      // adaptation, but we guarantee template sets that we will construct all
      // template queries in order.
      _stats.kept_element<0>(tmplt.template num_entities<0>());

      // Adapt child edges
      if constexpr(dim_ >= 1)
      {
        _adapt_sub_entities<dim_, 1>(depth, self, tmplt, levels);

        // Adapt child faces
        if constexpr(dim_ >= 2)
        {
          _adapt_sub_entities<dim_, 2>(depth, self, tmplt, levels);

          // Adapt child cells
          if constexpr(dim_ >= 3)
          {
            _adapt_sub_entities<dim_, 3>(depth, self, tmplt, levels);
          }
        }
      }
    }

  protected:
    /**
     * \brief Erases an element and all its children from the mesh
     *
     * \param[in] key Key pointin to root of subtree to erase
     *
     * Wrapper around MeshStorage::erase for stat-keeping purposes
     */
    template<int dim_>
    void _erase(Intern::ElementKey<dim_> key)
    {
      _stats.removed_elements(_storage.erase(key));
    }

    /**
     * \brief Adapts children of dimension dim_
     *
     * \tparam parent_dim_ Dimension of parent
     * \tparam child_dim_ Dimension of children
     * \tparam TQuery_ Template query type
     * \tparam MQuery_ Mesh query type
     *
     * \param[in] depth Current refinement tree depth
     * \param[in] self Current AdaptiveElement
     * \param[in] t_query Template query for current dimension
     * \param[in] m_query Mesh Query for self
     */
    template<int parent_dim_, int child_dim_>
    void _adapt_sub_entities(
        Index depth,
        AdaptiveElement<parent_dim_>& self,
        const RefinementTemplate<typename Shape::FaceTraits<ShapeType, parent_dim_>::ShapeType>& tmplt,
        const LevelTuple<parent_dim_>& levels)
    {
      // children is an std::array of slotmap keys
      auto& children = self.children.template by_dim<child_dim_>();
      auto& topo_templates = tmplt.template get_topologies<child_dim_>();

      for(Index i = 0; i < tmplt.template num_entities<child_dim_>(); i++)
      {
        // Determine child SubdivisionLevels and type
        LevelTuple<child_dim_> child_levels = _levels<parent_dim_, child_dim_>(topo_templates[i], levels);
        auto child_type = typename TemplateSet::template RefinementTypeByDim<child_dim_>(child_levels);

        if(_storage.type(children[i]) == child_type)
        {
          // Type is unchanged. Adapt children recursively.
          _adapt_children<child_dim_>(children[i], child_levels);
        }
        else
        {
          // Type is changed. Erase child and rebuild it.
          _erase(children[i]);

          ElementTopology<child_dim_> child_topology;
          _build_topology<parent_dim_, child_dim_>(topo_templates[i], child_topology, self.topology, self.children);

          children[i] = _build<child_dim_>(depth + 1, child_levels, child_topology);
        }
      }
    }

    /**
     * \brief Spread refinement levels to a child entity
     *
     * \tparam parent_dim_ Dimension of the parent
     * \tparam entity_dim_ Dimension of the child entity
     *
     * \param tmplt Topology template for the children
     * \param lvls Refinement field markings of the parent
     *
     * Calculates the vertex markings for a child element, based on its
     * topology template.
     */
    template<int parent_dim_, int entity_dim_>
    LevelTuple<entity_dim_> _levels(
        const TopologyTemplate<typename Shape::FaceTraits<ShapeType, entity_dim_>::ShapeType>& tmplt,
        const LevelTuple<parent_dim_>& lvls)
    {
      // OPTIMIZATION(mmuegge): Currenlty a new minimum is calculated for all
      // vertices. We can use that all sibling vertices share the same new
      // marking. As do all vertices on any of the boundary elements. Add a
      // MarkingSpreader class or something like that and cache the results.
      using ParentShape = typename Shape::FaceTraits<ShapeType, parent_dim_>::ShapeType;

      const auto& vertex_refs = tmplt.template get_references<0>();

      LevelTuple<entity_dim_> result;

      for(Index i(0); i < result.size; ++i)
      {
        result[i] = TemplateSet::template spread_refinement_field<ParentShape>(vertex_refs[i], lvls);
      }

      return result;
    }

    /**
     * \brief Construct a topology for a new child element
     *
     * \tparam parent_dim_ Dimension of the parent element
     * \tparam topo_dim Dimension of the topology
     * \tparam dim_ Current working dim
     *
     * \param[in] tmplt Topology template
     * \param[inout] topo Topology being constructed
     * \param[in] parent_topo Parent topology
     * \param[in] siblings Sibling elemets, aka children of the parent entity
     *
     * Implemented via recursion on dim_.
     */
    template<int parent_dim_, int topo_dim_, int dim_ = topo_dim_ - 1>
    void _build_topology(
        const TopologyTemplate<typename Shape::FaceTraits<ShapeType, topo_dim_>::ShapeType>& tmplt,
        ElementTopology<topo_dim_>& topo,
        const ElementTopology<parent_dim_>& parent_topo,
        const ElementChildren<parent_dim_>& siblings)
    {
      if constexpr(dim_ >= 0)
      {
        _build_topology<parent_dim_, topo_dim_, dim_ - 1>(tmplt, topo, parent_topo, siblings);

        using TopologyShape = typename Shape::FaceTraits<ShapeType, topo_dim_>::ShapeType;
        static constexpr int num_entities = Shape::FaceTraits<TopologyShape, dim_>::count;

        auto& target = topo.template by_dim<dim_>();
        for(int i(0); i < num_entities; i++)
        {
          EntityReference ref = tmplt.template get_reference<dim_>(i);
          target[i] = _resolve_entity_reference<parent_dim_, dim_>(parent_topo, siblings, ref);

          if constexpr(dim_ >= 1)
          {
            if(ref.source == EntitySource::BoundaryEdge || ref.source == EntitySource::BoundaryFace)
            {
              // Calc orientation
              using ResultShape = typename Shape::FaceTraits<ShapeType, dim_>::ShapeType;
              using FaceMapping = Intern::FaceIndexMapping<TopologyShape, dim_, 0>;
              using CongruencySampler = Intern::CongruencySampler<ResultShape>;

              static constexpr int num_vertices = Shape::FaceTraits<ResultShape, 0>::count;
              std::array<ElementRef<0>, num_vertices> local_keys;
              std::array<ElementRef<0>, num_vertices> foreign_keys;

              const auto& entity = _storage[target[i]];
              for(int vertex(0); vertex < num_vertices; vertex++)
              {
                Index local_vertex = FaceMapping::map(i, vertex);
                local_keys[vertex] = topo.template key_by_dim<0>(local_vertex).key;
                foreign_keys[vertex] = entity.topology.template key_by_dim<0>(vertex).key;
              }

              target[i].orientation = CongruencySampler::compare(local_keys.data(), foreign_keys.data());
            }
          }
        }
      }
    }

    /*
     * \brief Map a EntityReference to an OrientedElementRef
     *
     * \tparam entity_dim_ Dimension of the current entity
     * \tparam result_dim_ Dimension the reference refers to
     *
     * \param[in] topo Current topology
     * \param[in] siblings Current siblings
     * \param[in] ref EntityReference
     */
    template<int entity_dim_, int result_dim_>
    OrientedElementRef<result_dim_> _resolve_entity_reference(
        const ElementTopology<entity_dim_>& topo,
        const ElementChildren<entity_dim_>& siblings,
        const EntityReference& ref)
    {
      switch(ref.source)
      {
        case EntitySource::ParentTopology:
          return topo.template key_by_dim<result_dim_>(ref.index);
        case EntitySource::Sibling:
          return OrientedElementRef<result_dim_>(ref.orientation, siblings.template by_dim<result_dim_>()[ref.index]);
        case EntitySource::BoundaryEdge:
        {
          if constexpr (entity_dim_ >= 2 && result_dim_ <= 1)
          {
            return _resolve_boundary_entity_reference<entity_dim_, result_dim_, 1>(topo, ref);
          }
        }
        case EntitySource::BoundaryFace:
        {
          if constexpr (entity_dim_ >= 3 && result_dim_ <= 2)
          {
            return _resolve_boundary_entity_reference<entity_dim_, result_dim_, 2>(topo, ref);
          }
        }
        default:
          XABORTM("select failed");
      }
    }

    /// Helper-function for _resolve_entity_reference. Handles boundary entities of any dimension.
    template<int entity_dim_, int result_dim_, int boundary_dim_>
    OrientedElementRef<result_dim_> _resolve_boundary_entity_reference(
        const ElementTopology<entity_dim_>& topo,
        const EntityReference& ref)
    {
      ASSERT(ref.source == EntitySource::BoundaryEdge || ref.source == EntitySource::BoundaryFace);

      // Retrieve type of boundary element
      using BoundaryShape = typename Shape::FaceTraits<ShapeType, boundary_dim_>::ShapeType;
      OrientedElementRef<boundary_dim_> boundary_ref = topo.template key_by_dim<boundary_dim_>(ref.entity);
      typename TemplateSet::template RefinementTypeByDim<boundary_dim_> type = _storage.type(boundary_ref.key);

      // Use type and orientation to correct the index of the reference for orientation
      std::pair<Index, int> orientation_mapping =
        TemplateSet::template correct_for_orientation<BoundaryShape, result_dim_>(type, boundary_ref.orientation, ref.index);

      // Retrieve correct child based on orientation mapping
      AdaptiveElement<boundary_dim_>& boundary_element = _storage[boundary_ref];
      ElementRef<result_dim_> result_ref =
        boundary_element.children.template by_dim<result_dim_>()[orientation_mapping.first];

      // Return oriented reference
      return OrientedElementRef<result_dim_>(orientation_mapping.second, result_ref);
    }

    /**
     * \brief Helper for _garbage_collect. Deletes all entities of dimension dim_ that are not in the given set.
     *
     * \param[in, out] roots
     * MeshRoots to operate on
     * \param[in] set
     * Set of mesh indices to keep
     */
    template<int dim_>
    void _delete_unused(std::unordered_map<Index, ElementRef<dim_>>& roots, const std::set<Index>& elements)
    {
      auto iter = roots.begin();
      while(iter != roots.end())
      {
        if(elements.find(iter->first) == elements.end())
        {
          _erase(iter->second);
          iter = roots.erase(iter);
        }
        else
        {
          iter++;
        }
      }
    }

    /**
     * \brief Delete all elements from the storage and mesh roots that are not in the index set.
     *
     * \param[in, out] roots
     * MeshRoots to operate on
     * \param[in] set
     * Set of mesh indices to keep
     */
    template<int dim_ = shape_dim>
    void _garbage_collect(Intern::MeshRoots<shape_dim>& roots, const Intern::MeshIndexSet& set)
    {
      if constexpr(dim_ >= 0)
      {
        _garbage_collect<dim_ -1>(roots, set);
        _delete_unused<dim_>(roots.template by_dim<dim_>(), set.by_dim<dim_>());
      }
    }

    template<int topo_dim_, int dim_ = topo_dim_ - 1>
    bool _dbg_check_topology_unique(const ElementTopology<topo_dim_>& topology)
    {
      if constexpr(dim_ >= 0)
      {
        using TopologyShape = typename Shape::FaceTraits<ShapeType, topo_dim_>::ShapeType;

        bool all_unique = true;

        static constexpr int num_entities = Shape::FaceTraits<TopologyShape, dim_>::count;

        const auto& entities = topology.template by_dim<dim_>();
        for(int i(0); i < num_entities; i++)
        {
          for(int j(0); j < num_entities; j++)
          {
            all_unique = all_unique && (i == j || (entities[i] != entities[j]));
          }
        }

        return all_unique && _dbg_check_topology_unique<topo_dim_, dim_ - 1>(topology);
      }
      else
      {
        return true;
      }
    }

    template<int topo_dim_, int dim_ = topo_dim_ - 1>
    bool _dbg_check_topology_orientations(const ElementTopology<topo_dim_>& topology)
    {
      if constexpr(dim_ == 0)
      {
        return true;
      }
      else
      {
        bool is_oriented = true;

        using TopologyShape = typename Shape::FaceTraits<ShapeType, topo_dim_>::ShapeType;
        using EntityShape = typename Shape::FaceTraits<ShapeType, dim_>::ShapeType;

        using FaceMapping = Intern::FaceIndexMapping<TopologyShape, dim_, 0>;
        using CongruencySampler = Intern::CongruencySampler<EntityShape>;

        static constexpr int num_entities = Shape::FaceTraits<TopologyShape, dim_>::count;
        static constexpr int num_vertices = Shape::FaceTraits<EntityShape, 0>::count;

        const auto& entities = topology.template by_dim<dim_>();

        for(int entity_idx(0); entity_idx < num_entities; entity_idx++)
        {
          const Intern::OrientedElement<dim_> o_ref = entities[entity_idx];
          const auto& entity = _storage[o_ref];

          std::array<ElementRef<0>, num_vertices> local_keys;
          std::array<ElementRef<0>, num_vertices> foreign_keys;

          for(int vertex(0); vertex < num_vertices; vertex++)
          {
            Index local_vertex = FaceMapping::map(entity_idx, vertex);
            local_keys[vertex] = topology.template key_by_dim<0>(local_vertex).key;
            foreign_keys[vertex] = entity.topology.template key_by_dim<0>(vertex).key;
          }

          is_oriented = is_oriented && (CongruencySampler::compare(local_keys.data(), foreign_keys.data()) == o_ref.orientation);
        }

        return is_oriented && _dbg_check_topology_orientations<topo_dim_, dim_ - 1>(topology);
      }
    }

    template<int topo_dim_, int dim_ = topo_dim_ - 1>
    bool _dbg_is_topology_consistent(const ElementTopology<topo_dim_>& topology)
    {
      if constexpr(dim_ == 0)
      {
        return true;
      }
      else
      {
        bool is_consistent = true;

        using TopologyShape = typename Shape::FaceTraits<ShapeType, topo_dim_>::ShapeType;
        using EntityShape = typename Shape::FaceTraits<ShapeType, dim_>::ShapeType;

        using FaceMapping = Intern::FaceIndexMapping<TopologyShape, dim_, 0>;
        using CongruencyMapping = Intern::CongruencyMapping<EntityShape, 0>;

        static constexpr int num_entities = Shape::FaceTraits<TopologyShape, dim_>::count;
        static constexpr int num_vertices = Shape::FaceTraits<EntityShape, 0>::count;

        const auto& entities = topology.template by_dim<dim_>();
        for(int entity_idx(0); entity_idx < num_entities; entity_idx++)
        {
          const auto& entity = _storage[entities[entity_idx]];

          for(int vertex(0); vertex < num_vertices; vertex++)
          {
            ElementRef<0> v_topo = topology.template by_dim<0>()[FaceMapping::map(entity_idx, vertex)].key;
            int orientation = entities[entity_idx].orientation;
            Index mapped = CongruencyMapping::map(orientation, vertex);
            ElementRef<0> v_entity = entity.topology.template by_dim<0>()[mapped].key;
            is_consistent = is_consistent && (v_topo == v_entity);

            if(v_topo != v_entity)
            {
              std::cout << "Inconsistent topology. Vertex " << vertex << " of entity " << entity_idx << " of dimension << " << dim_ << " mismatches!\n";
            }
          }
        }

        return is_consistent && _dbg_is_topology_consistent<topo_dim_, dim_ - 1>(topology);
      }
    }

    template<int topo_dim_, int dim_ = topo_dim_ - 1>
    bool _dbg_is_topology_layering_consistent(Intern::Layer layer, const ElementTopology<topo_dim_>& topology)
    {
      if constexpr(dim_ == 0)
      {
        bool is_consistent = true;

        // Vertices are layered consistently as long as we pick no vertices from a higher layer
        // NOTE(mmuegge): This can be merged with the dim_ != 0 case, if we
        // mark vertices as permanent by default. Which we should do anyway,
        // because they are
        for(const auto& entity : topology.template by_dim<0>())
        {
          is_consistent = is_consistent && entity.key.layer <= layer;
        }

        return is_consistent;
      }
      else
      {
        bool is_consistent = _dbg_is_topology_layering_consistent<topo_dim_, dim_ - 1>(layer, topology);

        for(const auto& entity : topology.template by_dim<dim_>())
        {
          is_consistent = is_consistent && (entity.key.layer == layer || (entity.key.layer < layer && entity.key.is_permanent));
        }

        return is_consistent;
      }
    }
  }; // class AdaptiveMesh
} // namespace FEAT::Geometry

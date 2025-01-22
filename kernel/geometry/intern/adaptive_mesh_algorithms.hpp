// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/geometry/intern/refinement_field.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/base_header.hpp>
#include <kernel/geometry/intern/adaptive_mesh_storage.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/shape.hpp>

namespace FEAT::Geometry::Intern
{
#ifdef DOXYGEN
  /**
   * \brief Algorithm for exporting a layer of an AdaptiveMesh to a ConformalMesh
   *
   * This is the main declaration of the algorithm. The implementation is split
   * across specializations for individual dimension to allow the algorithm to
   * interact with meshes of any dimension.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename ConformalMesh_, typename Shape_>
  struct ConformalMeshWriter
  {
    /**
     * \brief Write entites of the current shape into the target ConformalMesh
     *
     * \param[out] target The ConformalMesh to write into
     * \param[in] source The AdaptiveMesh to write from
     * \param[in] layer The layer to export
     */
    static void write_to_mesh(ConformalMesh_& target, AdaptiveMesh_& source, Index layer)
    {
    }
  };
#else
  template<typename AdaptiveMesh_, typename ConformalMesh_, typename Shape_>
  struct ConformalMeshWriter;
#endif

  /**
   * \brief Algorithm for exporting a layer of an AdaptiveMesh to a ConformalMesh
   *
   * Specialization for vertices
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename ConformalMesh_>
  struct ConformalMeshWriter<AdaptiveMesh_, ConformalMesh_, Shape::Vertex>
  {
    /// \copdydoc ConformalMeshWriter::write_to_mesh
    static void write_to_mesh(ConformalMesh_& target, const AdaptiveMesh_& source, Layer layer)
    {
      auto& vertices = target.get_vertex_set();

      Index num_vertices = source.get_num_entities(layer, 0);
      for(Index i = 0; i < num_vertices; i++)
      {
        vertices[i] = source.vertex(layer, i);
      }
    }
  };

  /**
   * \brief Algorithm for exporting a layer of an AdaptiveMesh to a ConformalMesh
   *
   * Specialization for edges
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename ConformalMesh_>
  struct ConformalMeshWriter<AdaptiveMesh_, ConformalMesh_, Shape::Hypercube<1>>
  {
    /// \copdydoc ConformalMeshWriter::write_to_mesh
    static void write_to_mesh(ConformalMesh_& target, const AdaptiveMesh_& source, Layer layer)
    {
      // Write vertices
      ConformalMeshWriter<AdaptiveMesh_, ConformalMesh_, Shape::Vertex>::write_to_mesh(target, source, layer);

      auto& v_at_e = target.template get_index_set<1, 0>();

      Index num_edges = source.get_num_entities(layer, 1);
      for(Index i = 0; i < num_edges; i++)
      {
        v_at_e[i][0] = source.template get_face_index<1, 0>(layer, i, 0);
        v_at_e[i][1] = source.template get_face_index<1, 0>(layer, i, 1);
      }
    }
  };

  /**
   * \brief Algorithm for exporting a layer of an AdaptiveMesh to a ConformalMesh
   *
   * Specialization for faces
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename ConformalMesh_>
  struct ConformalMeshWriter<AdaptiveMesh_, ConformalMesh_, Shape::Hypercube<2>>
  {
    /// \copdydoc ConformalMeshWriter::write_to_mesh
    static void write_to_mesh(ConformalMesh_& target, const AdaptiveMesh_& source, Layer layer)
    {
      // Write edges
      ConformalMeshWriter<AdaptiveMesh_, ConformalMesh_, Shape::Hypercube<1>>::write_to_mesh(target, source, layer);

      auto& v_at_f = target.template get_index_set<2, 0>();
      auto& e_at_f = target.template get_index_set<2, 1>();

      Index num_faces = source.get_num_entities(layer, 2);
      for(Index i = 0; i < num_faces; i++)
      {
        for(int j = 0; j < v_at_f.get_num_indices(); j++)
        {
          v_at_f[i][j] = source.template get_face_index<2, 0>(layer, i, j);
        }

        for(int j = 0; j < e_at_f.get_num_indices(); j++)
        {
          e_at_f[i][j] = source.template get_face_index<2, 1>(layer, i, j);
        }
      }
    }
  };

  /**
   * \brief Algorithm for exporting a layer of an AdaptiveMesh to a ConformalMesh
   *
   * Specialization for cells
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename ConformalMesh_>
  struct ConformalMeshWriter<AdaptiveMesh_, ConformalMesh_, Shape::Hypercube<3>>
  {
    /// \copdydoc ConformalMeshWriter::write_to_mesh
    static void write_to_mesh(ConformalMesh_& target, const AdaptiveMesh_& source, Layer layer)
    {
      // Write faces
      ConformalMeshWriter<AdaptiveMesh_, ConformalMesh_, Shape::Hypercube<2>>::write_to_mesh(target, source, layer);

      auto& v_at_c = target.template get_index_set<3, 0>();
      auto& e_at_c = target.template get_index_set<3, 1>();
      auto& f_at_c = target.template get_index_set<3, 2>();

      Index num_cells = source.get_num_entities(layer, 3);
      for(Index i = 0; i < num_cells; i++)
      {
        for(int j = 0; j < v_at_c.get_num_indices(); j++)
        {
          v_at_c[i][j] = source.template get_face_index<3, 0>(layer, i, j);
        }
        for(int j = 0; j < e_at_c.get_num_indices(); j++)
        {
          e_at_c[i][j] = source.template get_face_index<3, 1>(layer, i, j);
        }
        for(int j = 0; j < f_at_c.get_num_indices(); j++)
        {
          f_at_c[i][j] = source.template get_face_index<3, 2>(layer, i, j);
        }
      }
    }
  };

  /**
   * \brief Visitor class for mesh part projection
   *
   * Collects all mesh entities that stem from the chosen root entity.
   *
   * \author Markus Muegge
   */
  template<typename MeshStorage_>
  struct MeshPartProjectorVisitor
  {
    MeshPartProjectorVisitor(const MeshStorage_& s, Layer l) : storage(s), layer(l)
    {
    }

    void operator()(Intern::VertexKey key)
    {
      if(storage[key].layer <= layer)
      {
        vertices.push_back(key);
      }
    }
    void operator()(Intern::EdgeKey key)
    {
      if(storage[key].layer == layer || storage[key].type.is_zero_refinement())
      {
        edges.push_back(key);
      }
    };
    void operator()(Intern::FaceKey key)
    {
      if(storage[key].layer == layer || storage[key].type.is_zero_refinement())
      {
        faces.push_back(key);
      }
    };
    void operator()(Intern::CellKey key)
    {
      if(storage[key].layer == layer || storage[key].type.is_zero_refinement())
      {
        cells.push_back(key);
      }
    };

    const MeshStorage_& storage;
    Layer layer;

    std::vector<Intern::VertexKey> vertices;
    std::vector<Intern::EdgeKey> edges;
    std::vector<Intern::FaceKey> faces;
    std::vector<Intern::CellKey> cells;
  };

#ifdef DOXYGEN
  /**
   * \brief Algorithm for projecting MeshParts of an underlying mesh onto its adaptive mesh.
   *
   * This is the main declaration of the algorithm. The implementation is split
   * across specializations for individual dimension to allow the algorithm to
   * interact with meshes of any dimension.
   *
   * The algorithm works in two phases. It first collects the references of all
   * mesh entities that stem from one of the entites in the source meshpart.
   * Then it writes a new meshpart using those references.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename TargetMeshType_, typename Shape_>
  struct MeshPartProjector
  {
    using FoundationMesh = typename AdaptiveMesh_::FoundationMeshType;
    using Roots = typename AdaptiveMesh_::MeshRoots;
    using Storage = typename AdaptiveMesh_::MeshStorage;
    using Visitor = MeshPartProjectorVisitor<Storage>;

    /**
     * \brief Collect mesh entities stemming from one of the entities in the source meshpart
     *
     * \param[inout] visitor MeshPartProjectorVisitor for collecting references
     * \param[in] source Source meshpart
     * \param[in] layer The layer to project to
     * \param[in] roots AdaptiveMesh roots
     * \param[in] storage AdaptiveMeshStorage for refinement trees access
     */
    static void collect_keys(
      Visitor& visitor,
      const MeshPart<FoundationMesh>& source,
      Layer layer,
      const Roots& roots,
      const Storage& storage)
    {
    }

    /**
     * \brief Write collected references into a new target meshpart
     *
     * \param[in] visitor Collected entity references
     * \param[out] target Target meshpart to write to
     * \param[in] storage AdaptiveMeshStorage for refinement tree access
     */
    static void write_meshpart(Visitor& visitor, MeshPart<TargetMeshType_>& target, const Storage& storage)
    {
    }
  };
#else
  template<typename AdaptiveMesh_, typename TargetMeshType_, typename Shape_>
  struct MeshPartProjector;
#endif

  /**
   * \brief Algorithm for projecting MeshParts of an underlying mesh onto its adaptive mesh.
   *
   * Specialization for vertices
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename TargetMeshType_>
  struct MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Vertex>
  {
    using FoundationMesh = typename AdaptiveMesh_::FoundationMeshType;
    using Roots = typename AdaptiveMesh_::MeshRoots;
    using Storage = typename AdaptiveMesh_::MeshStorage;
    using Visitor = MeshPartProjectorVisitor<Storage>;

    static void collect_keys(
      Visitor& visitor,
      const MeshPart<FoundationMesh>& source,
      Layer layer,
      const Roots& roots,
      const Storage& storage)
    {
      auto& reg_vertices = source.template get_target_set<0>();
      auto& root_vertices = roots.template by_dim<0>();
      for(Index i(0); i < reg_vertices.get_num_entities(); ++i)
      {
        auto iter = root_vertices.find(reg_vertices[i]);
        if(iter != root_vertices.end())
        {
          storage.walk_subtree(iter->second, visitor, layer);
        }
      }
    }

    static void write_meshpart(Visitor& visitor, MeshPart<TargetMeshType_>& target, const Storage& storage)
    {
      auto& result_verts = target.template get_target_set<0>();
      for(Index i(0); i < visitor.vertices.size(); ++i)
      {
        result_verts[i] = storage.get_index(visitor.vertices[i]);
      }
    }
  };

  /**
   * \brief Algorithm for projecting MeshParts of an underlying mesh onto its adaptive mesh.
   *
   * Specialization for edges
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename TargetMeshType_>
  struct MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Hypercube<1>>
  {
    using FoundationMesh = typename AdaptiveMesh_::FoundationMeshType;
    using Roots = typename AdaptiveMesh_::MeshRoots;
    using Storage = typename AdaptiveMesh_::MeshStorage;
    using Visitor = MeshPartProjectorVisitor<Storage>;

    static void collect_keys(
      Visitor& visitor,
      const MeshPart<FoundationMesh>& source,
      Layer layer,
      const Roots& roots,
      const Storage& storage)
    {
      MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Vertex>::collect_keys(
        visitor,
        source,
        layer,
        roots,
        storage);

      auto& reg_edges = source.template get_target_set<1>();
      auto& root_edges = roots.template by_dim<1>();
      for(Index i(0); i < reg_edges.get_num_entities(); ++i)
      {
        auto iter = root_edges.find(reg_edges[i]);
        if(iter != root_edges.end())
        {
          storage.walk_subtree(iter->second, visitor, layer);
        }
      }
    }

    static void write_meshpart(Visitor& visitor, MeshPart<TargetMeshType_>& target, const Storage& storage)
    {
      MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Vertex>::write_meshpart(visitor, target, storage);

      auto& result_edges = target.template get_target_set<1>();
      for(Index i(0); i < visitor.edges.size(); ++i)
      {
        result_edges[i] = storage.get_index(visitor.edges[i]);
      }
    }
  };

  /**
   * \brief Algorithm for projecting MeshParts of an underlying mesh onto its adaptive mesh.
   *
   * Specialization for faces
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename TargetMeshType_>
  struct MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Hypercube<2>>
  {
    using FoundationMesh = typename AdaptiveMesh_::FoundationMeshType;
    using Roots = typename AdaptiveMesh_::MeshRoots;
    using Storage = typename AdaptiveMesh_::MeshStorage;
    using Visitor = MeshPartProjectorVisitor<Storage>;

    static void collect_keys(
      Visitor& visitor,
      const MeshPart<FoundationMesh>& source,
      Layer layer,
      const Roots& roots,
      const Storage& storage)
    {
      MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Hypercube<1>>::collect_keys(
        visitor,
        source,
        layer,
        roots,
        storage);

      auto& reg_faces = source.template get_target_set<2>();
      auto& root_faces = roots.template by_dim<2>();
      for(Index i(0); i < reg_faces.get_num_entities(); ++i)
      {
        auto iter = root_faces.find(reg_faces[i]);
        if(iter != root_faces.end())
        {
          storage.walk_subtree(iter->second, visitor, layer);
        }
      }
    }

    static void write_meshpart(Visitor& visitor, MeshPart<TargetMeshType_>& target, const Storage& storage)
    {
      MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Hypercube<1>>::write_meshpart(visitor, target, storage);

      auto& result_faces = target.template get_target_set<2>();
      for(Index i(0); i < visitor.faces.size(); ++i)
      {
        result_faces[i] = storage.get_index(visitor.faces[i]);
      }
    }
  };

  /**
   * \brief Algorithm for projecting MeshParts of an underlying mesh onto its adaptive mesh.
   *
   * Specialization for cells
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMesh_, typename TargetMeshType_>
  struct MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Hypercube<3>>
  {
    using FoundationMesh = typename AdaptiveMesh_::FoundationMeshType;
    using Roots = typename AdaptiveMesh_::MeshRoots;
    using Storage = typename AdaptiveMesh_::MeshStorage;
    using Visitor = MeshPartProjectorVisitor<Storage>;

    static void collect_keys(
      Visitor& visitor,
      const MeshPart<FoundationMesh>& source,
      Layer layer,
      const Roots& roots,
      const Storage& storage)
    {
      MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Hypercube<2>>::collect_keys(
        visitor,
        source,
        layer,
        roots,
        storage);

      auto& reg_cells = source.template get_target_set<3>();
      auto& root_cells = roots.template by_dim<3>();
      for(Index i(0); i < reg_cells.get_num_entities(); ++i)
      {
        auto iter = root_cells.find(reg_cells[i]);
        if(iter != root_cells.end())
        {
          storage.walk_subtree(iter->second, visitor, layer);
        }
      }
    }

    static void write_meshpart(Visitor& visitor, MeshPart<TargetMeshType_>& target, const Storage& storage)
    {
      MeshPartProjector<AdaptiveMesh_, TargetMeshType_, Shape::Hypercube<2>>::write_meshpart(visitor, target, storage);

      auto& result_cells = target.template get_target_set<3>();
      for(Index i(0); i < visitor.cells.size(); ++i)
      {
        result_cells[i] = storage.get_index(visitor.cells[i]);
      }
    }
  };

  /**
   * \brief Set of mesh indices
   *
   * \author Markus Muegge
   */
  struct MeshIndexSet
  {
    std::set<Index> vertices;
    std::set<Index> edges;
    std::set<Index> faces;
    std::set<Index> cells;

    template<int dim_>
    const std::set<Index>& by_dim() const
    {
      if constexpr(dim_ == 0)
      {
        return vertices;
      }
      else if constexpr(dim_ == 1)
      {
        return edges;
      }
      else if constexpr(dim_ == 2)
      {
        return faces;
      }
      else
      {
        return cells;
      }
    }
  };

#ifdef DOXYGEN
  /**
   * \brief Alorithm for collecting mesh entities that must be adaptively refined
   *
   * This is the main declaration of the algorithm. The implementation is split
   * across specializations for individual dimension to allow the algorithm to
   * interact with meshes of any dimension.
   *
   * \author Markus Muegge
   */
  template<typename Shape_, typename TemplateSet, typename MeshType_>
  struct EntityCollector
  {
    /**
     * \brief Collect all entities of current dimension with non zero refinement type
     *
     * \param[inout] set Set of collected entities
     * \param[mesh] Underlying mesh
     * \param[in] levels Subdivision levels belonging to underlying mesh
     * \param[in] import_all If true, collect all entities irrespective of refinement type
     */
    static void collect(MeshIndexSet& set, const MeshType_& mesh, const SubdivisionLevels& levels, bool import_all)
    {
    }
  };
#else
  template<typename Shape_, typename TemplateSet, typename MeshType_>
  struct EntityCollector;
#endif

  /**
   * \brief Alorithm for collecting mesh entities that must be adaptively refined
   *
   * Specialization for faces
   *
   * \author Markus Muegge
   */
  template<typename TemplateSet, typename MeshType_>
  struct EntityCollector<Shape::Quadrilateral, TemplateSet, MeshType_>
  {
    /// \copydoc EntityCollector::collect()
    static void collect(MeshIndexSet& set, const MeshType_& mesh, const RefinementField<typename TemplateSet::VertexMarkerType>& levels, bool import_all)
    {
      auto& v_at_f = mesh.template get_index_set<2, 0>();
      auto& e_at_f = mesh.template get_index_set<2, 1>();

      // Collect either all faces (if import_all is true) or those faces with non-zero refinement type
      // Also collects all their sub-entities, i.e. surrounding vertices, and edges.
      for(Index face = 0; face < v_at_f.get_num_entities(); face++)
      {
        auto face_levels = levels.get_tuple(v_at_f[face]);

        if(import_all || ! typename TemplateSet::template RefinementTypeByDim<2>(face_levels).is_zero_refinement())
        {
          set.faces.insert(face);

          set.edges.insert(e_at_f[face][0]);
          set.edges.insert(e_at_f[face][1]);
          set.edges.insert(e_at_f[face][2]);
          set.edges.insert(e_at_f[face][3]);

          set.vertices.insert(v_at_f[face][0]);
          set.vertices.insert(v_at_f[face][1]);
          set.vertices.insert(v_at_f[face][2]);
          set.vertices.insert(v_at_f[face][3]);
        }
      }
    }
  };

  /**
   * \brief Alorithm for collecting mesh entities that must be adaptively refined
   *
   * Specialization for cells
   *
   * \author Markus Muegge
   */
  template<typename TemplateSet, typename MeshType_>
  struct EntityCollector<Shape::Hexahedron, TemplateSet, MeshType_>
  {
    /// \copydoc EntityCollector::collect()
    static void collect(MeshIndexSet& set, const MeshType_& mesh, const RefinementField<typename TemplateSet::VertexMarkerType>& levels, bool import_all)
    {
      auto& f_at_c = mesh.template get_index_set<3, 2>();
      auto& e_at_c = mesh.template get_index_set<3, 1>();
      auto& v_at_c = mesh.template get_index_set<3, 0>();

      // Collect either all cells (if import_all is true) or those cells with non-zero refinement type
      // Also collects all their sub-entities, i.e. surrounding vertices, edges, and faces.
      for(Index cell = 0; cell < v_at_c.get_num_entities(); cell++)
      {
        auto cell_levels = levels.get_tuple(v_at_c[cell]);

        if(import_all || ! typename TemplateSet::template RefinementTypeByDim<3>(cell_levels).is_zero_refinement())
        {
          set.cells.insert(cell);

          for(int i = 0; i < f_at_c.num_indices; i++)
          {
            set.faces.insert(f_at_c[cell][i]);
          }

          for(int i = 0; i < e_at_c.num_indices; i++)
          {
            set.edges.insert(e_at_c[cell][i]);
          }

          for(int i = 0; i < v_at_c.num_indices; i++)
          {
            set.vertices.insert(v_at_c[cell][i]);
          }
        }
      }
    }
  };

  /**
   * \brief Computes the orientation codes for mesh elements.
   *
   * Given a mesh and a cell computes the orientation
   * of a face (lower-dimensional element) of the cell as stored in
   * the mesh compared to how the cell sees that element.
   *
   * \author Markus Muegge
   */
  template<typename MeshType_, int cell_dim, int face_dim>
  int congruency(const MeshType_& mesh, Index cell, int face)
  {
    using ShapeType = typename MeshType_::ShapeType;
    using CellShape = typename Shape::FaceTraits<ShapeType, cell_dim>::ShapeType;
    using FaceShape = typename Shape::FaceTraits<ShapeType, face_dim>::ShapeType;
    using IndexMapping = Intern::FaceIndexMapping<CellShape, face_dim, 0>;
    using Sampler = Intern::CongruencySampler<FaceShape>;

    constexpr int face_verts = Shape::FaceTraits<FaceShape, 0>::count;

    auto& face_at_cell = mesh.template get_index_set<cell_dim, face_dim>();
    auto& v_at_face = mesh.template get_index_set<face_dim, 0>();
    auto& v_at_cell = mesh.template get_index_set<cell_dim, 0>();

    std::array<Index, face_verts> local;

    for(int i = 0; i < face_verts; i++)
    {
      local[i] = v_at_cell[cell][IndexMapping::map(face, i)];
    }
    //return Sampler::compare(v_at_face[face_at_cell[cell][face]], local);
    return Sampler::compare(local, v_at_face[face_at_cell(cell, face)]);
  }

#ifdef DOXYGEN
  /**
   * \brief Algorithm for collecting topologies from the foundation mesh
   *
   * \tparam AdaptiveMeshType_ Type of adaptive mesh to collect topologies for.
   * \tparam topology_dim_ Dimension of topology to collect.
   * \tparam collection_dim_ Current dimension considered by algorithm.
   *
   * This is the main declaration. See specializations for actual implementations.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_, int topology_dim_, int collection_dim_>
  struct FoundationTopologyCollector
  {
  };
#else
  template<typename AdaptiveMeshType_, int topology_dim_, int collection_dim_>
  struct FoundationTopologyCollector;
#endif

  /**
   * \brief Algorithm for collecting topologies from the foundation mesh
   *
   * This is the specialization for collecting vertices.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_, int topology_dim_>
  struct FoundationTopologyCollector<AdaptiveMeshType_, topology_dim_, 0>
  {
    using Topology = typename AdaptiveMeshType_::template ElementTopology<topology_dim_>;
    using TopologyShape = typename Shape::FaceTraits<typename AdaptiveMeshType_::ShapeType, topology_dim_>::ShapeType;
    static constexpr int num_vertices = Shape::FaceTraits<TopologyShape, 0>::count;

    /**
     * \brief Collect vertices of foundation mesh for topologies
     *
     * \param[inout] topology Topology to collect into
     * \param[in] entity Entity of foundation mesh to collect topology of
     * \param[in] a_mesh AdaptiveMesh to collect topology for.
     */
    static void collect(Topology& topology, Index entity, const AdaptiveMeshType_& a_mesh)
    {
      auto& foundation_entity_vertices = a_mesh._foundation_mesh.template get_index_set<topology_dim_, 0>();

      auto& topology_vertices = topology.template by_dim<0>();
      auto& root_vertices = a_mesh._roots.template by_dim<0>();
      for(int vert = 0; vert < num_vertices; vert++)
      {
        Index vertex_index = foundation_entity_vertices(entity, vert);
        topology_vertices[vert] = Intern::OrientedElement<0>(0, root_vertices.at(vertex_index));
      }
    }
  };

  /**
   * \brief Algorithm for collecting topologies from the foundation mesh
   *
   * This is the specialization for collecting edges.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_, int topology_dim_>
  struct FoundationTopologyCollector<AdaptiveMeshType_, topology_dim_, 1>
  {
    using Topology = typename AdaptiveMeshType_::template ElementTopology<topology_dim_>;
    using ShapeType = typename AdaptiveMeshType_::ShapeType;
    using TopologyShape = typename Shape::FaceTraits<ShapeType, topology_dim_>::ShapeType;
    static constexpr int num_edges = Shape::FaceTraits<TopologyShape, 1>::count;

    /**
     * \brief Collect edges of foundation mesh for topologies
     *
     * \param[inout] topology Topology to collect into
     * \param[in] entity Entity of foundation mesh to collect topology of
     * \param[in] a_mesh AdaptiveMesh to collect topology for.
     */
    static void collect(Topology& topology, Index entity, const AdaptiveMeshType_& a_mesh)
    {
      // Recursive call to lower dimension
      FoundationTopologyCollector<AdaptiveMeshType_, topology_dim_, 0>::collect(topology, entity, a_mesh);

      using FoundationMeshType = typename AdaptiveMeshType_::FoundationMeshType;
      auto& foundation_entity_edges = a_mesh._foundation_mesh.template get_index_set<topology_dim_, 1>();

      auto& topology_edges = topology.template by_dim<1>();
      auto& root_edges = a_mesh._roots.template by_dim<1>();
      for(int edge = 0; edge < num_edges; edge++)
      {
        int orientation = congruency<FoundationMeshType, topology_dim_, 1>(a_mesh._foundation_mesh, entity, edge);
        Index edge_index = foundation_entity_edges(entity, (int)edge);
        topology_edges[edge] = Intern::OrientedEdge(orientation, root_edges.at(edge_index));
      }
    }
  };

  /**
   * \brief Algorithm for collecting topologies from the foundation mesh
   *
   * This is the specialization for collecting faces.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_, int topology_dim_>
  struct FoundationTopologyCollector<AdaptiveMeshType_, topology_dim_, 2>
  {
    using Topology = typename AdaptiveMeshType_::template ElementTopology<topology_dim_>;
    using ShapeType = typename AdaptiveMeshType_::ShapeType;
    using TopologyShape = typename Shape::FaceTraits<ShapeType, topology_dim_>::ShapeType;
    static constexpr int num_faces = Shape::FaceTraits<TopologyShape, 2>::count;

    /**
     * \brief Collect faces of foundation mesh for topologies
     *
     * \param[inout] topology Topology to collect into
     * \param[in] entity Entity of foundation mesh to collect topology of
     * \param[in] a_mesh AdaptiveMesh to collect topology for.
     */
    static void collect(Topology& topology, Index entity, const AdaptiveMeshType_& a_mesh)
    {
      // Recursive call to lower dimension
      FoundationTopologyCollector<AdaptiveMeshType_, topology_dim_, 1>::collect(topology, entity, a_mesh);

      using FoundationMeshType = typename AdaptiveMeshType_::FoundationMeshType;
      auto& foundation_entity_faces = a_mesh._foundation_mesh.template get_index_set<topology_dim_, 2>();

      auto& topology_faces = topology.template by_dim<2>();
      auto& root_faces = a_mesh._roots.template by_dim<2>();
      for(int face = 0; face < num_faces; face++)
      {
        int orientation = congruency<FoundationMeshType, topology_dim_, 2>(a_mesh._foundation_mesh, entity, face);
        Index face_index = foundation_entity_faces(entity, face);
        topology_faces[face] = Intern::OrientedFace(orientation, root_faces.at(face_index));
      }
    }
  };
} // namespace FEAT::Geometry::Intern

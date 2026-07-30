// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/adjacency/base.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/fe_interpolator.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/base_header.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_part_ops.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/projection_filter.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/shape.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT::Assembly
{
  /**
   * \brief Factory for MeshParts describing features of the mesh boundary
   *
   * Features are parts of the boundary that should be left unchanged by
   * geometric algorithms because they describe high frequency detail. This
   * factory determines features by determining the angle between the normals of
   * adjacent mesh facets. Entities shared by facets where that angle is larger
   * than some critical angle are categorized as feature entities.
   *
   * On a 2D mesh, this factory compares the normals of adjacent boundary edges.
   * If the angle between the normals is larger than the critical angle, the
   * shared vertex is considered a feature vertex.
   *
   * On a 3D mesh, this factory compares the normals of adjacent boundary faces.
   * If the angle between the normals is larger than the critical angle, the
   * shared edge and the shared vertices are considered feature entities.
   *
   * \author Markus Muegge
   */
  template<typename Mesh_>
  class FeatureMeshPartFactory
  {
  public:
    ///MeshPart type
    using MeshPart = Geometry::MeshPart<Mesh_>;

    /**
     * \brief Creates a MeshPart whose target sets contain all feature entities of the given boundary
     *
     * \param[in] mesh
     * The mesh to operate on. Note that the mesh is not modified despite not
     * being marked as const. This is required because the factory needs to
     * define a transform on the mesh to determine normal vectors.
     *
     * \param[in] boundary
     * MeshPart describing the boundary that features should be determined on
     *
     * \param[in] critical_angle
     * Angle in radians. If the normals of adjacent facets differ by more than
     * this angle, their shared entities are considered feature entities.
     *
     * \returns A MeshPart whose target sets contain all feature entities of the given boundary.
     */
    static MeshPart make(Mesh_& mesh, const MeshPart& boundary, typename Mesh_::CoordType critial_angle = 0.34)
    {
      // Dimension of mesh shapes
      static constexpr int dim = Mesh_::ShapeType::dimension;

      // We are producing 0-dimensional features on 2d meshes and 1-dimensional features on 3d meshes
      static constexpr int feature_dim = dim - 2;

      // Dimension of facets whose normals we compare to determine features
      static constexpr int facet_dim = dim - 1;

      // Trafo type for determining surface normals
      using TrafoType = Trafo::Standard::Mapping<Mesh_>;

      // Normal computer type
      using NormalComputer = Assembly::OuterNormalComputer<TrafoType>;

      // Datatype
      using DataType = typename Mesh_::CoordType;

      // Entities of the feature meshpart
      std::array<std::set<Index>, dim + 1> entities;

      // Precompute orientations and facet_mask
      std::vector<DataType> orientations(mesh.get_num_entities(facet_dim), DataType(0));
      std::vector<Index> facet_mask(mesh.get_num_entities(facet_dim), ~Index(0));
      {
        const auto& facet_tgs = boundary.template get_target_set<facet_dim>();
        Index idx(0);
        for(Index i(0); i < boundary.get_num_entities(facet_dim); i++)
        {
          facet_mask[facet_tgs[i]] = idx;
          idx++;
        }
      }
      NormalComputer::compute_orientations(orientations.data(), facet_mask.data(), mesh);

      const TrafoType trafo(mesh);
      const auto& feature_tgs = boundary.template get_target_set<feature_dim>();

      const auto& feature_at_facet = mesh.template get_index_set<facet_dim, feature_dim>();
      const Adjacency::Graph facet_at_feature(Adjacency::RenderType::transpose, feature_at_facet);

      for(Index i(0); i < boundary.get_num_entities(feature_dim); i++)
      {
        const Index candidate = feature_tgs[i];

        // NOTE: We don't want to assume that the boundary meshpart has a topology.
        // We thus search through the mesh to find the 2 adjacent facets and re-use the facet mask
        // to check if a facet is part of the boundary.
        std::array<Index, 2> boundary_facets;
        Index idx(0);
        for(auto it = facet_at_feature.image_begin(candidate); it != facet_at_feature.image_end(candidate); ++it)
        {
          if(facet_mask[*it] != ~Index(0))
          {
            boundary_facets[idx++] = *it;
          }

        }

        XASSERTM(idx == 2, "FeatureMeshPartFactory: Too many adjacent facets at feature entity " + stringify(candidate));

        // Found adjacent facets, determine angle between normal vectors

        Index facet_a = boundary_facets[0];
        Index facet_b = boundary_facets[1];

        auto normals_a = NormalComputer::compute_oriented_normals(trafo, facet_a);
        auto normals_b = NormalComputer::compute_oriented_normals(trafo, facet_b);

        Tiny::Vector<DataType, dim> normal_a = normals_a[0] * orientations[facet_a];
        Tiny::Vector<DataType, dim> normal_b = normals_b[0] * orientations[facet_b];

        const DataType angle = Math::acos(Tiny::dot(normal_a, normal_b));

        if(angle > critial_angle) // 20 degrees
        {
          add_sub_entities<dim, feature_dim>(entities, mesh, candidate);
        }
      }

      // Create result meshpart by writing collected entities into target sets
      std::array<Index, 4> result_size{0, 0, 0, 0};
      for(Index i(0); i < Index(dim); i++)
      {
        result_size[i] = Index(entities[i].size());
      }
      MeshPart result(result_size.data());

      fill_target_sets<dim, dim>(entities, result);

      return result;
    }

  private:
    /// Write collected entities into proper target sets
    template<int mesh_dim, int dim>
    static void fill_target_sets(std::array<std::set<Index>, mesh_dim + 1>& entities, MeshPart& mesh_part)
    {
      if constexpr(dim >= 0)
      {
        auto& target_set = mesh_part.template get_target_set<dim>();
        Index i(0);
        for(Index entity : entities[dim])
        {
          target_set[i++] = entity;
        }

        fill_target_sets<mesh_dim, dim - 1>(entities, mesh_part);
      }
    }

    /// Add entity and all contained lower-dimensional entities into entity arrays
    template<int mesh_dim, int entity_dim>
    static void add_sub_entities(std::array<std::set<Index>, mesh_dim + 1>& entities, const Mesh_& mesh, Index entity)
    {
      entities[entity_dim].insert(entity);

      if constexpr(entity_dim > 0)
      {
        const auto& idx_set = mesh.template get_index_set<entity_dim, entity_dim - 1>();

        for(int i(0); i < idx_set.num_indices; i++)
        {
          add_sub_entities<mesh_dim, entity_dim - 1>(entities, mesh, idx_set(entity, i));
        }
      }
    }
  };

  /**
   * \brief ProjectionFilter assembler class template
   *
   * Utility class for managing the Neumann boundary and boundary features,
   * determining vertex kinds (Boundary, Feature, Fixed), and assembling
   * ProjectionFilters.
   *
   * \author Markus Muegge
   */
  template<typename Mesh_>
  class ProjectionFilterAssembler
  {
  public:
    /// Type of mesh to assemble filter on
    using MeshType = Mesh_;

  private:
    /// World dimension
    static constexpr int world_dim = MeshType::world_dim;
    /// Shape dimension
    static constexpr int shape_dim = MeshType::shape_dim;
    /// IndexSet type
    using IndexSet = std::set<Index>;
    /// dof-index sets for entities that must be filtered
    std::array<IndexSet, shape_dim + 1> _entities;
    /// dof-index sets for feature entities
    std::array<IndexSet, shape_dim + 1> _feature_entities;

    static constexpr std::array<Index, 4> empty_mesh_part_size = {{0, 0, 0, 0}};

    /// MeshPart describing Neumann boundary
    Geometry::MeshPart<Mesh_> _neuman_zero_boundary{empty_mesh_part_size.data()};
    /// MeshPart describing boundary features
    Geometry::MeshPart<Mesh_> _feature_network{empty_mesh_part_size.data()};

  public:
    ProjectionFilterAssembler() = default;
    ProjectionFilterAssembler(ProjectionFilterAssembler&&) = default;
    ProjectionFilterAssembler& operator=(ProjectionFilterAssembler&&) = default;
    ProjectionFilterAssembler(const ProjectionFilterAssembler&) = delete;
    ProjectionFilterAssembler& operator=(const ProjectionFilterAssembler&) = delete;
    virtual ~ProjectionFilterAssembler() = default;

    /**
     * \brief Adds the given MeshPart to the Neumann boundary
     *
     * \param[in] mesh_part MeshPart to add to Neumann boundary
     */
    void add_mesh_part(const Geometry::MeshPart<MeshType>& mesh_part)
    {
      Intern::UnitAsmWrapper<shape_dim>::merge(_entities.data(), mesh_part);
      _neuman_zero_boundary = Geometry::MeshPartOperations<Mesh_>::meshpart_union(
        _neuman_zero_boundary,
        mesh_part,
        [](auto a, auto /*b*/) { return a; });
    }

    /**
     * \brief Adds the given MeshPart to the boundary features
     *
     * \param[in] mesh_part MeshPart to add to the boundary features
     */
    void add_feature_mesh_part(const Geometry::MeshPart<MeshType>& mesh_part)
    {
      Intern::UnitAsmWrapper<shape_dim>::merge(_feature_entities.data(), mesh_part);
      _feature_network = Geometry::MeshPartOperations<Mesh_>::meshpart_union(
        _feature_network,
        mesh_part,
        [](auto a, auto /*b*/) { return a; });
    }

    /**
     * \brief Assemble a ProjectionFilter
     *
     * \param[out] filter Projection filter to assemble.
     * \param[in] space A \transient reference to the finite element space for which the filter is to be assembled.
     * \param[in] gate An optional pointer to a gate for synchronizing vertex kinds. Must be defined on \c space
     *
     * \warning Currently only implemented for the Space::Lagrange1::Element space.
     */
    template<typename DataType_, typename IndexType_, typename Space_, typename Gate_>
    void
    assemble(LAFEM::ProjectionFilter<MeshType, DataType_, IndexType_>& filter, Space_& space, const Gate_* gate = nullptr) const
    {
      using FilterType = LAFEM::ProjectionFilter<MeshType, DataType_, IndexType_>;

      // We need two things for assembling a ProjectionFilter
      // 1. The coordinates of the DoFs we want to filter
      // 2. A classification of the boundary DoFs into boundary, feature, and fixed DoFs

      // Determine DoF positions
      LAFEM::DenseVectorBlocked<DataType_, IndexType_, MeshType::world_dim> dof_positions =
        _dof_positions<DataType_, IndexType_>(space);

      // Determine valences of DoFs in feature network.
      // WARNING: Only works for Space_ = Space::Lagrange1::Element for now
      // To make this work for at least other Lagrange spaces, one needs to interpolate the valences similarly to the DoF positions
      LAFEM::DenseVector<DataType_, Index> feature_valences = _feature_network_valences<DataType_>(space.get_mesh());

      // Sync valences if possible/required
      if(gate)
      {
        gate->sync_0(feature_valences);
      }

      // Create new filter
      filter = FilterType(&space.get_mesh(), _neuman_zero_boundary.clone(), _feature_network.clone());

      // Assemble filter
      // WARNING: Only works for Space_ = Space::Lagrange1::Element for now
      // To make this work for other spaces, one needs to properly maps DoFs in this routine
      _assemble(filter, feature_valences, dof_positions);
    }

  protected:
    /// Inner assembly helper. Classifies DoFs based on feature network valences and adds DoFs to ProjectionFilter
    template<typename DataType_, typename IndexType_>
    void _assemble(
      LAFEM::ProjectionFilter<MeshType, DataType_, IndexType_>& filter,
      const LAFEM::DenseVector<DataType_, IndexType_>& feature_valences,
      const LAFEM::DenseVectorBlocked<DataType_, IndexType_, MeshType::world_dim>& dof_positions) const
    {
      const auto& target_set = _neuman_zero_boundary.template get_target_set<0>();

      Memory::TypedView<Tiny::Vector<DataType_, MeshType::world_dim>> dof_view = dof_positions.elements_view_r();
      Memory::TypedView<DataType_> valences = feature_valences.elements_view_r();

      for(Index i(0); i < target_set.get_num_entities(); i++)
      {
        // TODO: Handle DoF-mapping
        Index vertex = target_set[i];
        if(Index(valences(vertex)) == 0)
        {
          // A vertex with valence 0 lies somewhere on the boundary, but not on a feature.
          // Tell the filter it is a boundary DoF and can be freely projected
          filter.add_boundary_dof(vertex, dof_view(vertex));
        }
        else if(Index(valences(vertex)) == 2)
        {
          // A vertex with valence 2 is part of a run of feature edges.
          // Tell the filter is is a boundary DoF and needs to projected onto the feature network
          filter.add_feature_dof(vertex, dof_view(vertex));
        }
        else
        {
          // A vertex with any other valence is some special vertex in the feature network.
          // In 2D it is just a feature vertex.
          // In 3D it is either an end of a run of feature edges or an intersection between runs of feature edges.
          // Either way we want the filter to not allow any modification its position
          filter.add_fixed_dof(vertex, dof_view(vertex));
        }
      }
    }

    // Determine valences of vertices in the feature network
    template<typename DataType_>
    LAFEM::DenseVector<DataType_, Index> _feature_network_valences(const MeshType& mesh) const
    {
      LAFEM::DenseVector<DataType_, Index> valences(mesh.get_num_vertices(), 0);
      Memory::TypedView<DataType_> view = valences.elements_view_w();

      if constexpr(Mesh_::world_dim == 3)
      {
        // Feature network consists of edges
        const auto& target_set = _feature_network.template get_target_set<1>();
        const auto& index_set = mesh.template get_index_set<1, 0>();
        for(Index i(0); i < _feature_network.get_num_entities(1); i++)
        {
          view[index_set(target_set[i], 0)]++;
          view[index_set(target_set[i], 1)]++;
        }
      }

      // mesh.get_index_set<0, 0>() should be allowed and be an identity mapping
      // Then these if constexpr would not be needed
      if constexpr(Mesh_::world_dim == 2)
      {
        // Feature network consists of vertices
        const auto& target_set = _feature_network.template get_target_set<0>();
        for(Index i(0); i < _feature_network.get_num_entities(0); i++)
        {
          view[target_set[i]]++;
        }
      }

      return valences;
    }

    /// Determine position of DoFs by interpolating vertex positions from P1 to the given space.
    template<typename DataType_, typename IndexType_, typename Space_>
    LAFEM::DenseVectorBlocked<DataType_, Index, MeshType::world_dim> _dof_positions(Space_& space) const
    {
      using TrafoType = typename Space_::TrafoType;

      MeshType& mesh = space.get_mesh();

      // Copy vertices into a DenseVectorBlocked for interpolation
      LAFEM::DenseVectorBlocked<DataType_, IndexType_, MeshType::world_dim> vertices(mesh.get_num_vertices());
      vertices.format();
      // Copy mesh vertices into dense vector blocked
      {
        const auto& vtx = mesh.get_vertex_set();
        Memory::TypedView<Tiny::Vector<DataType_, MeshType::world_dim>> vertex_view = vertices.elements_view_w();
        for(Index i(0); i < mesh.get_num_vertices(); i++)
        {
          vertex_view[i] = vtx[i];
        }
      }

      // Interpolate from P1 to given space
      Space::Lagrange1::Element<TrafoType> vertex_space(space.get_trafo());
      LAFEM::DenseVectorBlocked<DataType_, IndexType_, MeshType::world_dim> dof_positions(space.get_num_dofs());
      Assembly::FEInterpolator<decltype(vertex_space), Space_>::interpolate(
        dof_positions,
        vertices,
        space,
        vertex_space);

      return dof_positions;
    }
  };
} // namespace FEAT::Assembly

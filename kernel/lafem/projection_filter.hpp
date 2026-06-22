// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/geometry/atlas/surface_mesh.hpp>
#include <kernel/geometry/cgal.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/vertex_set.hpp>
#include <kernel/lafem/base.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/shape.hpp>
#include <kernel/space/dof_mapping_renderer.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT::LAFEM
{
  namespace Intern
  {
    /// Helper functions for determining closest points on a segment and triangle
    template<typename DataType_>
    struct ClosestPointHelpers
    {
      template<int dim>
      using Vec = Tiny::Vector<DataType_, dim>;

      template<int dim_>
      static Vec<dim_> on_segment(const Vec<dim_>& a, const Vec<dim_>& b, const Vec<dim_>& p)
      {
        const Vec<dim_> ab = b - a;
        const Vec<dim_> ap = p - a;

        const DataType_ alpha = Math::clamp(Tiny::dot(ab, ap) / Tiny::dot(ab, ab), 0.0, 1.0);
        return a + alpha * ab;
      }

      static Vec<3> on_triangle(const Vec<3>& a, const Vec<3>& b, const Vec<3>& c, const Vec<3>& p)
      {
        // https://stackoverflow.com/a/74395029
        const Vec<3> ab = b - a;
        const Vec<3> ac = c - a;
        const Vec<3> ap = p - a;

        const DataType_ d1 = dot(ab, ap);
        const DataType_ d2 = dot(ac, ap);
        if (d1 <= 0.0 && d2 <= 0.0)
        {
          return a; //#1
        }

        const Vec<3> bp = p - b;
        const DataType_ d3 = dot(ab, bp);
        const DataType_ d4 = dot(ac, bp);
        if (d3 >= 0.0 && d4 <= d3)
        {
          return b; //#2
        }

        const Vec<3> cp = p - c;
        const DataType_ d5 = dot(ab, cp);
        const DataType_ d6 = dot(ac, cp);
        if (d6 >= 0.0 && d5 <= d6)
        {
          return c; //#3
        }

        const DataType_ vc = (d1 * d4) - (d3 * d2);
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
            const DataType_ v = d1 / (d1 - d3);
            return a + v * ab; //#4
        }

        const DataType_ vb = (d5 * d2) - (d1 * d6);
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        {
            const DataType_ v = d2 / (d2 - d6);
            return a + v * ac; //#5
        }

        const DataType_ va = (d3 * d6) - (d5 * d4);
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
        {
            const DataType_ v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + v * (c - b); //#6
        }

        const DataType_ denom = 1.0 / (va + vb + vc);
        const DataType_ v = vb * denom;
        const DataType_ w = vc * denom;
        return a + v * ab + w * ac; //#0
      }
    };

    /// Closest point calculator
    template<typename FacetShape_, typename MeshType_>
    struct ClosestPoint
    {
      using MeshType = MeshType_;

      using VertexType = typename MeshType_::VertexType;

      static VertexType on_facet(const MeshType& mesh, Index facet, const VertexType& p);
    };

    template<typename MeshType_>
    struct ClosestPoint<Shape::Vertex, MeshType_>
    {
      using MeshType = MeshType_;

      using VertexType = typename MeshType_::VertexType;

      static VertexType on_facet(const MeshType& mesh, Index facet, const VertexType& /*p*/)
      {
        return mesh.get_vertex_set()[facet];
      }
    };

    template<typename MeshType_>
    struct ClosestPoint<Shape::Hypercube<1>, MeshType_>
    {
      using MeshType = MeshType_;

      using VertexType = typename MeshType_::VertexType;

      static VertexType on_facet(const MeshType& mesh, Index facet, const VertexType& p)
      {
        const auto& vertices = mesh.get_vertex_set();
        const auto& indices = mesh.template get_index_set<1, 0>();
        const VertexType a = vertices[indices(facet, 0)];
        const VertexType b = vertices[indices(facet, 1)];
        return ClosestPointHelpers<typename MeshType::CoordType>::on_segment(a, b, p);
      }
    };

    template<typename MeshType_>
    struct ClosestPoint<Shape::Simplex<1>, MeshType_>
    {
      using MeshType = MeshType_;

      using VertexType = typename MeshType_::VertexType;

      static VertexType on_facet(const MeshType& mesh, Index facet, const VertexType& p)
      {
        const auto& vertices = mesh.get_vertex_set();
        const auto& indices = mesh.template get_index_set<1, 0>();
        const VertexType a = vertices[indices(facet, 0)];
        const VertexType b = vertices[indices(facet, 1)];
        return ClosestPointHelpers<typename MeshType::CoordType>::on_segment(a, b, p);
      }
    };

    template<typename MeshType_>
    struct ClosestPoint<Shape::Simplex<2>, MeshType_>
    {
      using MeshType = MeshType_;

      using VertexType = typename MeshType_::VertexType;

      static VertexType on_facet(const MeshType& mesh, Index facet, const VertexType& p)
      {
        const auto& vertices = mesh.get_vertex_set();
        const auto& indices = mesh.template get_index_set<2, 0>();
        const VertexType a = vertices[indices(facet, 0)];
        const VertexType b = vertices[indices(facet, 1)];
        const VertexType c = vertices[indices(facet, 2)];
        return ClosestPointHelpers<typename MeshType::CoordType>::on_triangle(a, b, c, p);
      }
    };

    template<typename MeshType_>
    struct ClosestPoint<Shape::Hypercube<2>, MeshType_>
    {
      using MeshType = MeshType_;

      using VertexType = typename MeshType_::VertexType;

      static VertexType on_facet(const MeshType& mesh, Index facet, const VertexType& p)
      {
        const auto& vertices = mesh.get_vertex_set();
        const auto& indices = mesh.template get_index_set<2, 0>();
        const VertexType a = vertices[indices(facet, 0)];
        const VertexType b = vertices[indices(facet, 1)];
        const VertexType c = vertices[indices(facet, 2)];
        const VertexType d = vertices[indices(facet, 3)];

        auto p1 = ClosestPointHelpers<typename MeshType::CoordType>::on_triangle(a, b, c, p);
        auto p2 = ClosestPointHelpers<typename MeshType::CoordType>::on_triangle(c, b, d, p);

        if((p1 - p).norm_euclid_sqr() < (p2 - p).norm_euclid_sqr())
        {
          return p1;
        }
        else
        {
          return p2;
        }
      }
    };

    /// Given a point \p, finds point on the given MeshPart closest to \c p
    template<typename MeshType_, int dim_ = MeshType_::world_dim - 1>
    typename MeshType_::VertexType project_onto_mesh_part(
      const MeshType_& mesh,
      const Geometry::MeshPart<MeshType_>& mesh_part,
      const typename MeshType_::VertexType& p)
    {
      using DataType = typename MeshType_::CoordType;
      using VertexType = typename MeshType_::VertexType;

      using FacetShape = typename Shape::FaceTraits<typename MeshType_::ShapeType, dim_>::ShapeType;

      const auto& vertices = mesh.get_vertex_set();

      VertexType closest_point = vertices[0];
      if constexpr (dim_ > 0)
      {
        closest_point = project_onto_mesh_part<MeshType_, dim_ - 1>(mesh, mesh_part, p);
      }

      const auto& target_set = mesh_part.template get_target_set<dim_>();
      DataType distance_to_closest = (p - closest_point).norm_euclid_sqr();
      for(Index i(0); i < mesh_part.get_num_entities(dim_); i++)
      {
        VertexType candidate = ClosestPoint<FacetShape, MeshType_>::on_facet(mesh, target_set[i], p);
        DataType distance = (p - candidate).norm_euclid_sqr();
        if(distance < distance_to_closest)
        {
          closest_point = candidate;
          distance_to_closest = distance;
        }
      }

      return closest_point;
    }
  }

  /// Enforces homogenous Neumann boundary conditions in normal direction by projecting onto the boundary
  template<typename MeshType_, typename DT_, typename IT_>
  class ProjectionFilter
  {
  public:
    /// Data type
    using DataType = DT_;

    /// Index type
    using IndexType = IT_;

    /// Mesh type this filter works on
    using MeshType = MeshType_;

    /// Surface mesh chart type
    using MeshPartType = Geometry::MeshPart<MeshType>;

    static constexpr int BlockSize = MeshType::world_dim;

    static constexpr std::array<Index, 4> empty_mesh_part_size = {{0, 0, 0, 0}};

    /// Value type
    using ValueType = Tiny::Vector<DataType, BlockSize>;

    /// Supported vector type
    using VectorType = DenseVectorBlocked<DataType, IndexType, BlockSize>;

    static_assert(BlockSize > 1, "BlockSize has to be >= 2 in ProjectionFilter!");

    /// Type of filter entry
    enum class EntryType
    {
      Boundary, //< DoF must be projected onto boundary
      Feature, //< DoF must be projected onto feature
      Fixed, //< DoF must be set to zero
    };

    /// Filter entry. Stores information on how to project a single DoF
    struct FilterEntry
    {
      /// DoF to filter
      Index dof;
      /// Position of this DoF. Used as origin in projections
      ValueType position;
      /// The type of this DoF, determines how DoF is projected
      EntryType type;

      /// Constructor
      FilterEntry(Index d, ValueType p, EntryType t) : dof(d), position(p), type(t)
      {
      }
    };

  private:
    const MeshType* _mesh;

    /// Surface mesh chart to project onto
    Geometry::MeshPart<MeshType_> _neumann_zero_boundary{empty_mesh_part_size.data()};

    /// Features to project onto
    Geometry::MeshPart<MeshType_> _feature_network{empty_mesh_part_size.data()};

    /// Filter entries
    std::vector<FilterEntry> _entries;

  public:

    /// Default constructor
    ProjectionFilter() = default;

    /// Constructor
    ProjectionFilter(const MeshType* m, MeshPartType&& neumann_boundary, MeshPartType&& feature_network) :
      _mesh(m),
      _neumann_zero_boundary(std::move(neumann_boundary)),
      _feature_network(std::move(feature_network))
    {
    }

    /// Move constructor
    ProjectionFilter(ProjectionFilter&& other) = default;

    /// Move-assignment
    ProjectionFilter& operator=(ProjectionFilter&& other) = default;

    /// Copy constructor
    ProjectionFilter(const ProjectionFilter& other) = delete;

    /// Copy constructor
    ProjectionFilter& operator=(const ProjectionFilter& other) = delete;

    virtual ~ProjectionFilter() = default;

    /**
     * \brief Create a clone of this filter
     */
    ProjectionFilter clone(CloneMode clone_mode = CloneMode::Deep)
    {
      ProjectionFilter other;
      other.clone(*this, clone_mode);
      return other;
    }

    /// Clone a ProjectionFiler. This is always a deep clone.
    void clone(const ProjectionFilter& other, CloneMode /*clone_mode*/ = CloneMode::Deep)
    {
      if(*this != other)
      {
        _mesh = other._mesh;
        _neumann_zero_boundary = other._neumann_zero_boundary.clone();
        _feature_network = other._feature_network.clone();
        _entries = other._entries;
      }
    }

    /// Size of this filter in bytes
    std::size_t bytes() const
    {
      return _feature_network.bytes() + _neumann_zero_boundary.bytes() + _entries.size() * sizeof(FilterEntry);
    }

    /**
     * \brief Add a DoF to the filter
     *
     * \param[in] dof Index of the DoF
     * \param[in] position Position of the DoF
     * \param[in] type Type of the DoF
     *
     * DoF will get projected differently based on their type
     * - EntryType::Boundary will cause DoFs to be projected to the closest point anywhere on the Neumann boundary,
     * - EntryType::Feature will cause DoFs to be projected to the closest point anywhere on a feature,
     * - EntryType::Fixes will cause DoFs to be set to zero.
     */
    void add(Index dof, ValueType position, EntryType type)
    {
      _entries.emplace_back(dof, position, type);
    }

    /// Add a boundary DoF. See ProjectionFilter::add for details.
    void add_boundary_dof(Index dof, ValueType position)
    {
      _entries.emplace_back(dof, position, EntryType::Boundary);
    }

    /// Add a feature DoF. See ProjectionFilter::add for details.
    void add_feature_dof(Index dof, ValueType position)
    {
      _entries.emplace_back(dof, position, EntryType::Feature);
    }

    /// Add a fixed DoF. See ProjectionFilter::add for details.
    void add_fixed_dof(Index dof, ValueType position)
    {
      _entries.emplace_back(dof, position, EntryType::Fixed);
    }

    /// \copydoc UnitFilter::filter_rhs
    void filter_rhs(VectorType& vector) const
    {
      _filter(vector);
    }

    /// \copydoc UnitFilter::filter_sol
    void filter_sol(VectorType& vector) const
    {
      _filter(vector);
    }

    /// \copydoc UnitFilter::filter_def
    void filter_def(VectorType& vector) const
    {
      _filter(vector);
    }

    /// \copydoc UnitFilter::filter_cor
    void filter_cor(VectorType& vector) const
    {
      _filter(vector);
    }

  private:
    void _filter(VectorType& vector) const
    {
      Memory::TypedView<ValueType> view = vector.elements_view_rw();
      for(auto& [dof, position, type] : _entries)
      {
        ValueType velocity = view[dof];

        DataType speed = velocity.norm_euclid();

        DataType scale = 0.01 / speed;

        // Add velocity to location
        ValueType displaced = position + scale * velocity;

        // Project onto surface
        ValueType projected;
        switch(type)
        {
          case EntryType::Boundary:
            projected = Intern::project_onto_mesh_part(*_mesh, _neumann_zero_boundary, displaced);
            break;
          case EntryType::Feature:
            projected = Intern::project_onto_mesh_part(*_mesh, _feature_network, displaced);
            break;
          case EntryType::Fixed:
            projected = position;
            break;
        }

        // Determine new velocity
        view[dof] = (projected - position) * (1.0 / scale);
      }
    }

  };
} // namespace FEAT::LAFEM

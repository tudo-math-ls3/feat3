// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/adjacency/adjactor.hpp"
#include "kernel/adjacency/base.hpp"
#include "kernel/adjacency/graph.hpp"
#include "kernel/shape.hpp"
#include "kernel/util/math.hpp"
#include "kernel/util/tiny_algebra.hpp"
#include <kernel/base_header.hpp>

#include <limits>
#include <optional>

namespace FEAT::Geometry
{
  /**
   * \brief Geometric ray
   *
   * \tparam Vector_ Point and vector type
   *
   * Describes the geometric ray given by
   * \f$r(t) = \text{origin} + t * \text{direction}\f$
   */
  template<typename Vector_>
  struct Ray
  {
    /// Ray origin
    Vector_ origin;

    /// Ray direction
    Vector_ direction;
  };

  template<typename Vector_>
  std::ostream& operator<<(std::ostream& out, const Ray<Vector_>& ray)
  {
    out << "Ray{ origin: " << ray.origin << ", direction: " << ray.direction << "}";
    return out;
  }

  /**
   * \brief Contains information about a ray-primitive intersection
   *
   * \tparam DT_ Data type
   *
   * Contains the distance along the ray of the point(s) where a ray
   * enters and exists some primitive shape.
   *
   * Entry and exit points can differ. For example in case of a coplanar ray-triagle pair.
   * In that case the part of the ray that overlaps the triangle is counted as
   * an intersection and \c t_entry gives the first point of overlap and \c t_exit gives
   * the last point of overlap.
   *
   * In case of a clean intersection, for example between a ray and a non coplanar triangle,
   * \c t_entry and \c t_exit will be identical.
   *
   * The exit point \c t_exit is guaranteed to be identical to \c t_entry or further along the ray.
   */
  template<typename DT_>
  struct IntersectionData
  {
    /// Constructor
    explicit IntersectionData(DT_ t) : t_entry(t), t_exit(t)
    {
    }

    /// Constructor
    IntersectionData(DT_ t1, DT_ t2) : t_entry(t1), t_exit(t2)
    {
      ASSERT(t1 <= t2);
    }

    /// Distance along ray of entry point. Relative to ray direction.
    DT_ t_entry;

    /// Distance along ray of exit point. Relative to ray direction.
    DT_ t_exit;
  };

  template<typename DT_>
  std::ostream& operator<<(std::ostream& out, const IntersectionData<DT_>& d)
  {
    out << "IntersectionData{ t_entry: " << d.t_entry << ", t_exit: " << d.t_exit << "}";
    return out;
  }

  /**
   * \brief Merge two intersection datas
   *
   * Computes the closest entry point and furthest exit point.
   * Use this if you want to treat several raycast across connected shapes as a single raycast.
   */
  template<typename DT_>
  IntersectionData<DT_> merge(const IntersectionData<DT_>& a, const IntersectionData<DT_>& b)
  {
    return {Math::min(a.t_entry, b.t_entry), Math::max(a.t_exit, b.t_exit)};
  };

  /// Normal merge function lifted to optionals
  template<typename DT_>
  std::optional<IntersectionData<DT_>>
  merge(const std::optional<IntersectionData<DT_>>& a, const std::optional<IntersectionData<DT_>>& b)
  {
    if(!a && !b)
    {
      return std::nullopt;
    }
    else if(a && !b)
    {
      return a;
    }
    else if(!a && b)
    {
      return b;
    }
    else
    {
      return std::make_optional(merge(a.value(), b.value()));
    }
  };

  /**
   * \brief Intersection tests between Rays and primitives.
   *
   * \tparam DT_ Data type
   *
   * This class contains intersection test between rays and primitives in various dimensions.
   * Each test is self-contained and implemented as a static member function.
   * The class serves mainly as a container of commonly used types, to keep
   * method signatures neat.
   */
  template<typename DT_>
  class RayIntersectionPrimitives
  {
  public:
    /// Data type
    using DataType = DT_;

    /// Intersection data type
    using IntersectionDataType = IntersectionData<DataType>;

    /// Intersection result type
    using IntersectionResultType = std::optional<IntersectionDataType>;

    /// 2D vector/point type
    using Vector2D = Tiny::Vector<DataType, 2>;

    /// 3D vector/point type
    using Vector3D = Tiny::Vector<DataType, 3>;

    /// Precision
    static constexpr DataType eps = Math::eps<DataType>();

    /**
     * \brief 2D ray-segment intersection
     *
     * \param[in] r Ray
     * \param[in] a First point of segment
     * \param[in] b Second point of segment
     *
     * \returns Entry and exit points, as distances along the ray, if the given ray and segment intersect.
     *
     * Determine intersection between ray \c r and segment \c a, \c b.
     *
     * Correctly detects overlap betweeen colinear rays and segments.
     */
    static IntersectionResultType ray_segment_intersection(const Ray<Vector2D>& r, const Vector2D& a, const Vector2D& b)
    {
      const Vector2D ab = b - a;

      // https://stackoverflow.com/a/565282
      //
      // Solve r.origin + t * r.direction = a + u * ab
      // With x as the 2D cross product:
      // => (r.origin + t * r.direction) x ab = (a + u* ab) x ab
      // => t(r.direction x ab) = (a - r.origin) x ab
      // => t = (a - r.origin) x ab / (r.direction x ab)
      // Analoguous
      // u = (r.origin - a) x r.direction / (ab x r.direction)
      // with v x w = -w x v
      // u = (a - r.origin) x r.direction / (r.direction x ab)

      DataType denom = cross(r.direction, ab);
      DataType t_num = cross(a - r.origin, ab);
      DataType u_num = cross(a - r.origin, r.direction);

      // Are directions of ray and segment parallel
      bool dir_parallel = Math::abs(denom) < eps;

      // Are ray and segment colinear
      bool aligned = Math::abs(u_num) < eps;

      if(dir_parallel && aligned)
      {
        // Solve for overlap between ray and segment

        // Solve r.origin + t * r.direction = a
        // => t * r.direction = a - r.origin
        // => t = (a - r.origin) * r.direction / (r.direction * r.direction)

        // Solve r.origin + t * r.direction = b
        // => t = (a - r.origin) * r.direction / (r.direction * r.direction)

        const DataType rr = Tiny::dot(r.direction, r.direction);
        DataType t1 = Tiny::dot(a - r.origin, r.direction) / rr;
        DataType t2 = Tiny::dot(b - r.origin, r.direction) / rr;

        if(t1 > t2)
        {
          std::swap(t1, t2);
        }

        if(t2 < DataType(0.0))
        {
          // Exit point is behind ray start => no overlap
          return std::nullopt;
        }

        // Discard overlap behind ray origin
        t1 = Math::max(t1, DataType(0.0));
        t2 = Math::max(t2, DataType(0.0));

        return std::make_optional<IntersectionDataType>(t1, t2);
      }
      else if(dir_parallel && !aligned)
      {
        // Lines are parallel and non-intersecting
        return std::nullopt;
      }
      else
      {
        const DataType t = t_num / denom;
        const DataType u = u_num / denom;

        if(t >= DT_(0.0) && DT_(0.0) <= u && u <= DT_(1.0))
        {
          // Intersection in front of ray
          return std::make_optional<IntersectionDataType>(t);
        }
        else
        {
          // Either no intersection or intersection behind ray
          return std::nullopt;
        }
      }
    }

    /**
     * Ray-Triangle intersection
     *
     * \param[in] r Ray
     * \param[in] a First vertex of triangle
     * \param[in] b Second vertex of triangle
     * \param[in] c Third vertex of triangle
     *
     * \returns A ray triangle intersection, if the given ray intersects the given triangle
     *
     * Implementation based on:
     * Tomas Moeller, Ben Trumbore
     * Fast, minimum storage ray-triangle intersection
     * https://doi.org/10.1145/1198555.1198746
     */
    static IntersectionResultType
    ray_triangle_intersection(const Ray<Vector3D>& r, const Vector3D& a, const Vector3D& b, const Vector3D& c)
    {
      // Solve system
      // r.origin + t * r.direction = (1 - u - v)a + ub + uc
      // with u >= 0, v >= 0, u + v <= 1.0
      // Rearranging yields
      // [-r.direction, b - a, c - a] [t, u, v]' = ray.origin - a
      // The below code solves the system using Cramer's rule

      const Vector3D e1 = b - a;
      const Vector3D e2 = c - a;

      Vector3D normal;
      Tiny::cross(normal, r.direction, e2);

      DataType det = Tiny::dot(e1, normal);

      if(det > -eps && det < eps)
      {
        // Ray is parallel to triangle

        if(Math::abs(Tiny::dot(normal, r.origin - a)) > eps)
        {
          // Ray is not coplanar => no intersection
          return std::nullopt;
        }

        // Ray is coplanar. There might be an interval of intersection
        // Reduce to 2D polygon intersection test

        // Using a as origin and e1, e2 as axes for the 2D coordinate system
        // The triangle is then always (0, 0), (1, 0), (0, 1)

        // Project ray into this coordinate system

        // https://math.stackexchange.com/a/1307635
        const auto world_to_plane = [&](const Vector3D& v)
        {
          const DataType alpha = Tiny::dot(normal, cross(e2, v)) / Tiny::dot(normal, cross(e2, e1));
          const DataType beta = Tiny::dot(normal, cross(e1, v)) / Tiny::dot(normal, cross(e1, e2));
          return Vector2D{alpha, beta};
        };

        const auto plane_to_world = [&](const Vector2D& v) { return a + v[0] * e1 + v[1] * e2; };

        Vector2D origin = world_to_plane(r.origin - a);
        Vector2D direction = world_to_plane(r.direction);

        auto i1 = ray_segment_intersection(Ray<Vector2D>{origin, direction}, Vector2D{0.0, 0.0}, Vector2D{1.0, 0.0});
        auto i2 = ray_segment_intersection(Ray<Vector2D>{origin, direction}, Vector2D{1.0, 0.0}, Vector2D{0.0, 1.0});
        auto i3 = ray_segment_intersection(Ray<Vector2D>{origin, direction}, Vector2D{0.0, 1.0}, Vector2D{0.0, 0.0});

        DataType t1(std::numeric_limits<DataType>::max());
        DataType t2(0.0);

        const DataType ray_length = r.direction.norm_euclid();

        auto i = merge(merge(i1, i2), i3);

        if(i)
        {
          const Vector2D p1 = origin + i.value().t_entry * direction;
          const Vector2D p2 = origin + i.value().t_exit * direction;

          t1 = Math::min(t1, (plane_to_world(p1) - r.origin).norm_euclid() / ray_length);
          t2 = Math::max(t2, (plane_to_world(p2) - r.origin).norm_euclid() / ray_length);

          return std::make_optional<IntersectionDataType>(t1, t2);
        }
        else
        {
          return std::nullopt;
        }
      }

      DataType inv_det = DataType(1.0) / det;

      Vector3D tvec = r.origin - a;

      // Calculate u
      const DataType u = Tiny::dot(tvec, normal) * inv_det;
      if(u < DataType(0.0) || u > DataType(1.0))
      {
        return std::nullopt;
      }

      Vector3D qvec;
      Tiny::cross(qvec, tvec, e1);

      // Calculate v
      const DataType v = Tiny::dot(r.direction, qvec) * inv_det;
      if(v < DataType(0.0) || u + v > DataType(1.0))
      {
        return std::nullopt;
      }

      const DataType t = Tiny::dot(e2, qvec) * inv_det;

      if(t > eps)
      {
        return std::make_optional<IntersectionDataType>(t);
      }
      else
      {
        return std::nullopt;
      }
    }

  protected:
    static DataType cross(const Vector2D& a, const Vector2D& b)
    {
      return (a[0] * b[1]) - (a[1] * b[0]);
    };

    static Vector3D cross(const Vector3D& a, const Vector3D& b)
    {
      Vector3D result;
      Tiny::cross(result, a, b);
      return result;
    };
  };

  template<typename MeshType_>
  std::optional<IntersectionData<typename MeshType_::CoordType>>
  facet_raycast(const MeshType_& mesh, const Ray<typename MeshType_::VertexType>& r, Index facet)
  {
    using DataType = typename MeshType_::CoordType;
    using IntersectionData = IntersectionData<DataType>;
    using RIP = RayIntersectionPrimitives<DataType>;

    using ShapeType = typename MeshType_::ShapeType;

    using FacetShapeType = typename Shape::FaceTraits<ShapeType, ShapeType::dimension - 1>::ShapeType;
    static constexpr int facet_dim = FacetShapeType::dimension;
    static constexpr int world_dim = MeshType_::world_dim;

    const auto& v_at_f = mesh.template get_index_set<facet_dim, 0>();
    const auto& vtx = mesh.get_vertex_set();

    if constexpr(std::is_same_v<FacetShapeType, Shape::Vertex> && world_dim == 1)
    {
      // Raycasting on 1d mesh
      DataType origin = r.origin[0];
      DataType direction = r.direction[0];

      DataType facet_coord = mesh.get_vertex_set()[facet][0];

      if(direction > 0 && facet >= origin)
      {
        return std::make_optional<IntersectionData>((facet_coord - origin) / direction);
      }
      else if(direction < 0 && facet <= origin)
      {
        return std::make_optional<IntersectionData>((origin - facet_coord) / direction);
      }
      else
      {
        return std::nullopt;
      }
    }

    if constexpr(std::is_same_v<FacetShapeType, Shape::Hypercube<1>> && world_dim == 2)
    {
      return RIP::ray_segment_intersection(r, vtx[v_at_f(facet, 0)], vtx[v_at_f(facet, 1)]);
    }

    if constexpr(std::is_same_v<FacetShapeType, Shape::Simplex<1>> && world_dim == 2)
    {
      return RIP::ray_segment_intersection(r, vtx[v_at_f(facet, 0)], vtx[v_at_f(facet, 1)]);
    }

    if constexpr(std::is_same_v<FacetShapeType, Shape::Simplex<2>> && world_dim == 3)
    {
      return RIP::ray_triangle_intersection(r, vtx[v_at_f(facet, 0)], vtx[v_at_f(facet, 1)], vtx[v_at_f(facet, 2)]);
    }

    if constexpr(std::is_same_v<FacetShapeType, Shape::Hypercube<2>> && world_dim == 3)
    {
      // Split face into two triangles
      auto i1 = RIP::ray_triangle_intersection(r, vtx[v_at_f(facet, 0)], vtx[v_at_f(facet, 1)], vtx[v_at_f(facet, 2)]);
      auto i2 = RIP::ray_triangle_intersection(r, vtx[v_at_f(facet, 2)], vtx[v_at_f(facet, 1)], vtx[v_at_f(facet, 3)]);

      return merge(i1, i2);
    }

    XABORTM("Unsupported raycast case");
  }

  /**
   * \brief Utility class for casting rays originating from a mesh vertex
   */
  template<typename MeshType_>
  class VertexRaycaster
  {
  public:
    /// Mesh type
    using MeshType = MeshType_;

    /// Data type
    using DataType = typename MeshType::CoordType;

    /// Intersection data type
    using IntersectionDataType = IntersectionData<DataType>;

    /// Intersection result type
    using IntersectionResultType = std::optional<IntersectionDataType>;

    /// Vector type
    using Vector = typename MeshType::VertexType;

    /// Dimension of the mesh shape
    static constexpr int shape_dim = MeshType::shape_dim;

    /// Dimension of the facets to test against
    static constexpr int facet_dim = shape_dim - 1;

  protected:
    /// Reference to mesh
    const MeshType& _mesh;

    /// Precomputes mapping from vertices to adjacent shape_dim-cells.
    Adjacency::Graph _c_at_v;

  public:
    /// Constructor
    explicit VertexRaycaster(const MeshType& mesh) :
      _mesh(mesh),
      _c_at_v(Adjacency::RenderType::transpose, _mesh.template get_index_set<shape_dim, 0>())
    {
    }

    /**
     * \brief Cast a ray from the given vertex in the given direction
     *
     * \param[in] vertex Index of origin vertex
     * \param[in] direction Direction to cast ray in
     *
     * \returns Either nothing or IntersectionData for the closest intersection point by exit point.
     * If multiple intersections with the same closest exit point exits, it chooses the intersection with the closest
     * entry point.
     *
     * Itersections are checked with all facets of all cells adjacent to the given vertex.
     *
     * \note Because the ray always originates from a vertex the returned entry point is extremely unstable.
     * Use the exit point \c t_exit to determine how far a ray reached into the mesh.
     */
    IntersectionResultType cast(const Index vertex, const Vector& direction) const
    {
      const DataType eps = 1e-8;
      const Ray<Vector> r{_mesh.get_vertex_set()[vertex], direction};

      const auto& f_at_c = _mesh.template get_index_set<shape_dim, facet_dim>();
      const Adjacency::CompositeAdjactor adj(_c_at_v, f_at_c);

      IntersectionResultType result(std::nullopt);

      for(auto face = adj.image_begin(vertex); face != adj.image_end(vertex); ++face)
      {
        const IntersectionResultType i = facet_raycast(_mesh, r, *face);

        if(!i)
        {
          continue;
        }

        const IntersectionDataType& hit = i.value();

        if(!result && hit.t_exit > eps)
        {
          result = i;
        }
        else if(result)
        {
          const IntersectionDataType& res = result.value();

          // Update the result if we either found an earlier exit point,
          // or if we found an earlier entry point with the same exit distance.
          if(
            hit.t_exit > eps &&
            (hit.t_exit < res.t_exit || (Math::abs(hit.t_exit - res.t_exit) < eps && hit.t_entry < res.t_entry)))
          {
            result = hit;
          }
        }
      }

      return result;
    }
  };
} // namespace FEAT::Geometry

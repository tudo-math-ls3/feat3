// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/base_header.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/tiny_algebra.hpp>

#if defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)

#include <cstddef>

namespace FEAT
{
  namespace Geometry
  {

    enum class CGALFileMode
    {
      fm_off = 0, /**< OFF */
      fm_obj /**< OBJ */
    };

    /**
     * \brief A feature is an edge-path on a surface mesh, stored as a list of vertex indices.
     *
     * Features generally describe paths of "sharp" edges, i.e. edges with a large difference in surface normal
     * between their adjacent faces. These paths are of special interest when trying to reproduce a surface mesh exactly.
     * Features do not self-intersect. Note that this does not mean that the vertices of a feature acre distinct.
     * Circular features contain one vertex twice to account for all edges of the path.
     */
    using CGALFeature = std::vector<std::uint32_t>;

    /**
     * \brief A FeatureNetwork is a list of features.
     *
     * Each feature is an edge-path. Together these features describe a network of paths across the surface mesh.
     */
    using CGALFeatureNetwork = std::vector<CGALFeature>;

    /**
     * \brief Wrapper for the CGAL Library
     *
     * \tparam DT_ The dataype tp be used.
     *
     * This class acts as a wrapper for a small portion of the CGAL Library.
     * As of now, only the interaction with 3d tetrahedral OFF files is implemented.
     *
     * \note https://www.cgal.org/
     */
    template<typename DT_>
    class CGALWrapper
    {
    public:
      typedef DT_ DataType;
      typedef Tiny::Vector<DataType, 3> PointType;
      typedef Tiny::Matrix<DataType, 3, 3> TransformMatrix;

    private:
      void * _cgal_data;

      /// read in stream in prescibed file format and preprocess search tree for in/out test
      void _parse_mesh(std::istream & file, CGALFileMode file_mode);

    public:
      /// rule of five
      CGALWrapper(const CGALWrapper&) = delete;
      CGALWrapper& operator=(const CGALWrapper&) = delete;
      CGALWrapper(CGALWrapper&&) noexcept;
      CGALWrapper& operator=(CGALWrapper&& other) noexcept;
      virtual ~CGALWrapper();

      /// Create a new CGALWrapper Instance and open the provided file in given format.
      explicit CGALWrapper(const String & filename, CGALFileMode file_mode);

      /// Create a new CGALWrapper Instance and open the provided file stream in given format.
      explicit CGALWrapper(std::istream & file, CGALFileMode file_mode);

      /**
       * \brief Create a new CGALWrapper Instance from vertices and faces
       *
       * \param[in] vertices Array of vertices
       * \param[in] faces Array of faces
       *
       * \returns A CGALWrapper for a surface mesh containing the given vertices and faces
       *
       * \warning CGAL expects faces to be wound counter-clockwise.
       * This constructor will abort the program if that is not the case.
       */
      explicit CGALWrapper(const std::vector<PointType>& vertices, const std::vector<std::array<Index, 3>>& faces);

      /// Check whether a point is inside the Polyhedron defined at objects' construction.
      bool point_inside(DataType x, DataType y, DataType z) const;

      /// Check whether a point is inside or on the boundary of the Polyhedron defined at objects' construction.
      bool point_not_outside(DataType x, DataType y, DataType z) const;

      /// Returns the minimun squared distance between the query point and all input primitives defined at objects' construction.
      DataType squared_distance(DataType x, DataType y, DataType z) const;

      /// Returns the minimum squared distance between the query point and the given feature
      DataType squared_distance_to_feature(const CGALFeature& f, const PointType& p) const;

      /// Returns the point of the surface mesh with the given index
      PointType point(std::uint32_t idx) const;

      /// Returns the nearest point regarding on all input primitives defined at objects' construction.
      PointType closest_point(const PointType& point) const;

      /// Returns the nearest point regarding on all input primitives defined at objects' construction.
      PointType closest_point(DataType x, DataType y, DataType z) const;

      PointType closest_point(const PointType& point, PointType& primitive_grad) const;

      /// Returns the closest point closest to p on the given feature
      PointType closest_point_on_feature(const CGALFeature& f, const PointType& p) const;

      /// Returns a vector of all feature vertices defined at objects' construction.
      CGALFeatureNetwork detect_features(DataType critical_angle);

      /// Applies an affine transformation to the underlying polyhedron and reinitializes the AABB tree
      /// see https://doc.cgal.org/5.5.3/Kernel_23/classCGAL_1_1Aff__transformation__3.html for definition
      /// of transformation matrix and translation
      void transform(const TransformMatrix& trafo_mat, const PointType& translation, DataType scale = DataType(1));

      /// Returns the size in bytes
      std::size_t bytes() const;


    private:
      /// Delete tree, which also requires to delete the inside tester
      void _delete_tree();
      /// initializes tree and inside tester with already initialized polyhedron
      void _init_wrapper();

    }; // class CGALWrapper<typename DT_>

    template<typename MeshType_>
    static CGALWrapper<typename MeshType_::CoordType> cgal_wrapper_from_mesh(const MeshType_& mesh, const MeshPart<MeshType_>& part)
    {
      using MeshType = MeshType_;
      using DT = typename MeshType::CoordType;
      using ShapeType = typename MeshType::ShapeType;
      using FaceType = typename Shape::FaceTraits<ShapeType, 2>::ShapeType;
      using PointType = typename CGALWrapper<DT>::PointType;

      constexpr int vertices_per_face = Shape::FaceTraits<FaceType, 0>::count;

      static_assert(ShapeType::dimension == 3);

      std::vector<PointType> vertices;
      vertices.reserve(part.get_num_entities(0));

      // Prepare vertices
      const auto& vertex_set = mesh.get_vertex_set();
      const TargetSet& vts = part.template get_target_set<0>();
      for(Index i(0); i < part.get_num_entities(0); i++)
      {
        vertices.push_back(vertex_set[vts[i]]);
      }

      // Prepare faces
      std::vector<std::array<Index, 3>> faces;

      const auto& v_at_f = mesh.template get_index_set<2, 0>();
      const TargetSet& fts = part.template get_target_set<2>();
      const Index num_mesh_faces = part.get_num_entities(2);

      if constexpr(std::is_same_v<ShapeType, Shape::Hexahedron>)
      {
        // Quad faces get split into triangles
        faces.reserve(2 * num_mesh_faces);
      }
      if constexpr(std::is_same_v<ShapeType, Shape::Tetrahedron>)
      {
        // Triangle faces do not need to be split
        faces.reserve(num_mesh_faces);
      }

      std::array<Index, vertices_per_face> vs;
      for(Index i(0); i < num_mesh_faces; i++)
      {
        for(int j(0); j < vertices_per_face; j++)
        {
          // Get mesh index
          Index v = v_at_f(fts[i], j);

          // Find corresponding index in vertex target set
          for(Index k(0); k < vts.get_num_entities(); k++)
          {
            if(vts[k] == v)
            {
              vs[j] = k;
              break;
            }
          }
        }

        if constexpr(std::is_same_v<ShapeType, Shape::Hexahedron>)
        {
          faces.push_back({vs[0], vs[1], vs[2]});
          faces.push_back({vs[3], vs[2], vs[1]});
        }
        if constexpr(std::is_same_v<ShapeType, Shape::Tetrahedron>)
        {
          faces.push_back({vs[0], vs[1], vs[2]});
        }
      }

      return CGALWrapper<DT>(vertices, faces);
    }
  } // namespace Geometry
} // namespace FEAT
#endif //defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)

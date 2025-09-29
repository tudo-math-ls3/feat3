// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>

#ifdef FEAT_HAVE_CGAL

#define CGAL_HEADER_ONLY 1

#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(disable: 4244)
#endif // FEAT_COMPILER_MICROSOFT

FEAT_DISABLE_WARNINGS
#define CGAL_NO_DEPRECATION_WARNINGS 1
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h> // deprecated in CGAL 6.x
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
FEAT_RESTORE_WARNINGS

// check if cgal threading support is enabled
#ifdef FEAT_HAVE_OMP
#ifndef CGAL_HAS_THREADS
static_assert(false, "No cgal multithreading support");
#endif

#ifndef BOOST_HAS_THREADS
static_assert(false, "Boost has no threads");
#endif
#endif


#ifdef FEAT_COMPILER_INTEL
#pragma warning disable 3280
#pragma warning disable 1224
#endif //FEAT_COMPILER_INTEL

#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(default: 4244)
#endif // FEAT_COMPILER_MICROSOFT

#undef FEAT_EICKT

#include <kernel/geometry/cgal.hpp>

#ifdef FEAT_HAVE_HALFMATH
#include <kernel/util/half.hpp>
#endif

namespace FEAT::Geometry
{
  template<typename DT_>
  struct CGALTypeWrapper
  {
    typedef DT_ CGALDT_;
    typedef typename CGAL::Simple_cartesian<DT_> K;
    typedef typename K::Point_3 Point_;
    typedef typename CGAL::Surface_mesh<Point_> Polyhedron_;
    typedef typename CGAL::AABB_face_graph_triangle_primitive<Polyhedron_> Primitive_;
    typedef typename CGAL::AABB_traits<K, Primitive_> Traits_;
    typedef typename CGAL::AABB_tree<Traits_> Tree_;
    typedef typename CGAL::Side_of_triangle_mesh<Polyhedron_, K> Point_inside_;
    typedef typename CGAL::Aff_transformation_3<K> Transformation_;
    typedef typename CGAL::AABB_traits<K, Primitive_>::Point_and_primitive_id Point_and_primitive_id_;
  };

  template<typename DT_>
  struct CGALWrapperData
  {
    typedef CGALTypeWrapper<DT_> TW_;
    typename TW_::Polyhedron_ * _polyhedron;
    typename TW_::Tree_ * _tree;
    typename TW_::Point_inside_ * _inside_tester;

    CGALWrapperData() :
      _polyhedron(nullptr),
      _tree(nullptr),
      _inside_tester(nullptr)
    {
    }

    CGALWrapperData(const CGALWrapperData&) = delete;
    CGALWrapperData& operator=(const CGALWrapperData&) = delete;

    void swap(CGALWrapperData& other)
    {
      std::swap(this->_polyhedron, other._polyhedron);
      std::swap(this->_tree, other._tree);
      std::swap(this->_inside_tester, other._inside_tester);
    }

    CGALWrapperData(CGALWrapperData&& other) noexcept
    : _polyhedron(nullptr), _tree(nullptr), _inside_tester(nullptr)
    {
      this->swap(other);
    }

    CGALWrapperData& operator=(CGALWrapperData&& other) noexcept
    {
      this->swap(other);
      return *this;
    }


    ~CGALWrapperData()
    {
      if(_inside_tester)
        delete _inside_tester;
      if(_tree)
        delete _tree;
      if(_polyhedron)
        delete _polyhedron;
    }
  };

  #ifdef FEAT_HAVE_QUADMATH
  template<>
  struct CGALTypeWrapper<__float128>
  {
    typedef double CGALDT_;
    typedef typename CGAL::Simple_cartesian<double> K;
    typedef typename K::Point_3 Point_;
    typedef typename CGAL::Surface_mesh<Point_> Polyhedron_;
    typedef typename CGAL::AABB_face_graph_triangle_primitive<Polyhedron_> Primitive_;
    typedef typename CGAL::AABB_traits<K, Primitive_> Traits_;
    typedef typename CGAL::AABB_tree<Traits_> Tree_;
    typedef typename CGAL::Side_of_triangle_mesh<Polyhedron_, K> Point_inside_;
    typedef typename CGAL::Aff_transformation_3<K> Transformation_;
    typedef typename CGAL::AABB_traits<K, Primitive_>::Point_and_primitive_id Point_and_primitive_id_;
  };
  #endif

  template<typename DT_>
  using PointTypeAlias = typename FEAT::Geometry::template CGALWrapper<DT_>::PointType;

  template<typename DT_>
  using TransformMatrixAlias = typename FEAT::Geometry::template CGALWrapper<DT_>::TransformMatrix;

  template<typename DT_>
  CGALWrapper<DT_>::CGALWrapper(std::istream & file, CGALFileMode file_mode) :
    _cgal_data(nullptr)
  {
    _parse_mesh(file, file_mode);
  }

  template<typename DT_>
  CGALWrapper<DT_>::CGALWrapper(const String & filename, CGALFileMode file_mode) :
    _cgal_data(nullptr)
  {
    std::ifstream file(filename.c_str());
    XASSERTM(file.is_open(), "CGALWrapper: file read error: " + filename + " !");
    _parse_mesh(file, file_mode);
    file.close();
  }

  template<typename MeshType_>
  static CGALWrapperData<typename MeshType_::CoordType>* wrapper_from_mesh(const MeshType_& mesh, const MeshPart<MeshType_>& part)
  {
    static_assert(MeshType_::world_dim == 3);

    using DT = typename MeshType_::CoordType;
    using ShapeType = typename MeshType_::ShapeType;
    using FeatVertexType = typename MeshType_::VertexType;

    using PointType = typename CGALTypeWrapper<DT>::Point_;
    using SurfaceMesh = typename CGALTypeWrapper<DT>::Polyhedron_;
    using VertexIndex = typename SurfaceMesh::Vertex_index;
    using FaceIndex = typename SurfaceMesh::Face_index;

    CGALWrapperData<DT>* result = new CGALWrapperData<DT>;
    result->_polyhedron = new SurfaceMesh;

    // Copy vertices of mesh part into surface mesh
    std::vector<VertexIndex> vertices(part.get_num_entities(0));
    const auto& vertex_set = mesh.get_vertex_set();
    const TargetSet& vertex_target_set = part.template get_target_set<0>();
    const auto& v_at_f = mesh.template get_index_set<2, 0>();
    const TargetSet& face_target_set = part.template get_target_set<2>();

    for(Index i(0); i < part.get_num_entities(0); i++)
    {
      const FeatVertexType v = vertex_set[vertex_target_set[i]];
      vertices[i] = result->_polyhedron->add_vertex(PointType(v[0], v[1], v[2]));
    }

    if constexpr(std::is_same_v<ShapeType, Shape::Hexahedron>)
    {
      std::array<VertexIndex, 4> indices;
      for(Index i(0); i < part.get_num_entities(2); i++)
      {
        for(int j(0); j < 4; j++)
        {
          indices[j] = vertices[vertex_target_set[v_at_f(face_target_set[i], j)]];
        }

        FaceIndex f0 = result->_polyhedron->add_face(indices[0], indices[1], indices[2]);
        FaceIndex f1 = result->_polyhedron->add_face(indices[3], indices[2], indices[1]);

        ASSERTM(f0 != SurfaceMesh::null_face(), "Failed to add face to surface mesh!");
        ASSERTM(f1 != SurfaceMesh::null_face(), "Failed to add face to surface mesh!");
      }
    }

    if constexpr(std::is_same_v<ShapeType, Shape::Tetrahedron>)
    {
      std::array<VertexIndex, 3> indices;

      for(Index i(0); i < part.get_num_entities(2); i++)
      {
        for(int j(0); j < 3; j++)
        {
          indices[j] = vertices[vertex_target_set[v_at_f(face_target_set[i], j)]];
        }

        FaceIndex f0 = result->_polyhedron->add_face(indices[0], indices[1], indices[2]);

        ASSERTM(f0 != SurfaceMesh::null_face(), "Failed to add face to surface mesh!");
      }
    }

    return result;
  }

  template<typename DT_>
  CGALWrapper<DT_>::CGALWrapper(const ConformalMesh<Shape::Hexahedron, 3, DT_>& mesh, const MeshPart<ConformalMesh<Shape::Hexahedron, 3, DT_>>& part) :
    _cgal_data(wrapper_from_mesh(mesh, part))
  {
    _init_wrapper();
  }

  template<typename DT_>
  CGALWrapper<DT_>::CGALWrapper(const ConformalMesh<Shape::Tetrahedron, 3, DT_>& mesh, const MeshPart<ConformalMesh<Shape::Tetrahedron, 3, DT_>>& part) :
    _cgal_data(wrapper_from_mesh(mesh, part))
  {
    _init_wrapper();
  }

  template<typename DT_>
  CGALWrapper<DT_>::CGALWrapper::~CGALWrapper()
  {
    if(_cgal_data)
      delete (CGALWrapperData<DT_>*)_cgal_data;
  }


  template<typename DT_>
  bool CGALWrapper<DT_>::point_inside(DT_ x, DT_ y, DT_ z) const
  {
    typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
    typename CGALTypeWrapper<DT_>::Point_ query{IDT_(x), IDT_(y), IDT_(z)};

    // Determine the side and return true if inside!
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    return (*(cd->_inside_tester))(query) == CGAL::ON_BOUNDED_SIDE;
  }

  template<typename DT_>
  bool CGALWrapper<DT_>::point_not_outside(DT_ x, DT_ y, DT_ z) const
  {
    typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
    typename CGALTypeWrapper<DT_>::Point_ query{IDT_(x), IDT_(y), IDT_(z)};

    // Determine the side and return true if inside!
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    CGAL::Bounded_side test_result = (*(cd->_inside_tester))(query);
    return test_result != CGAL::ON_UNBOUNDED_SIDE;
  }

  template<typename DT_>
  DT_ CGALWrapper<DT_>::squared_distance(DT_ x, DT_ y, DT_ z) const
  {
    typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
    typename CGALTypeWrapper<DT_>::Point_ query{IDT_(x), IDT_(y), IDT_(z)};

    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    return DT_(cd->_tree->squared_distance(query));
  }

  template<typename DT_>
  DT_ CGALWrapper<DT_>::squared_distance_to_feature(const CGALFeature& f, const PointType& p) const
  {
    return (p - closest_point_on_feature(f, p)).norm_euclid_sqr();
  }


  template<typename DT_>
  typename CGALWrapper<DT_>::PointType CGALWrapper<DT_>::closest_point(const PointType& point) const
  {
    return closest_point(point[0], point[1], point[2]);
  }

  template<typename DT_>
  typename CGALWrapper<DT_>::PointType CGALWrapper<DT_>::closest_point(DT_ x, DT_ y, DT_ z) const
  {
    typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
    typename CGALTypeWrapper<DT_>::Point_ query{IDT_(x), IDT_(y), IDT_(z)};
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    typename CGALTypeWrapper<DT_>::Point_ n_point = cd->_tree->closest_point(query);
    return PointType{{DT_(n_point[0]), DT_(n_point[1]), DT_(n_point[2])}};
  }

  template<typename DT_>
  typename CGALWrapper<DT_>::PointType CGALWrapper<DT_>::closest_point(const PointType& point, PointType& primitive_grad) const
  {
    typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
    typename CGALTypeWrapper<DT_>::Point_ query{IDT_(point[0]), IDT_(point[1]), IDT_(point[2])};
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    const typename CGALTypeWrapper<DT_>::Point_and_primitive_id_& cl_p_query = cd->_tree->closest_point_and_primitive(query);
    typename CGALTypeWrapper<DT_>::Polyhedron_::Face_index f = cl_p_query.second;
    auto v = CGAL::Polygon_mesh_processing::compute_face_normal(f,*(cd->_polyhedron));
    for(int i = 0; i < 3; ++i)
      primitive_grad[i] = DT_(v[i]);
    return PointType{{DT_(cl_p_query.first[0]), DT_(cl_p_query.first[1]), DT_(cl_p_query.first[2])}};
  }

  template<typename DT_>
  typename CGALWrapper<DT_>::PointType
  CGALWrapper<DT_>::closest_point_on_feature(const CGALFeature& f, const PointType& query) const
  {
    using PolyhedronType = typename CGALTypeWrapper<DT_>::Polyhedron_;
    using VertexIndex = typename PolyhedronType::Vertex_index;

    CGALWrapperData<DT_>* cd = (CGALWrapperData<DT_>*)_cgal_data;
    PolyhedronType* poly = cd->_polyhedron;

    const auto to_feat_point = [](const auto& p)
    {
      return PointType{p.x(), p.y(), p.z()};
    };

    PointType closest = to_feat_point(poly->point(VertexIndex(f[0])));
    DT_ distance_to_closest = (query - closest).norm_euclid_sqr();

    const auto closest_point_on_segment = [&](PointType start, PointType end, PointType q)
    {
      PointType segment = end - start;
      PointType segment_dir = segment.normalize();

      // Orthogonal projection of (query - start) onto the segment
      DT_ coefficient = Math::clamp(dot(segment_dir, q - start), DT_(0.0), DT_(1.0));

      return start + coefficient * segment_dir;
    };

    for(std::size_t i(0); i < f.size() - 1; i++)
    {
      PointType segment_start = to_feat_point(poly->point(VertexIndex(f[i])));
      PointType segment_end = to_feat_point(poly->point(VertexIndex(f[i + 1])));

      PointType closest_on_segment = closest_point_on_segment(segment_start, segment_end, query);

      DT_ distance_to_segment = (query - closest_on_segment).norm_euclid_sqr();
      if(distance_to_segment < distance_to_closest)
      {
        closest = closest_on_segment;
        distance_to_closest = distance_to_segment;
      }
    }

    return closest;
  }

  template<typename DT_>
  typename CGALWrapper<DT_>::PointType CGALWrapper<DT_>::point(std::uint32_t idx) const
  {
    using PolyhedronType = typename CGALTypeWrapper<DT_>::Polyhedron_;
    using VertexIndex = typename PolyhedronType::Vertex_index;

    const auto to_feat_point = [](const auto& p)
    {
      return PointType{p.x(), p.y(), p.z()};
    };

    CGALWrapperData<DT_>* cd = (CGALWrapperData<DT_>*)_cgal_data;

    return to_feat_point(cd->_polyhedron->point(VertexIndex(idx)));
  }

  template<typename DT_>
  void CGALWrapper<DT_>::_parse_mesh(std::istream & file, CGALFileMode file_mode)
  {
    delete (CGALWrapperData<DT_>*)_cgal_data;

    _cgal_data = new CGALWrapperData<DT_>;
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;

    cd->_polyhedron = new typename CGALTypeWrapper<DT_>::Polyhedron_;

    bool status(false);
    switch(file_mode)
    {
      case CGALFileMode::fm_off:
        status = CGAL::IO::read_OFF(file, *(cd->_polyhedron));
        XASSERTM(status == true, "CGAL::IO::read_OFF: read/parse error !");
        break;
      case CGALFileMode::fm_obj:
        status = CGAL::IO::read_OBJ(file, *(cd->_polyhedron));
        XASSERTM(status == true, "CGAL::IO::read_OBJ: read/parse error !");
        break;
      default:
        XASSERTM(false, "CGAL FileMode not supported!");
    }

    _init_wrapper();

  }

  template<typename DT_>
  CGALFeatureNetwork CGALWrapper<DT_>::detect_features(DT_ critical_angle)
  {
    using PolyhedronType = typename CGALTypeWrapper<DT_>::Polyhedron_;
    using VertexIndex = typename PolyhedronType::Vertex_index;
    using EdgeIndex = typename PolyhedronType::Edge_index;
    using HalfedgeIndex = typename PolyhedronType::Halfedge_index;
    using EdgeIsSharpMapType = typename PolyhedronType::template Property_map<EdgeIndex, bool>;

    CGALFeatureNetwork result;

    auto* cd = static_cast<CGALWrapperData<DT_>*>(_cgal_data);
    PolyhedronType* poly = cd->_polyhedron;

    ///////////////////////
    // Detect sharp edges
    ///////////////////////

    std::pair<EdgeIsSharpMapType, bool> map = poly->template add_property_map<EdgeIndex, bool>("e:is_sharp");

    if(map.second)
    {
      // Property was newly created. Fill it with proper data.
      CGAL::Polygon_mesh_processing::detect_sharp_edges(*poly, critical_angle, map.first);
    }

    /////////////////////////////////
    // Find feature starting points
    /////////////////////////////////

    std::vector<VertexIndex> seeds;

    for(VertexIndex v : poly->vertices())
    {
      Index incident_feature_edges(0);

      for(HalfedgeIndex he : poly->halfedges_around_target(poly->halfedge(v)))
      {
        EdgeIndex edge = poly->edge(he);
        if(map.first[edge])
        {
          incident_feature_edges++;
        }
      }

      // At least one incident feature that is not just a feature passing through this vertex
      if(incident_feature_edges != 0 && incident_feature_edges != 2)
      {
        seeds.push_back(v);
      }
    }

    /////////////////////////////////
    // Flood fill for feature lines
    /////////////////////////////////

    std::vector<bool> visited(poly->number_of_edges(), false);
    const auto flood_feature = [&](VertexIndex v, HalfedgeIndex h)
    {
      CGALFeature feature;
      feature.push_back(static_cast<std::uint32_t>(v));
      VertexIndex current = poly->source(h);
      feature.push_back(static_cast<std::uint32_t>(current));

      visited[poly->edge(h)] = true;

      while(current != v && std::find(seeds.begin(), seeds.end(), current) == seeds.end())
      {
        for(HalfedgeIndex he : poly->halfedges_around_target(poly->halfedge(current)))
        {
          EdgeIndex e = poly->edge(he);
          if(!map.first[e] || visited[e])
          {
            // Edge is either not sharp or has already been visited.
            // Either way it is not the way forward
            continue;
          }

          current = poly->source(he);
          feature.push_back(static_cast<std::uint32_t>(current));
          visited[e] = true;
        }
      }

      return feature;
    };

    /////////////////////
    // Fill in features
    /////////////////////

    for(VertexIndex seed : seeds)
    {
      for(HalfedgeIndex he : poly->halfedges_around_target(poly->halfedge(seed)))
      {
        EdgeIndex e = poly->edge(he);
        if(map.first[e] && !visited[e])
        {
          // Edge is sharp and unvisited. It must be the start of a new feature.
          result.push_back(flood_feature(seed, he));
        }
      }
    }

    // On isolated circular features, all vertices have exactly two incident sharp edges.
    // They thus contain to seed vertices and have to be handled specially.

    for(EdgeIndex e : poly->edges())
    {
      if(map.first[e] && !visited[e])
      {
        HalfedgeIndex he = poly->halfedge(e);
        VertexIndex start = poly->target(he);

        result.push_back(flood_feature(start, he));
      }
    }

    return result;
  }

  template<typename DT_>
  void CGALWrapper<DT_>::transform(const TransformMatrix& scale_rot, const PointType& translation, DT_ scale)
  {
    typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    //delete tree and _inside_tester
    _delete_tree();
    //build affine transformation
    typename CGALTypeWrapper<DT_>::Transformation_ trafo_mat{(IDT_)(scale_rot[0][0]), (IDT_)(scale_rot[0][1]), (IDT_)(scale_rot[0][2]), (IDT_)(translation[0]),
                              (IDT_)(scale_rot[1][0]), (IDT_)(scale_rot[1][1]), (IDT_)(scale_rot[1][2]), (IDT_)(translation[1]),
                              (IDT_)(scale_rot[2][0]), (IDT_)(scale_rot[2][1]), (IDT_)(scale_rot[2][2]), (IDT_)(translation[2]),
                              IDT_(scale)};
    //apply affine transformation on polyhedron
    CGAL::Polygon_mesh_processing::transform(trafo_mat, *(cd->_polyhedron));
    //create new trees
    _init_wrapper();


  }

  template<typename DT_>
  std::size_t CGALWrapper<DT_>::bytes() const
  {
    auto* cd = static_cast<CGALWrapperData<DT_>*>(_cgal_data);

    // CGALs surface mesh uses std::uint32_t to store it indices
    constexpr std::size_t index_size = sizeof(std::uint32_t);

    // CGAL stores an halfedge-index, the vertex coordinates, and a removed flag for each vertex
    constexpr std::size_t bytes_per_vertex = index_size + sizeof(typename CGALTypeWrapper<DT_>::Point_) + sizeof(bool);

    // CGAL stores 4 indices per halfedge
    constexpr std::size_t bytes_per_halfedge = 4 * index_size;

    // CGAL stores just the removed flag for edges
    constexpr std::size_t bytes_per_edge = sizeof(bool);

    // CGAL stores 1 index per face and the removed flag
    constexpr std::size_t bytes_per_face = index_size + sizeof(bool);

    std::size_t result = 0;

    if(cd->_polyhedron)
    {
      result += cd->_polyhedron->number_of_vertices() * bytes_per_vertex +
                cd->_polyhedron->number_of_halfedges() * bytes_per_halfedge +
                cd->_polyhedron->number_of_edges() * bytes_per_edge +
                cd->_polyhedron->number_of_faces() * bytes_per_face;
    }

    if(cd->_tree != nullptr)
    {
      // See https://doc.cgal.org/latest/AABB_tree/index.html, section 4.2 Memory
      constexpr std::size_t bytes_per_primitive = 150;

      result += cd->_tree->size() * bytes_per_primitive;
    }

    // Inside tester stores no meaningful own data

    return result;
  }

  template<typename DT_>
  void CGALWrapper<DT_>::_delete_tree()
  {
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    XASSERTM(cd->_polyhedron != nullptr, "ERROR: Polyhedron is not initialized!");
    delete cd->_inside_tester;
    cd->_inside_tester = nullptr;
    delete cd->_tree;
    cd->_tree = nullptr;
  }

  template<typename DT_>
  void CGALWrapper<DT_>::_init_wrapper()
  {
    CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
    XASSERTM(cd->_polyhedron != nullptr, "ERROR: Polyhedron is not initialized!");
    XASSERTM(cd->_tree == nullptr && cd->_inside_tester == nullptr, "ERROR: Tree or Inside Tester are already initialized");

      // Construct AABB tree with a KdTree
    cd->_tree = new typename CGALTypeWrapper<DT_>::Tree_(faces(*(cd->_polyhedron)).first, faces(*(cd->_polyhedron)).second, *(cd->_polyhedron));
    cd->_tree->accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    cd->_inside_tester = new typename CGALTypeWrapper<DT_>::Point_inside_(*(cd->_tree));
  }

  template<typename DT_>
  CGALWrapper<DT_>::CGALWrapper(CGALWrapper<DT_>&& other) noexcept
  : _cgal_data(nullptr)
  {
    std::swap(this->_cgal_data, other._cgal_data);
  }

  template<typename DT_>
  CGALWrapper<DT_>& CGALWrapper<DT_>::operator=(CGALWrapper<DT_>&& other) noexcept
  {
    std::swap(this->_cgal_data, other._cgal_data);
    return *this;
  }

  // explicitly instantiate templates for all sensible datatype

  template class CGALWrapper<double>;
  template class CGALWrapper<float>;
  #ifdef FEAT_HAVE_HALFMATH
  template class CGALWrapper<Half>;
  #endif
  #ifdef FEAT_HAVE_QUADMATH
  template class CGALWrapper<__float128>;
  #endif
} // namespace FEAT::Geometry

#elif !defined(DOXYGEN)

// Dummy class instance to silence ipo linker optimization warnings about empty libgeometry
class ipo_foobar_geometry_cgal
{
public:
  int i;
  ipo_foobar_geometry_cgal() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo_geometry_cgal;

#endif // FEAT_HAVE_CGAL

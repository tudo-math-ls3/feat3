// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <iostream>
#include <fstream>

#ifdef FEAT_HAVE_CGAL

#define CGAL_HEADER_ONLY 1

#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(disable: 4244)
#endif // FEAT_COMPILER_MICROSOFT

FEAT_DISABLE_WARNINGS
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
FEAT_RESTORE_WARNINGS


#ifdef FEAT_COMPILER_INTEL
#pragma warning disable 3280
#pragma warning disable 1224
#endif //FEAT_COMPILER_INTEL

#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(default: 4244)
#endif // FEAT_COMPILER_MICROSOFT

#include <kernel/geometry/cgal.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_;
typedef CGAL::Polyhedron_3<K> Polyhedron_;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_> Primitive_;
typedef CGAL::AABB_traits<K, Primitive_> Traits_;
typedef CGAL::AABB_tree<Traits_> Tree_;
typedef CGAL::Side_of_triangle_mesh<Polyhedron_, K> Point_inside_;
typedef CGAL::Aff_transformation_3<K> Transformation_;
typedef CGAL::AABB_traits<K, Primitive_>::Point_and_primitive_id Point_and_primitive_id_;

struct CGALWrapperData
{
  Polyhedron_ * _polyhedron;
  Tree_ * _tree;
  Point_inside_ * _inside_tester;

  CGALWrapperData() :
    _polyhedron(nullptr),
    _tree(nullptr),
    _inside_tester(nullptr)
  {
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

typedef typename FEAT::Geometry::CGALWrapper::PointType PointType;
typedef typename FEAT::Geometry::CGALWrapper::TransformMatrix TransformMatrix;

FEAT::Geometry::CGALWrapper::CGALWrapper(std::istream & file, CGALFileMode file_mode) :
  _cgal_data(nullptr)
{
  _parse_mesh(file, file_mode);
}

FEAT::Geometry::CGALWrapper::CGALWrapper(const String & filename, CGALFileMode file_mode) :
  _cgal_data(nullptr)
{
  std::ifstream file(filename.c_str());
  XASSERTM(file.is_open(), "CGALWrapper: file read error: " + filename + " !");
  _parse_mesh(file, file_mode);
  file.close();
}

FEAT::Geometry::CGALWrapper::~CGALWrapper()
{
  if(_cgal_data)
    delete (CGALWrapperData*)_cgal_data;
}

bool FEAT::Geometry::CGALWrapper::point_inside(double x, double y, double z) const
{
  Point_ query(x, y, z);

  // Determine the side and return true if inside!
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;
  return (*(cd->_inside_tester))(query) == CGAL::ON_BOUNDED_SIDE;
}

double FEAT::Geometry::CGALWrapper::squared_distance(double x, double y, double z) const
{
  Point_ query(x, y, z);

  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;
  return (double)cd->_tree->squared_distance(query);
}

PointType FEAT::Geometry::CGALWrapper::closest_point(const PointType& point) const
{
  return closest_point(point[0], point[1], point[2]);
}

PointType FEAT::Geometry::CGALWrapper::closest_point(double x, double y, double z) const
{
  Point_ query(x, y, z);
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;
  Point_ n_point = cd->_tree->closest_point(query);
  return PointType{{(double)n_point[0], (double)n_point[1], (double)n_point[2]}};
}

PointType FEAT::Geometry::CGALWrapper::closest_point(const PointType& point, PointType& primitive_grad) const
{
  Point_ query{(double)point[0], (double)point[1], (double)point[2]};
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;
  const Point_and_primitive_id_& cl_p_query = cd->_tree->closest_point_and_primitive(query);
  Polyhedron_::Face_handle f = cl_p_query.second;
  auto v = CGAL::Polygon_mesh_processing::compute_face_normal(f,*(cd->_polyhedron));
  for(int i = 0; i < 3; ++i)
    primitive_grad[i] = (double)v[i];
  return PointType{{(double)cl_p_query.first[0], (double)cl_p_query.first[1], (double)cl_p_query.first[2]}};
}

void FEAT::Geometry::CGALWrapper::_parse_mesh(std::istream & file, CGALFileMode file_mode)
{
  delete (CGALWrapperData*)_cgal_data;

  _cgal_data = new CGALWrapperData;
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;

  cd->_polyhedron = new Polyhedron_;

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

void FEAT::Geometry::CGALWrapper::transform(const TransformMatrix& scale_rot, const PointType& translation, double scale)
{
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;
  //delete tree and _inside_tester
  _delete_tree();
  //build affine transformation
  Transformation_ trafo_mat(scale_rot[0][0], scale_rot[0][1], scale_rot[0][2], translation[0],
                            scale_rot[1][0], scale_rot[1][1], scale_rot[1][2], translation[1],
                            scale_rot[2][0], scale_rot[2][1], scale_rot[2][2], translation[2],
                                                                                                      scale);
  //apply affine transformation on polyhedron
  CGAL::Polygon_mesh_processing::transform(trafo_mat, *(cd->_polyhedron));
  //create new trees
  _init_wrapper();


}

/// TODO: Implement me
std::size_t FEAT::Geometry::CGALWrapper::bytes() const
{
  return 0;
}

void FEAT::Geometry::CGALWrapper::_delete_tree()
{
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;
  XASSERTM(cd->_polyhedron != nullptr, "ERROR: Polyhedron is not initialized!");
  delete cd->_inside_tester;
  cd->_inside_tester = nullptr;
  delete cd->_tree;
  cd->_tree = nullptr;
}

void FEAT::Geometry::CGALWrapper::_init_wrapper()
{
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;
  XASSERTM(cd->_polyhedron != nullptr, "ERROR: Polyhedron is not initialized!");
  XASSERTM(cd->_tree == nullptr && cd->_inside_tester == nullptr, "ERROR: Tree or Inside Tester are already initialized");

    // Construct AABB tree with a KdTree
  cd->_tree = new Tree_(faces(*(cd->_polyhedron)).first, faces(*(cd->_polyhedron)).second, *(cd->_polyhedron));
  cd->_tree->accelerate_distance_queries();
  // Initialize the point-in-polyhedron tester
  cd->_inside_tester = new Point_inside_(*(cd->_tree));
}

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

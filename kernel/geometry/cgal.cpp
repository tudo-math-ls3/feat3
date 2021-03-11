// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <iostream>
#include <fstream>

#ifdef FEAT_HAVE_CGAL

#define CGAL_HEADER_ONLY 1

FEAT_DISABLE_WARNINGS
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>
FEAT_RESTORE_WARNINGS


#ifdef FEAT_COMPILER_INTEL
#pragma warning disable 3280
#pragma warning disable 1224
#endif //FEAT_COMPILER_INTEL

#include <kernel/geometry/cgal.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_;
typedef CGAL::Polyhedron_3<K> Polyhedron_;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_> Primitive_;
typedef CGAL::AABB_traits<K, Primitive_> Traits_;
typedef CGAL::AABB_tree<Traits_> Tree_;
typedef CGAL::Side_of_triangle_mesh<Polyhedron_, K> Point_inside_;

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
    delete _inside_tester;
    delete _tree;
    delete _polyhedron;
  }
};

FEAT::Geometry::CGALWrapper::CGALWrapper(String filename) :
  _cgal_data(nullptr)
{
  std::ifstream file(filename.c_str());
  XASSERTM(file.is_open(), "CGAL::read_off: file read error!");

  _parse_off_data(file);
  file.close();
}

FEAT::Geometry::CGALWrapper::CGALWrapper(std::istream & file) :
  _cgal_data(nullptr)
{
  _parse_off_data(file);
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

void FEAT::Geometry::CGALWrapper::_parse_off_data(std::istream & file)
{
  delete (CGALWrapperData*)_cgal_data;

  _cgal_data = new CGALWrapperData;
  CGALWrapperData * cd = (CGALWrapperData*)_cgal_data;

  cd->_polyhedron = new Polyhedron_;
  bool status = read_off(file, *(cd->_polyhedron));
  XASSERTM(status == true, "CGAL::read_off: read error!");

  // Construct AABB tree with a KdTree
  cd->_tree = new Tree_(faces(*(cd->_polyhedron)).first, faces(*(cd->_polyhedron)).second, *(cd->_polyhedron));
  cd->_tree->accelerate_distance_queries();
  // Initialize the point-in-polyhedron tester
  cd->_inside_tester = new Point_inside_(*(cd->_tree));
}

FEAT::Geometry::CGALWrapper::~CGALWrapper()
{
  delete (CGALWrapperData*)_cgal_data;
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

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <iostream>
#include <fstream>

#ifdef FEAT_HAVE_CGAL

#ifdef FEAT_COMPILER_INTEL
#pragma warning disable 3280
#pragma warning disable 1224
#endif //FEAT_COMPILER_INTEL

#include <kernel/geometry/cgal.hpp>

using namespace CGAL;

FEAT::Geometry::CGALWrapper::CGALWrapper(String filename) :
  _polyhedron(nullptr),
  _tree(nullptr),
  _inside_tester(nullptr)
{
  std::ifstream file(filename.c_str());
  XASSERTM(file.is_open(), "CGAL::read_off: file read error!");

  _parse_off_data(file);
  file.close();
}

FEAT::Geometry::CGALWrapper::CGALWrapper(std::istream & file) :
  _polyhedron(nullptr),
  _tree(nullptr),
  _inside_tester(nullptr)
{
  _parse_off_data(file);
}

bool FEAT::Geometry::CGALWrapper::point_inside(double x, double y, double z) const
{
  Point query(x, y, z);

  // Determine the side and return true if inside!
  return (*_inside_tester)(query) == CGAL::ON_BOUNDED_SIDE;
}

double FEAT::Geometry::CGALWrapper::squared_distance(double x, double y, double z) const
{
  Point query(x, y, z);

  return (double)_tree->squared_distance(query);
}

void FEAT::Geometry::CGALWrapper::_parse_off_data(std::istream & file)
{
  delete _polyhedron;
  delete _tree;
  delete _inside_tester;

  _polyhedron = new Polyhedron;
  bool status = read_off(file, *_polyhedron);
  XASSERTM(status == true, "CGAL::read_off: read error!");

  // Construct AABB tree with a KdTree
  _tree = new Tree(faces(*_polyhedron).first, faces(*_polyhedron).second, *_polyhedron);
  _tree->accelerate_distance_queries();
  // Initialize the point-in-polyhedron tester
  _inside_tester = new Point_inside(*_tree);
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

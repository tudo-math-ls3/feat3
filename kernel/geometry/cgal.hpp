#pragma once
#ifndef KERNEL_GEOMETRY_CGAL_HPP
#define KERNEL_GEOMETRY_CGAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

#if defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)

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

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Wrapper for the CGAL Library
     *
     *
     * This class acts as a wrapper for a small portion of the CGAL Library.
     * As of now, only the interaction with 3d tetrahedral OFF files is implemented.
     *
     * \note https://www.cgal.org/
     */
    class CGALWrapper
    {
      private:
      typedef CGAL::Simple_cartesian<double> K;
      typedef K::Point_3 Point;
      typedef CGAL::Polyhedron_3<K> Polyhedron;
      typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
      typedef CGAL::AABB_traits<K, Primitive> Traits;
      typedef CGAL::AABB_tree<Traits> Tree;
      typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;
      Polyhedron * _polyhedron;
      Tree * _tree;
      Point_inside * _inside_tester;

      /// read in stream in off file format and preprocess search tree for in/out test
      void _parse_off_data(std::istream & file);

      public:
      /// Create a new CGALWrapper Instance and open the provided off file.
      CGALWrapper(String filename);

      /// Create a new CGALWrapper Instance and open the provided off (file-) stream.
      CGALWrapper(std::istream & file);

      /// Check whether a point is inside the Polyhedron defined by the aforementioned off file.
      bool point_inside(double x, double y, double z) const;

      /// Destructor
      ~CGALWrapper()
      {
        delete _inside_tester;
        delete _tree;
        delete _polyhedron;
      }
    };

  } // namespace Geometry
} // namespace FEAT
#endif //defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)
#endif //KERNEL_GEOMETRY_CGAL_HPP

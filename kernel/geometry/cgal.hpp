// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include <cstddef>
#ifndef KERNEL_GEOMETRY_CGAL_HPP
#define KERNEL_GEOMETRY_CGAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/tiny_algebra.hpp>

#if defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)

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
    public:
      typedef Tiny::Vector<double, 3> PointType;
      typedef Tiny::Matrix<double, 3, 3> TransformMatrix;
    private:
      void * _cgal_data;

      /// read in stream in prescibed file format and preprocess search tree for in/out test
      void _parse_mesh(std::istream & file, CGALFileMode file_mode);

    public:
      /// Create a new CGALWrapper Instance and open the provided file in given format.
      explicit CGALWrapper(const String & filename, CGALFileMode file_mode);

      /// Create a new CGALWrapper Instance and open the provided file stream in given format.
      explicit CGALWrapper(std::istream & file, CGALFileMode file_mode);

      /// Destructor
      virtual ~CGALWrapper();

      /// Check whether a point is inside the Polyhedron defined at objects' construction.
      bool point_inside(double x, double y, double z) const;

      /// Returns the minimun squared distance between the query point and all input primitives defined at objects' construction.
      double squared_distance(double x, double y, double z) const;

      /// Returns the nearest point regarding on all input primitives defined at objects' construction.
      PointType closest_point(const PointType& point) const;

      /// Returns the nearest point regarding on all input primitives defined at objects' construction.
      PointType closest_point(double x, double y, double z) const;

      PointType closest_point(const PointType& point, PointType& primitive_grad) const;

      /// Applies an affine transformation to the underlying polyhedron and reinitializes the AABB tree
      /// see https://doc.cgal.org/5.5.3/Kernel_23/classCGAL_1_1Aff__transformation__3.html for definition
      /// of transformation matrix and translation
      void transform(const TransformMatrix& trafo_mat, const PointType& translation, double scale = double(1));

      /// Returns the size in bytes, //TODO: Implement me
      std::size_t bytes() const;


    private:
      /// Delete tree, which also requires to delete the inside tester
      void _delete_tree();
      /// initializes tree and inside tester with already initialized polyhedron
      void _init_wrapper();

    }; // class CGALWrapper
  } // namespace Geometry
} // namespace FEAT
#endif //defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)
#endif //KERNEL_GEOMETRY_CGAL_HPP

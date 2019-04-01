// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_CGAL_HPP
#define KERNEL_GEOMETRY_CGAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

#if defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)

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
      void * _cgal_data;

      /// read in stream in off file format and preprocess search tree for in/out test
      void _parse_off_data(std::istream & file);

      public:
      /// Create a new CGALWrapper Instance and open the provided off file.
      CGALWrapper(String filename);

      /// Create a new CGALWrapper Instance and open the provided off (file-) stream.
      CGALWrapper(std::istream & file);

      /// Check whether a point is inside the Polyhedron defined at objects' construction.
      bool point_inside(double x, double y, double z) const;

      /// Returns the minimun squared distance between the query point and all input primitives defined at objects' construction.
      double squared_distance(double x, double y, double z) const;

      /// Destructor
      ~CGALWrapper();
    };

  } // namespace Geometry
} // namespace FEAT
#endif //defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)
#endif //KERNEL_GEOMETRY_CGAL_HPP

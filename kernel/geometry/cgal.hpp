// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
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
    }; // class CGALWrapper
  } // namespace Geometry
} // namespace FEAT
#endif //defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)
#endif //KERNEL_GEOMETRY_CGAL_HPP

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_FACE_REF_TRAFO_HPP
#define KERNEL_GEOMETRY_INTERN_FACE_REF_TRAFO_HPP 1

#include <kernel/shape.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Face-Reference trafo class template.
       *
       * \todo detailed description
       *
       * \author Peter Zajac
       */
      template<typename Shape_, int face_dim_>
#ifndef DOXYGEN
      class FaceRefTrafo;
#else
      class FaceRefTrafo
      {
      public:
        /**
         * \brief Computes the transformation matrix for an local face.
         *
         * \param[out] a
         * The transformation matrix for the orientation code.
         *
         * \param[out] b
         * The transformation vector for the orientation code.
         *
         * \param[in] face
         * The index of the local face.
         */
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, m, n, sm_, sn_>& a, Tiny::Vector<DataType_, n, sl_>& b, int face);
      };
#endif // DOXYGEN

      /**
       * \brief Face reference trafo implementation for Hypercube<2> shape and 1-face
       */
      template<>
      class FaceRefTrafo<Shape::Hypercube<2>, 1>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 2, 1, sm_, sn_>& a, Tiny::Vector<DataType_, 2, sl_>& b, int face)
        {
          switch(face)
          {
          case 0:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) = -DataType_(1);
            break;

          case 1:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) = +DataType_(1);
            break;

          case 2:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            b(Index(0)) = -DataType_(1);
            b(Index(1)) =  DataType_(0);
            break;

          case 3:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            b(Index(0)) = +DataType_(1);
            b(Index(1)) =  DataType_(0);
            break;
          }
        }
      };

      /**
       * \brief Face reference trafo implementation for Hypercube<3> shape and 1-face
       */
      template<>
      class FaceRefTrafo<Shape::Hypercube<3>, 1>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 3, 1, sm_, sn_>& a, Tiny::Vector<DataType_, 3, sl_>& b, int face)
        {
          switch(face)
          {
          case 0:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) = -DataType_(1);
            b(Index(2)) = -DataType_(1);
            break;

          case 1:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(1);
            b(Index(2)) = -DataType_(1);
            break;

          case 2:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) = -DataType_(1);
            b(Index(2)) =  DataType_(1);
            break;

          case 3:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(1);
            b(Index(2)) =  DataType_(1);
            break;

          case 4:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) = -DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) = -DataType_(1);
            break;

          case 5:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) = -DataType_(1);
            break;

          case 6:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) = -DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(1);
            break;

          case 7:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(1);
            break;

          case 8:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(1);
            b(Index(0)) = -DataType_(1);
            b(Index(1)) = -DataType_(1);
            b(Index(2)) =  DataType_(0);
            break;

          case 9:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(1);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) = -DataType_(1);
            b(Index(2)) =  DataType_(0);
            break;

          case 10:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(1);
            b(Index(0)) = -DataType_(1);
            b(Index(1)) =  DataType_(1);
            b(Index(2)) =  DataType_(0);
            break;

          case 11:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(1);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(1);
            b(Index(2)) =  DataType_(0);
            break;
          }
        }
      };

      /**
       * \brief Face reference trafo implementation for Hypercube<3> shape and 2-face
       */
      template<>
      class FaceRefTrafo<Shape::Hypercube<3>, 2>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 3, 2, sm_, sn_>& a, Tiny::Vector<DataType_, 3, sl_>& b, int face)
        {
          switch(face)
          {
          case 0:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(1);
            a(Index(2),Index(1)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) = -DataType_(1);
            break;

          case 1:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(1);
            a(Index(2),Index(1)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(1);
            break;

          case 2:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(0);
            a(Index(2),Index(1)) = DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) = -DataType_(1);
            b(Index(2)) =  DataType_(0);
            break;

          case 3:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(0);
            a(Index(2),Index(1)) = DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(1);
            b(Index(2)) =  DataType_(0);
            break;

          case 4:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(0);
            a(Index(2),Index(1)) = DataType_(1);
            b(Index(0)) = -DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 5:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(0);
            a(Index(2),Index(1)) = DataType_(1);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;
          }
        }
      };

      /**
       * \brief Face reference trafo implementation for Simplex<2> shape and 1-face
       */
      template<>
      class FaceRefTrafo<Shape::Simplex<2>, 1>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 2, 1, sm_, sn_>& a, Tiny::Vector<DataType_, 2, sl_>& b, int face)
        {
          switch(face)
          {
          case 0:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) =  DataType_(1);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(0);
            break;

          case 1:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = -DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(1);
            break;

          case 2:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            break;
          }
        }
      };

      /**
       * \brief Face reference trafo implementation for Simplex<3> shape and 1-face
       */
      template<>
      class FaceRefTrafo<Shape::Simplex<3>, 1>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 3, 1, sm_, sn_>& a, Tiny::Vector<DataType_, 3, sl_>& b, int face)
        {
          switch(face)
          {
          case 0:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 1:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 2:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 3:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) =  DataType_(1);
            a(Index(2),Index(0)) =  DataType_(0);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 4:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) =  DataType_(0);
            a(Index(2),Index(0)) =  DataType_(1);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 5:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = -DataType_(1);
            a(Index(2),Index(0)) =  DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(1);
            b(Index(2)) =  DataType_(0);
            break;
          }
        }
      };

      /**
       * \brief Face reference trafo implementation for Simplex<3> shape and 2-face
       */
      template<>
      class FaceRefTrafo<Shape::Simplex<3>, 2>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 3, 2, sm_, sn_>& a, Tiny::Vector<DataType_, 3, sl_>& b, int face)
        {
          switch(face)
          {
          case 0:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) =  DataType_(1);
            a(Index(2),Index(0)) =  DataType_(0);
            a(Index(0),Index(1)) = -DataType_(1);
            a(Index(1),Index(1)) =  DataType_(0);
            a(Index(2),Index(1)) =  DataType_(1);
            b(Index(0)) =  DataType_(1);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 1:
            a(Index(0),Index(0)) = DataType_(0);
            a(Index(1),Index(0)) = DataType_(1);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(0);
            a(Index(2),Index(1)) = DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 2:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(0);
            a(Index(2),Index(1)) = DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;

          case 3:
            a(Index(0),Index(0)) = DataType_(1);
            a(Index(1),Index(0)) = DataType_(0);
            a(Index(2),Index(0)) = DataType_(0);
            a(Index(0),Index(1)) = DataType_(0);
            a(Index(1),Index(1)) = DataType_(1);
            a(Index(2),Index(1)) = DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            b(Index(2)) =  DataType_(0);
            break;
          }
        }
      };

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_FACE_REF_TRAFO_HPP

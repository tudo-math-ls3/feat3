// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTER_SHAPE_CONVERT_TRAITS_HPP
#define KERNEL_GEOMETRY_INTER_SHAPE_CONVERT_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct OtherShape;

      template<>
      struct OtherShape<Shape::Vertex>
      {
        typedef Shape::Vertex Type;
      };

      template<int dim_>
      struct OtherShape<Shape::Simplex<dim_> >
      {
        typedef Shape::Hypercube<dim_> Type;
      };

      template<int dim_>
      struct OtherShape<Shape::Hypercube<dim_> >
      {
        typedef Shape::Simplex<dim_> Type;
      };

      template<
        typename Shape_,
        int face_dim_>
      struct ShapeConvertTraits;

      /**
       * \brief ShapeConvertTraits implementation for Vertex shape.
       *
       * \author Peter Zajac
       */
      template<int face_dim_>
      struct ShapeConvertTraits<Shape::Vertex, face_dim_>
      {
        // validate dimensions
        static_assert(face_dim_ == 0, "invalid face dimension");

        /// Shape type
        typedef Shape::Vertex ShapeType;

        /// cell dimension (=0)
        static constexpr int cell_dim = ShapeType::dimension; // = 0

        /// face dimension (=0)
        static constexpr int face_dim = face_dim_; // = 0

        /// Number of faces generated upon refinement
        static constexpr int count = 1;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Vertex," + stringify(face_dim_) + ">";
        }
      }; // struct ShapeConvertTraits<...,Vertex,...>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<1> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<1>, 0>
      {
        /// Shape type
        typedef Shape::Hypercube<1> ShapeType;
        /// Face type
        typedef Shape::Vertex FaceType;
        /// cell dimension
        static constexpr int cell_dim = 1;
        /// face dimension
        static constexpr int face_dim = 0;
        // number of faces generated
        static constexpr int count = 0;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<1>,0>";
        }
      }; // struct ShapeConvertTraits<Hypercube<1>,0>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<1> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<1>, 1>
      {
        /// Shape type
        typedef Shape::Hypercube<1> ShapeType;
        /// Face type
        typedef Shape::Simplex<1> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 1;
        /// face dimension
        static constexpr int face_dim = 1;
        // number of faces generated
        static constexpr int count = 1;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<1>,1>";
        }
      }; // struct ShapeConvertTraits<Hypercube<1>,1>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<2> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<2>, 0>
      {
        /// Shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// Face type
        typedef Shape::Vertex FaceType;
        /// cell dimension
        static constexpr int cell_dim = 2;
        /// face dimension
        static constexpr int face_dim = 0;
        // number of faces generated
        static constexpr int count = 1;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<2>,0>";
        }
      }; // struct ShapeConvertTraits<Hypercube<2>,0>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<2> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<2>, 1>
      {
        /// Shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// Face type
        typedef Shape::Simplex<1> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 2;
        /// face dimension
        static constexpr int face_dim = 1;
        // number of faces generated
        static constexpr int count = 4;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<2>,1>";
        }
      }; // struct ShapeConvertTraits<Hypercube<2>,1>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<2> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<2>, 2>
      {
        /// Shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// Face type
        typedef Shape::Simplex<2> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 2;
        /// face dimension
        static constexpr int face_dim = 2;
        // number of faces generated
        static constexpr int count = 4;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<2>,2>";
        }
      }; // struct ShapeConvertTraits<Hypercube<2>,2>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<3> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<3>, 0>
      {
        /// Shape type
        typedef Shape::Hypercube<3> ShapeType;
        /// Face type
        typedef Shape::Vertex FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 0;
        // number of faces generated
        static constexpr int count = 1;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<3>,0>";
        }
      }; // struct ShapeConvertTraits<Hypercube<3>,0>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<3> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<3>, 1>
      {
        /// Shape type
        typedef Shape::Hypercube<3> ShapeType;
        /// Face type
        typedef Shape::Simplex<1> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 1;
        // number of faces generated
        static constexpr int count = 14;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<3>,1>";
        }
      }; // struct ShapeConvertTraits<Hypercube<3>,1>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<3> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<3>, 2>
      {
        /// Shape type
        typedef Shape::Hypercube<3> ShapeType;
        /// Face type
        typedef Shape::Simplex<2> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 2;
        // number of faces generated
        static constexpr int count = 36;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<3>,2>";
        }
      }; // struct ShapeConvertTraits<Hypercube<3>,2>

      /**
       * \brief ShapeConvertTraits implementation for Hypercube<3> shape.
       *
       * \author Peter Zajac
       */
      template<>
      struct ShapeConvertTraits<Shape::Hypercube<3>, 3>
      {
        /// Shape type
        typedef Shape::Hypercube<3> ShapeType;
        /// Face type
        typedef Shape::Simplex<3> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 3;
        // number of faces generated
        static constexpr int count = 24;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Hypercube<3>,3>";
        }
      }; // struct ShapeConvertTraits<Hypercube<3>,3>

      /* ******************************************************************** */
      /* ******************************************************************** */
      /* ******************************************************************** */

      /**
       * \brief ShapeConvertTraits implementation for Simplex<1> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<1>, 0>
      {
        /// Shape type
        typedef Shape::Simplex<1> ShapeType;
        /// Face type
        typedef Shape::Vertex FaceType;
        /// cell dimension
        static constexpr int cell_dim = 1;
        /// face dimension
        static constexpr int face_dim = 0;
        // number of faces generated
        static constexpr int count = 1;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<1>,0>";
        }
      }; // struct ShapeConvertTraits<Simplex<1>,0>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<1> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<1>, 1>
      {
        /// Shape type
        typedef Shape::Simplex<1> ShapeType;
        /// Face type
        typedef Shape::Hypercube<1> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 1;
        /// face dimension
        static constexpr int face_dim = 1;
        // number of faces generated
        static constexpr int count = 2;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<1>,1>";
        }
      }; // struct ShapeConvertTraits<Simplex<1>,1>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<2> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<2>, 0>
      {
        /// Shape type
        typedef Shape::Simplex<2> ShapeType;
        /// Face type
        typedef Shape::Vertex FaceType;
        /// cell dimension
        static constexpr int cell_dim = 2;
        /// face dimension
        static constexpr int face_dim = 0;
        // number of faces generated
        static constexpr int count = 1;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<2>,0>";
        }
      }; // struct ShapeConvertTraits<Simplex<2>,0>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<2> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<2>, 1>
      {
        /// Shape type
        typedef Shape::Simplex<2> ShapeType;
        /// Face type
        typedef Shape::Hypercube<1> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 2;
        /// face dimension
        static constexpr int face_dim = 1;
        // number of faces generated
        static constexpr int count = 3;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<2>,1>";
        }
      }; // struct ShapeConvertTraits<Simplex<2>,1>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<2> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<2>, 2>
      {
        /// Shape type
        typedef Shape::Simplex<2> ShapeType;
        /// Face type
        typedef Shape::Hypercube<2> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 2;
        /// face dimension
        static constexpr int face_dim = 2;
        // number of faces generated
        static constexpr int count = 3;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<2>,2>";
        }
      }; // struct ShapeConvertTraits<Simplex<2>,2>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<3> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<3>, 0>
      {
        /// Shape type
        typedef Shape::Simplex<3> ShapeType;
        /// Face type
        typedef Shape::Vertex FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 0;
        // number of faces generated
        static constexpr int count = 1;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<3>,0>";
        }
      }; // struct ShapeConvertTraits<Simplex<3>,0>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<3> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<3>, 1>
      {
        /// Shape type
        typedef Shape::Simplex<3> ShapeType;
        /// Face type
        typedef Shape::Hypercube<1> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 1;
        // number of faces generated
        static constexpr int count = 4;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<3>,1>";
        }
      }; // struct ShapeConvertTraits<Simplex<3>,1>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<3> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<3>, 2>
      {
        /// Shape type
        typedef Shape::Simplex<3> ShapeType;
        /// Face type
        typedef Shape::Hypercube<2> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 2;
        // number of faces generated
        static constexpr int count = 6;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<3>,2>";
        }
      }; // struct ShapeConvertTraits<Simplex<3>,2>

      /**
       * \brief ShapeConvertTraits implementation for Simplex<3> shape.
       *
       * \author Stefan Wahlers, Constantin Christof
       */
      template<>
      struct ShapeConvertTraits<Shape::Simplex<3>, 3>
      {
        /// Shape type
        typedef Shape::Simplex<3> ShapeType;
        /// Face type
        typedef Shape::Hypercube<3> FaceType;
        /// cell dimension
        static constexpr int cell_dim = 3;
        /// face dimension
        static constexpr int face_dim = 3;
        // number of faces generated
        static constexpr int count = 4;

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<3>,3>";
        }
      }; // struct ShapeConvertTraits<Simplex<3>,3>

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTER_SHAPE_CONVERT_TRAITS_HPP

#pragma once
#ifndef KERNEL_GEOMETRY_INTER_SHAPE_CONVERT_TRAITS_HPP
#define KERNEL_GEOMETRY_INTER_SHAPE_CONVERT_TRAITS_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
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

        /// dummy enumeration
        enum
        {
          /// cell dimension (=0)
          cell_dim = ShapeType::dimension, // = 0

          /// face dimension (=0)
          face_dim = face_dim_, // = 0

          /// Number of faces generated upon refinement
          count = 1
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 1,
          /// face dimension
          face_dim = 0,
          // number of faces generated
          count = 0
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 1,
          /// face dimension
          face_dim = 1,
          // number of faces generated
          count = 1
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 2,
          /// face dimension
          face_dim = 0,
          // number of faces generated
          count = 1
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 2,
          /// face dimension
          face_dim = 1,
          // number of faces generated
          count = 4
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 2,
          /// face dimension
          face_dim = 2,
          // number of faces generated
          count = 4
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 0,
          // number of faces generated
          count = 1
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 1,
          // number of faces generated
          count = 14
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 2,
          // number of faces generated
          count = 36
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 3,
          // number of faces generated
          count = 24
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 1,
          /// face dimension
          face_dim = 0,
          // number of faces generated
          count = 1
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 1,
          /// face dimension
          face_dim = 1,
          // number of faces generated
          count = 2
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 2,
          /// face dimension
          face_dim = 0,
          // number of faces generated
          count = 1
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 2,
          /// face dimension
          face_dim = 1,
          // number of faces generated
          count = 3
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 2,
          /// face dimension
          face_dim = 2,
          // number of faces generated
          count = 3
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 0,
          // number of faces generated
          count = 1
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 1,
          // number of faces generated
          count = 4
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 2,
          // number of faces generated
          count = 6
        };

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
        /// dummy enumeration
        enum
        {
          /// cell dimension
          cell_dim = 3,
          /// face dimension
          face_dim = 3,
          // number of faces generated
          count = 4
        };

        /// \brief Returns the name of the class.
        static inline String name()
        {
          return "ShapeConvertTraits<Simplex<3>,3>";
        }
      }; // struct ShapeConvertTraits<Simplex<3>,3>

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INTER_SHAPE_CONVERT_TRAITS_HPP

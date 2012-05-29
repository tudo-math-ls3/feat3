#pragma once
#ifndef KERNEL_GEOMETRY_SHAPE_HPP
#define KERNEL_GEOMETRY_SHAPE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/factorial.hpp>

namespace FEAST
{
  /**
   * \brief Geometry namespace
   */
  namespace Geometry
  {
    /**
     * \brief Shape namespace
     */
    namespace Shape
    {
      /**
       * \brief Vertex shape tag struct
       *
       * \author Peter Zajac
       */
      struct Vertex
      {
        /// dummy enumeration
        enum
        {
          /// Vertex dimension
          dimension = 0
        };

        /// Returns the name of the class as a String.
        static String name()
        {
          return "Vertex";
        }
      }; // struct Vertex

      /**
       * \brief Simplex shape tag struct template
       *
       * \author Peter Zajac
       */
      template<int dimension_>
      struct Simplex
      {
        static_assert(dimension_ > 0, "parameter dimension_ must be greater than 0");

        /// dummy enumeration
        enum
        {
          /// Simplex dimension
          dimension = dimension_,
        };

        /// Returns the name of the class as a String.
        static String name()
        {
          return "Simplex<" + stringify(dimension_) + ">";
        }
      }; // class Simplex

      /**
       * \brief Hypercube shape tag struct template
       *
       * \author Peter Zajac
       */
      template<int dimension_>
      struct Hypercube
      {
        static_assert(dimension_ > 0, "parameter dimension_ must be greater than 0");

        /// dummy enumeration
        enum
        {
          /// Hypercube dimension
          dimension = dimension_,
        };

        /// Returns the name of the class as a String.
        static String name()
        {
          return "Hypercube<" + stringify(dimension_) + ">";
        }
      }; // struct Hypercube

      /// 2-Simplex: Triangle
      typedef Simplex<2> Triangle;
      /// 3-Simplex: Tetrahedron
      typedef Simplex<3> Tetrahedron;

      /// 2-Hypercube: Quadrilateral
      typedef Hypercube<2> Quadrilateral;
      /// 3-Hypercube: Hexahedron
      typedef Hypercube<3> Hexahedron;

      /**
       * \brief Face traits tag struct template
       *
       * \tparam Shape_
       * A shape tag class whose face information is to be determined.
       *
       * \tparam face_dim_
       * The dimension of the faces whose information is to be determined.
       * Must be 0 <= \p face_dim_ <= \p Shape_::dimension.
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        int face_dim_>
#ifndef DOXYGEN
      struct FaceTraits;
#else
      struct FaceTraits
      {
        /// Shape type of the face
        typedef ... ShapeType;

        /// dummy enum
        enum
        {
          /// Number of faces of dimension \p face_dim_
          count = ...
        };
      };
#endif // DOXYGEN

      /// \cond internal
      /**
       * \brief explicit FaceTraits specialisation for Vertex shape
       */
      template<>
      struct FaceTraits<Vertex, 0>
      {
        /// Shape type of face
        typedef Vertex ShapeType;

        /// dummy enumeration
        enum
        {
          /// Number of faces per cell
          count = 1
        };
      }; // struct FaceTraits<Vertex,0>

      /**
       * \brief partial FaceTraits specialisation for Simplex shape
       *
       * \author Peter Zajac
       */
      template<
        int cell_dim_,
        int face_dim_>
      struct FaceTraits<Simplex<cell_dim_>, face_dim_>
      {
        static_assert(face_dim_ > 0, "parameter face_dim_ must be greater than 0");
        static_assert(cell_dim_ >= face_dim_, "face_dim_ must not be greater than cell_dim_");

        /// Shape type of face
        typedef Simplex<face_dim_> ShapeType;

        /// dummy enumeration
        enum
        {
          /**
           * \brief Number of faces per cell
           *
           * For an <em>n</em>-Simplex, the number of <em>m</em>-faces is given by
           * \f[ {n+1\choose m+1} \f]
           */
          count = Binomial<cell_dim_ + 1, face_dim_ + 1>::value
        };
      }; // struct FaceTraits<Simplex<...>, ...>

      /**
       * \brief partial FaceTraits specialisation for Simplex shape and Vertex faces
       *
       * \author Peter Zajac
       */
      template<int cell_dim_>
      struct FaceTraits<Simplex<cell_dim_>, 0>
      {
        static_assert(cell_dim_ > 0, "parameter cell_dim_ must be greater than 0");

        /// Shape type of vertex
        typedef Vertex ShapeType;

        /// dummy enumeration
        enum
        {
          /** \brief Number of vertices per cell */
          count = cell_dim_ + 1
        };
      }; // struct FaceTraits<Simplex<...>, 0>

      /**
       * \brief FaceTraits specialisation for Hypercube shape
       *
       * \author Peter Zajac
       */
      template<
        int cell_dim_,
        int face_dim_>
      struct FaceTraits<Hypercube<cell_dim_>, face_dim_>
      {
        static_assert(face_dim_ > 0, "parameter face_dim_ must be greater than 0");
        static_assert(cell_dim_ >= face_dim_, "face_dim_ must not be greater than cell_dim_");

        ///  Shape type of face
        typedef Hypercube<face_dim_> ShapeType;

        /// dummy enumeration
        enum
        {
          /**
           * \brief Number of faces per cell
           *
           * For an <em>n</em>-Hypercube, the number of <em>m</em>-faces is given by
           * \f[ 2^{(n-m)}\cdot {n\choose m} \f]
           */
          count = (1 << (cell_dim_ - face_dim_)) * Binomial<cell_dim_, face_dim_>::value
        };
      }; // struct FaceTraits<Hypercube<...>, ...>

      /**
       * \brief partial FaceTraits specialisation for Hypercube shape and Vertex faces
       *
       * \author Peter Zajac
       */
      template<int cell_dim_>
      struct FaceTraits<Hypercube<cell_dim_>, 0>
      {
        static_assert(cell_dim_ > 0, "parameter cell_dim_ must be greater than 0");

        /// Shape type of vertex
        typedef Vertex ShapeType;

        /// dummy enumeration
        enum
        {
          /** \brief Number of vertices per cell */
          count = (1 << cell_dim_)
        };
      }; // struct FaceTraits<HyperCube<...>, 0>
      /// \endcond
    } // namespace Shape
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_SHAPE_HPP

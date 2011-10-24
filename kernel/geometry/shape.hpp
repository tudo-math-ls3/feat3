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

          /// number of vertices per cell
          num_verts = dimension_+1
        };

        /**
         * \brief Returns the number of faces.
         *
         * \param[in] face_dim
         * The dimension of the faces whose count is to be returned.
         *
         * \returns
         * The number of faces of dimension \p face_dim per simplex.
         */
        static inline int num_faces(int face_dim)
        {
          CONTEXT("Simplex::num_faces()");
          ASSERT(face_dim <= dimension, "invalid face_dim parameter");
          ASSERT(face_dim >= 0, "invalid face_dim parameter");
          return int(binomial(dimension+1, face_dim+1));
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

          /// number of vertices per cell
          num_verts = (1 << dimension_) // = 2^n
        };

        /**
         * \brief Returns the number of faces.
         *
         * \param[in] face_dim
         * The dimension of the faces whose count is to be returned.
         *
         * \returns
         * The number of faces of dimension \p face_dim per hypercube.
         */
        static inline int num_faces(int face_dim)
        {
          CONTEXT("Hypercube::num_faces()");
          ASSERT(face_dim <= dimension, "invalid face_dim parameter");
          ASSERT(face_dim >= 0, "invalid face_dim parameter");
          return (1 << (dimension - face_dim)) * int(binomial(dimension, face_dim));
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
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        int face_dim_>
      struct FaceTraits;

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
          count = Binomial<cell_dim_+1, face_dim_+1>::value
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
          count = Simplex<cell_dim_>::num_verts
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
          count = (1 << (cell_dim_ - face_dim_)) * Binomial<cell_dim_,face_dim_>::value
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
          count = Hypercube<cell_dim_>::num_verts
        };
      }; // struct FaceTraits<HyperCube<...>, 0>
    } // namespace Shape
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_SHAPE_HPP

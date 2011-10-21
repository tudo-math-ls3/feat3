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
     * \brief Simplex shape class template
     *
     * \author Peter Zajac
     */
    template<int dimension_>
    struct Simplex
    {
      static_assert(dimension_ > 0, "parameter dimension_ must be greater than 0");
      enum
      {
        /// Simplex dimension
        dimension = dimension_,

        /// number of vertices per cell
        num_verts = dimension_+1
      };
    }; // class Simplex

    /**
     * \brief Hypercube shape class template
     *
     * \author Peter Zajac
     */
    template<int dimension_>
    struct Hypercube
    {
      static_assert(dimension_ > 0, "parameter dimension_ must be greater than 0");
      enum
      {
        /// Hypercube dimension
        dimension = dimension_,

        /// number of vertices per cell
        num_verts = (1 << dimension_) // = 2^n
      };
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
     * \brief Face traits class template
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      int face_dim_>
    struct FaceTraits;

    /**
     * \brief FaceTraits specialisation for Simplex shape
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
    }; // struct FaceTraits<Simplex>

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
    }; // struct FaceTraits<Hypercube>

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_SHAPE_HPP

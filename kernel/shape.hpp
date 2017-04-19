#pragma once
#ifndef KERNEL_SHAPE_HPP
#define KERNEL_SHAPE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/meta_math.hpp>
#include <kernel/util/string.hpp>

namespace FEAT
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
      /// Vertex dimension
      static constexpr int dimension = 0;

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

      /// Simplex dimension
      static constexpr int dimension = dimension_;

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

      /// Hypercube dimension
      static constexpr int dimension = dimension_;

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

      /// Number of faces of dimension \p face_dim_
      static constexpr int count = ...
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

      /// Number of faces per cell
      static constexpr int count = 1;
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

      /**
       * \brief Number of faces per cell
       *
       * For an <em>n</em>-Simplex, the number of <em>m</em>-faces is given by
       * \f[ {n+1\choose m+1} \f]
       */
      static constexpr int count = MetaMath::Binomial<cell_dim_ + 1, face_dim_ + 1>::value;
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

      /// \brief Number of vertices per cell
      static constexpr int count = cell_dim_ + 1;
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

      /**
       * \brief Number of faces per cell
       *
       * For an <em>n</em>-Hypercube, the number of <em>m</em>-faces is given by
       * \f[ 2^{(n-m)}\cdot {n\choose m} \f]
       */
      static constexpr int count = (1 << (cell_dim_ - face_dim_)) * MetaMath::Binomial<cell_dim_, face_dim_>::value;
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

      /// Number of vertices per cell
      static constexpr int count = (1 << cell_dim_);
    }; // struct FaceTraits<HyperCube<...>, 0>
    /// \endcond

    /**
     * \brief Reference cell traits structure.
     *
     * This class template contains information about the reference cell of the corresponding shape,
     * especially the coordinates of the reference cell's vertices.
     *
     * \tparam Shape_
     * The shape whose reference cell is to be determined.
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
#ifndef DOXYGEN
    struct ReferenceCell;
#else
    struct ReferenceCell
    {
      /**
       * \brief Returns the coordinate of a reference cell vertex.
       *
       * \tparam T_
       * The desired type of the coordinate.
       *
       * \param[in] vertex_idx
       * The index of the vertex whose coordinate is to be returned.
       *
       * \param[in] coord_idx
       * The index of the coordinate of the vertex that is to be returned.
       *
       * \returns
       * The desired coordinate of the reference cell vertex.
       */
      template<typename T_>
      static T_ vertex(int vertex_idx, int coord_idx);

      /**
       * \brief Returns the coordinate of the reference cell centre.
       *
       * \tparam T_
       * The desired type of the coordinate.
       *
       * \param[in] coord_idx
       * The index of the coordinate of the centre that is to be returned.
       *
       * \returns
       * The desired coordinate of the reference cell centre.
       */
      template<typename T_>
      static T_ centre(int coord_idx);

      /**
       * \brief Returns the volume of the reference cell.
       *
       * \tparam T_
       * The desired type of the volume.
       *
       * \returns
       * The volume of the reference cell.
       */
      template<typename T_>
      static T_ volume();

      /**
       * \brief Returns the orientation of a facet.
       *
       * \note For simplices, all facets are positively oriented (i.e. the edges that make up the reference
       * triangle), but for hypercubes this is not the case.
       *
       * \param[in] facet_index
       * The number of the facet to return the orientation for.
       *
       * \returns
       * 1 for a positively or -1 for a negatively oriented facet.
       */
      static int orientation(int facet_index);
    };
#endif // DOXYGEN

    /// \cond internal
    /**
     * \brief ReferenceCell specialisation for Vertex shape
     *
     * \author Peter Zajac
     */
    template<>
    struct ReferenceCell<Shape::Vertex>
    {
      template<typename T_>
      static T_ vertex(int, int)
      {
        return T_(0);
      }

      template<typename T_>
      static T_ centre(int)
      {
        return T_(0);
      }

      template<typename T_>
      static T_ volume()
      {
        return T_(0);
      }

      // orientation is not implemented due to being nonsensical

    };

    /**
     * \brief Partial ReferenceCell specialisation for Simplex shape
     *
     * \author Peter Zajac
     */
    template<int dim_>
    struct ReferenceCell< Shape::Simplex<dim_> >
    {
      template<typename T_>
      static T_ vertex(int vertex_idx, int coord_idx)
      {
        return (coord_idx + 1) == vertex_idx ? T_(1) : T_(0);
      }

      template<typename T_>
      static T_ centre(int)
      {
        return T_(1) / T_(dim_+1);
      }

      template<typename T_>
      static T_ volume()
      {
        return T_(1) / T_(MetaMath::Factorial<dim_>::value);
      }

      /**
       * \brief Facet orientation for simplices
       *
       * For Simplex<1>, the left vertex is negatively oriented (the outer normal being -1), the right vertex is
       * positively oriented (the outer normal being 1).
       * Simplex<2> and Simplex<3> only have positively oriented facets.
       *
       * \param[in] facet_index
       * The number of the facet to return the orientation for.
       *
       * \returns
       * 1 for a positively or -1 for a negatively oriented facet.
       *
       */
      static int orientation(int facet_index)
      {
        ASSERTM((facet_index >= 0), "facet index out of range!");
        ASSERTM((facet_index < FaceTraits<Simplex<dim_>,dim_-1>::count), "facet index out of range!");
#ifndef DEBUG
        (void) facet_index;
#endif
        return 1 - ( dim_ > 1 ? 0 : (( facet_index & 1 ) ^ 1 ) << 1);
      }
    };

    /**
     * \brief Partial ReferenceCell specialisation for Hypercube shape
     *
     * \author Peter Zajac
     */
    template<int dim_>
    struct ReferenceCell< Shape::Hypercube<dim_> >
    {
      template<typename T_>
      static T_ vertex(int vertex_idx, int coord_idx)
      {
        return T_((((vertex_idx >> coord_idx) & 1) << 1) - 1);
      }

      template<typename T_>
      static T_ centre(int)
      {
        return T_(0);
      }

      template<typename T_>
      static T_ volume()
      {
        return T_(1 << dim_); // = 2^dim
      }

      /**
       * \brief Facet orientation for Hypercubes
       *
       * For Hypercube<1>, the facets are vertices, left is oriented negatively and right is oriented positively so
       * that the normals are outer normals: -1 +1
       * For Hypercube<2>, edges 0 and 3 are positively oriented: +1 -1 -1 +1
       * For Hypercube<3>, faces 0, 3, 4 are negatively oriented: -1 +1 +1 -1 -1 +1
       *
       * \param[in] facet_index
       * The number of the facet to return the orientation for.
       *
       * \returns
       * 1 for a positively or -1 for a negatively oriented facet.
       *
       */
      static int orientation(int facet_index)
      {
        ASSERTM((facet_index >= 0), "facet index out of range!");
        ASSERTM((facet_index < FaceTraits<Hypercube<dim_>,dim_-1>::count), "facet index out of range!");

        return 1 - (( ((facet_index >> 1) ^ facet_index ^ dim_) & 1 ) << 1);
      }
    };

    /// \endcond
  } // namespace Shape
} // namespace FEAT

#endif // KERNEL_SHAPE_HPP

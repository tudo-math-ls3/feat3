// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_VERTEX_SET_HPP
#define KERNEL_GEOMETRY_VERTEX_SET_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/adjacency/permutation.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Fixed-Sized Vertex Set class template
     *
     * This class implements a vertex set whose number of coordinates is specified at compile-time.
     *
     * \tparam num_coords_
     * The number of coordinates per vertex. Must be > 0.
     *
     * \tparam Coord_
     * The data type for a vertex coordinate.
     *
     * \author Peter Zajac
     */
    template<
      int num_coords_,
      typename Coord_ = Real>
    struct VertexSet
    {
      static_assert(num_coords_ > 0, "invalid num_coords_ parameter");

    public:
      /// number of coordinates per vertex
      static constexpr int num_coords = num_coords_;

      /// vertex coordinate type
      typedef Coord_ CoordType;

      /// vertex type
      typedef Tiny::Vector<CoordType, num_coords> VertexType;

    protected:
      /// vertex vector
      std::vector<VertexType> _vertices;

      /// \cond internal
      template<int sm_, int sn_, int sv_>
      static void _aux_rot_mat(
        Tiny::Matrix<CoordType, 1, 1, sm_, sn_>& r,
        const Tiny::Vector<CoordType, 1, sv_>&)
      {
        r.set_identity();
      }

      template<int sm_, int sn_, int sv_>
      static void _aux_rot_mat(
        Tiny::Matrix<CoordType, 2, 2, sm_, sn_>& r,
        const Tiny::Vector<CoordType, 2, sv_>& a)
      {
        r.set_rotation_2d(a[0]);
      }

      template<int sm_, int sn_, int sv_>
      static void _aux_rot_mat(
        Tiny::Matrix<CoordType, 3, 3, sm_, sn_>& r,
        const Tiny::Vector<CoordType, 3, sv_>& a)
      {
        r.set_rotation_3d(a[0], a[1], a[2]);
      }

      /// \endcond

      /// internal clone constructor
      explicit VertexSet(const std::vector<VertexType>& vtx) :
        _vertices(vtx)
      {
      }

    public:
      /// standard constructor
      VertexSet() :
        _vertices()
      {
      }

      /**
       * \brief Constructor.
       *
       * \param[in] num_vertices
       * The number of vertices to be allocated in the vertex set.
       */
      explicit VertexSet(Index num_vertices) :
        _vertices(std::size_t(num_vertices))
      {
      }

      /// move constructor
      VertexSet(VertexSet&& other) :
        _vertices(std::forward<std::vector<VertexType>>(other._vertices))
      {
      }

      /// move-assignment operator
      VertexSet& operator=(VertexSet&& other)
      {
        // avoid self-move
        if(this == &other)
          return *this;
        _vertices = std::forward<std::vector<VertexType>>(other._vertices);
        return *this;
      }

      /// virtual destructor
      virtual ~VertexSet()
      {
      }

      /// \returns An independent clone of this vertex set object.
      VertexSet clone() const
      {
        return VertexSet(this->_vertices);
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return _vertices.size() * sizeof(VertexType);
      }

      /// Returns the number of coordinates per vertex.
      int get_num_coords() const
      {
        return num_coords_;
      }

      /// Returns the number of vertices in the vertex set.
      Index get_num_vertices() const
      {
        return Index(_vertices.size());
      }

      /**
       * \brief Returns a reference to a vertex.
       *
       * \param[in] i
       * The index of the vertex to be returned.
       *
       * \returns
       * A reference to the vertex.
       */
      VertexType& operator[](Index i)
      {
        ASSERT(i < get_num_vertices());
        return _vertices[i];
      }

      /** \copydoc operator[]() */
      const VertexType& operator[](Index i) const
      {
        ASSERT(i < get_num_vertices());
        return _vertices[i];
      }

      /**
       * \brief Applies a "proper rigid" transformation onto the vertex set.
       *
       * Let \e v denote the \p origin world point, \e w the \p offset world point and \e R
       * the rotation matrix corresponding to the \p angles, then this function applies the
       * following transformation for any vertex \e x of the vertex set:
       *
       *   \f[ x \mapsto w + R\cdot (x - v) \f]
       *
       * \param[in] origin
       * The origin of the transformation. This is subtracted from any vertex before applying the
       * rotation.
       *
       * \param[in] angles
       * The angles of the rotation matrix.
       * - 2D: the rotation angle in radians is stored as:
       *   - angles(0): rotation angle
       *   - angles(1): \e ignored
       * - 3D: the rotation angles in radians stored as:
       *   - angles(0): yaw angle
       *   - angles(1): pitch angle
       *   - angles(2): roll angle
       *
       * \param[in] offset
       * The offset of the transformation. This is added to any vertex after applying the rotation.
       */
      void transform(const VertexType& origin, const VertexType& angles, const VertexType& offset)
      {
        // create rotation matrix
        Tiny::Matrix<CoordType, num_coords_, num_coords_> rot;
        _aux_rot_mat(rot, angles);

        // transform all vertices
        for(auto& v : _vertices)
        {
          v.set_mat_vec_mult(rot, v - origin) += offset;
        }
      }

      /**
       * \brief Applies a permutation onto this vertex set.
       *
       * \param[in] perm
       * The permutation that is to be applied.
       *
       * \param[in] invert
       * Specifies whether to apply the inverse permutation instead of the forward one.
       */
      void permute(const Adjacency::Permutation& perm, bool invert = false)
      {
        if(!perm.empty() && !_vertices.empty())
        {
          XASSERT(Index(_vertices.size()) == perm.size());
          perm.apply(_vertices.data(), invert);
        }
      }

      /// Returns the name of the class.
      static String name()
      {
        return "VertexSet<...>";
      }
    }; // class VertexSet<...>

  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_VERTEX_SET_HPP

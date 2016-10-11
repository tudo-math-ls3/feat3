#pragma once
#ifndef KERNEL_GEOMETRY_VERTEX_SET_HPP
#define KERNEL_GEOMETRY_VERTEX_SET_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>
#include <kernel/util/tiny_algebra.hpp>

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
     * \tparam stride_
     * The stride for the vertex. Must be >= num_coords_.
     *
     * \tparam Coord_
     * The data type for a vertex coordinate.
     *
     * \author Peter Zajac
     */
    template<
      int num_coords_,
      int stride_ = num_coords_,
      typename Coord_ = Real>
    struct VertexSet
    {
      static_assert(num_coords_ > 0, "invalid num_coords_ parameter");
      static_assert(stride_ >= num_coords_, "invalid stride_ parameter");

    public:
      /// number of coordinates per vertex
      static constexpr int num_coords = num_coords_;
      /// vertex stride
      static constexpr int stride = stride_;

      /// vertex coordinate type
      typedef Coord_ CoordType;

      /// vertex type
      //typedef CoordType VertexType[stride_];
      typedef Tiny::Vector<CoordType, num_coords, stride> VertexType;

      /// vertex reference type
      typedef VertexType& VertexReference;

      /// const vertex reference type
      typedef const VertexType& ConstVertexReference;

    protected:
      /// number of vertices
      Index _num_vertices;

      /// vertex array
      VertexType* _vertices;

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

    private:
      VertexSet& operator=(const VertexSet&);

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_vertices
       * The number of vertices to be allocated in the vertex set.
       */
      explicit VertexSet(Index num_vertices) :
        _num_vertices(num_vertices),
        _vertices(nullptr)
      {
        if(num_vertices > 0)
        {
          _vertices = new VertexType[num_vertices];
        }
      }

      /// \cond internal
      explicit VertexSet(
        Index num_vertices,
        int /*num_coords*/, // ignored
        int /*stride*/)     // ignored
         :
        _num_vertices(num_vertices),
        _vertices(nullptr)
      {
        if(num_vertices > 0)
        {
          _vertices = new VertexType[num_vertices];
        }
      }
      /// \endcond

      /**
       * \brief Copy Constructor
       *
       * \param[in] other
       * The vertex set that is to be copied.
       */
      //template<int stride2_, typename Coord2_>
      //VertexSet(const VertexSet<num_coords_, stride2_, Coord2_>& other) :
      VertexSet(const VertexSet& other) :
        _num_vertices(other.get_num_vertices()),
        _vertices(nullptr)
      {
        if(_num_vertices > 0)
        {
          _vertices = new VertexType[_num_vertices];
          for(Index i(0); i < _num_vertices; ++i)
          {
            for(int j(0); j < num_coords_; ++j)
            {
              _vertices[i][j] = Coord_(other[i][j]);
            }
          }
        }
      }

      /// virtual destructor
      virtual ~VertexSet()
      {
        if(_vertices != nullptr)
        {
          delete [] _vertices;
        }
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return std::size_t(_num_vertices) * sizeof(VertexType);
      }

      /// Returns the number of coordinates per vertex.
      int get_num_coords() const
      {
        return num_coords_;
      }

      /// Returns the vertex stride.
      int get_stride() const
      {
        return stride_;
      }

      /// Returns the number of vertices in the vertex set.
      Index get_num_vertices() const
      {
        return _num_vertices;
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
      VertexReference operator[](Index i)
      {
        ASSERT(_vertices != nullptr);
        ASSERT(i < _num_vertices);
        return _vertices[i];
      }

      /** \copydoc operator[]() */
      ConstVertexReference operator[](Index i) const
      {
        ASSERT(_vertices != nullptr);
        ASSERT(i < _num_vertices);
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

        // loop over all vertices
        VertexType tmp;
        for(Index i(0); i < _num_vertices; ++i)
        {
          tmp = _vertices[i] - origin;
          _vertices[i].set_mat_vec_mult(rot, tmp) += offset;
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

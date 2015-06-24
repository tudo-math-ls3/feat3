#pragma once
#ifndef KERNEL_GEOMETRY_VERTEX_SET_HPP
#define KERNEL_GEOMETRY_VERTEX_SET_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
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
      typedef CoordType VertexType[stride_];

      /// vertex reference type
      typedef VertexType& VertexReference;

      /// const vertex reference type
      typedef const VertexType& ConstVertexReference;

    protected:
      /// number of vertices
      Index _num_vertices;

      /// vertex array
      VertexType* _vertices;

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
        CONTEXT(name() + "::VertexSet()");
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
        CONTEXT(name() + "::VertexSet()");
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
        CONTEXT(name() + "::~VertexSet()");
        if(_vertices != nullptr)
        {
          delete [] _vertices;
        }
      }

      /// Returns the number of coordinates per vertex.
      int get_num_coords() const
      {
        CONTEXT(name() + "::get_num_coords()");
        return num_coords_;
      }

      /// Returns the vertex stride.
      int get_stride() const
      {
        CONTEXT(name() + "::get_num_stride()");
        return stride_;
      }

      /// Returns the number of vertices in the vertex set.
      Index get_num_vertices() const
      {
        CONTEXT(name() + "::get_num_vertices()");
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
        CONTEXT(name() + "::operator[]()");
        ASSERT_(_vertices != nullptr);
        ASSERT_(i < _num_vertices);
        return _vertices[i];
      }

      /** \copydoc operator[]() */
      ConstVertexReference operator[](Index i) const
      {
        CONTEXT(name() + "::operator[]() [const]");
        ASSERT_(_vertices != nullptr);
        ASSERT_(i < _num_vertices);
        return _vertices[i];
      }

      /// Returns the name of the class.
      static String name()
      {
        return "VertexSet<...>";
      }
    }; // class VertexSet<...>

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_VERTEX_SET_HPP

#pragma once
#ifndef KERNEL_GEOMETRY_VERTEX_SET_HPP
#define KERNEL_GEOMETRY_VERTEX_SET_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
{
  namespace Geometry
  {
    // forward declaration
    template<typename Coord_ = Real>
    class VertexSetVariable;

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
     * \see VertexSetVariable
     *
     * \author Peter Zajac
     */
    template<
      int num_coords_,
      int stride_ = num_coords_,
      typename Coord_ = Real>
    struct VertexSetFixed
    {
      static_assert(num_coords_ > 0, "invalid num_coords_ parameter");
      static_assert(stride_ >= num_coords_, "invalid stride_ parameter");

    public:
      /// dummy enum
      enum
      {
        /// number of coordinates per vertex
        num_coords = num_coords_,
        /// vertex stride
        stride = stride_
      };

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
      VertexSetFixed& operator=(const VertexSetFixed&);

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_vertices
       * The number of vertices to be allocated in the vertex set.
       */
      explicit VertexSetFixed(Index num_vertices) :
        _num_vertices(num_vertices),
        _vertices(nullptr)
      {
        CONTEXT(name() + "::VertexSetFixed()");
        if(num_vertices > 0)
        {
          _vertices = new VertexType[num_vertices];
        }
      }

      /// \cond internal
      explicit VertexSetFixed(
        Index num_vertices,
        int /*num_coords*/, // ignored
        int /*stride*/)     // ignored
         :
        _num_vertices(num_vertices),
        _vertices(nullptr)
      {
        CONTEXT(name() + "::VertexSetFixed()");
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
      template<
        int stride2_,
        typename Coord2_>
      VertexSetFixed(const VertexSetFixed<num_coords_, stride2_, Coord2_>& other) :
        _num_vertices(other.get_num_vertices()),
        _vertices(nullptr)
      {
        if(_num_vertices > 0)
        {
          _vertices = new VertexType(_num_vertices);
          for(Index i(0); i < _num_vertices; ++i)
          {
            for(int j(0); j < num_coords_; ++j)
            {
              _vertices[i][j] = Coord_(other[i][j]);
            }
          }
        }
      }

      /**
       * \brief Conversion Constructor.
       *
       * \param[in] other
       * The variable vertex set that is to be copied.
       */
      template<typename Coord2_>
      explicit VertexSetFixed(const VertexSetVariable<Coord2_>& other) :
        _num_vertices(other.get_num_vertices()),
        _vertices(nullptr)
      {
        ASSERT(other.get_num_coords() == num_coords_, "Coordinate count mismatch!");
        if(_num_vertices > 0)
        {
          _vertices = new VertexType(_num_vertices);
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
      virtual ~VertexSetFixed()
      {
        CONTEXT(name() + "::~VertexSetFixed()");
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
        return "VertexSetFixed<...>";
      }
    }; // class VertexSetFixed<...>

    /**
     * \brief Variable-Sized Vertex-Set class template
     *
     * This class template implements a vertex set whose number of coordinates is specified at runtime.
     *
     * \tparam Coord_
     * The datatype for a vertex coordinate.
     *
     * \see VertexSetFixed
     *
     * \author Peter Zajac
     */
    template<typename Coord_>
    class VertexSetVariable
    {
    public:
      /// coordinate type
      typedef Coord_ CoordType;

      /// vertex reference type
      typedef CoordType* VertexReference;

      /// const vertex reference type
      typedef const CoordType* ConstVertexReference;

    protected:
      /// number of vertices
      Index _num_vertices;

      /// number of coordinates per vertex
      int _num_coords;

      /// vertex stride
      int _stride;

      /// vertex array
      CoordType* _vertices;

    private:
      VertexSetVariable& operator=(const VertexSetVariable&);

    public:
      /**
       * \brief Constructor
       *
       * \param[in] num_vertices
       * The number of vertices to be allocated.
       *
       * \param[in] num_coords
       * The number of coordinates per vertex. Must be >= 0.
       *
       * \param[in] stride
       * The vertex stride. If = 0, then num_coords will be used for the vertex stride.
       */
      explicit VertexSetVariable(
        Index num_vertices,
        int num_coords,
        int stride = 0)
         :
        _num_vertices(num_vertices),
        _num_coords(num_coords),
        _stride(stride > 0 ? stride : num_coords),
        _vertices(nullptr)
      {
        CONTEXT(name() + "::VertexSetVariable()");
        ASSERT_(_num_coords >= 0);
        ASSERT_(_stride >= _num_coords);
        if((_num_vertices > 0) && (_num_coords > 0))
        {
          _vertices = new CoordType[_num_vertices * _stride];
        }
      }

      /**
       * \brief Copy Constructor
       *
       * \param[in] other
       * The vertex set that is to be copied.
       */
      template<typename Coord2_>
      VertexSetVariable(const VertexSetVariable<Coord2_>& other) :
        _num_vertices(other.get_num_vertices()),
        _num_coords(other.get_num_coords()),
        _stride(other.get_stride()),
        _vertices(nullptr)
      {
        if((_num_vertices > 0) && (_num_coords > 0))
        {
          _vertices = new CoordType[_num_vertices * _stride];
          for(Index i(0); i < _num_vertices; ++i)
          {
            for(int j(0); j < _num_coords; ++j)
            {
              _vertices[i * _stride + j] = CoordType(other[i][j]);
            }
            for(int j(_num_coords); j < _stride; ++j)
            {
              _vertices[i * _stride + j] = CoordType(0);
            }
          }
        }
      }

      /**
       * \brief Conversion Constructor
       *
       * \param[in] other
       * The fixed vertex set that is to be copied.
       *
       * \param[in] stride
       * The stride for the vertex set. If set to 0, the stride of the fixed vertex set is used.
       */
      template<int num_coords_, int stride_, typename Coord2_>
      explicit VertexSetVariable(
        const VertexSetFixed<num_coords_, stride_, Coord2_>& other,
        int stride = 0) :
        _num_vertices(other.get_num_vertices()),
        _num_coords(num_coords_),
        _stride(stride == 0 ? stride_ : stride),
        _vertices(nullptr)
      {
        ASSERT_(_stride >= _num_coords);
        if((_num_vertices > 0) && (_num_coords > 0))
        {
          _vertices = new CoordType[_num_vertices * _stride];
          for(Index i(0); i < _num_vertices; ++i)
          {
            for(int j(0); j < _num_coords; ++j)
            {
              _vertices[i * _stride + j] = CoordType(other[i][j]);
            }
            for(int j(_num_coords); j < _stride; ++j)
            {
              _vertices[i * _stride + j] = CoordType(0);
            }
          }
        }
      }

      /// virtual destructor
      virtual ~VertexSetVariable()
      {
        CONTEXT(name() + "::~VertexSetVariable()");
        if(_vertices != nullptr)
        {
          delete [] _vertices;
        }
      }

      /**
       * \brief Returns the number of coordinates per vertex.
       *
       * \warning This function nay return 0.
       */
      int get_num_coords() const
      {
        CONTEXT(name() + "::get_num_coords()");
        return _num_coords;
      }

      /**
       * \brief Returns the vertex stride.
       *
       * \warning This function may return 0.
       */
      int get_stride() const
      {
        CONTEXT(name() + "::get_stride()");
        return _stride;
      }

      /// Returns the number of vertices
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
        return &_vertices[_stride * i];
      }

      /** \copydoc operator[]() */
      ConstVertexReference operator[](Index i) const
      {
        CONTEXT(name() + "::operator[]() [const]");
        ASSERT_(_vertices != nullptr);
        ASSERT_(i < _num_vertices);
        return &_vertices[_stride * i];
      }

      /// Returns the name of the class.
      static String name()
      {
        return "VertexSetVariable<...>";
      }
    }; // class VertexSetVariable<...>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_VERTEX_SET_HPP

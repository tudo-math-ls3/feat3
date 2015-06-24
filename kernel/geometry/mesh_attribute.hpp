#pragma once
#ifndef KERNEL_GEOMETRY_MESH_ATTRIBUTE
#define KERNEL_GEOMETRY_MESH_ATTRIBUTE 1

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Container for saving data related to mesh entities
     *
     * This class template implements a container for saving attributes of vertices/edges/faces/cells and is used in
     * MeshPart. The attributes do not need to be scalar.
     *
     * This used to be VertexSetVariable and the interface is still compatible with VertexSetFixed, even though some
     * variable and interface names are nonsensical now. Through this, classes like the StandardVertexRefiner (which)
     * is used to refine MeshAttributes referring to vertices) can still be applied to both classes.
     *
     * \tparam Coord_
     * The datatype for a vertex coordinate.
     *
     * \see VertexSetFixed
     *
     * \author Peter Zajac
     */
    template<typename Coord_>
    class MeshAttribute
    {
      public:
        /// Type for the data in the attributes
        typedef Coord_ CoordType;
        /// Data reference type
        typedef CoordType* VertexReference;
        /// Const data reference type
        typedef const CoordType* ConstVertexReference;

      protected:
        /// Number of values
        Index _num_vertices;
        /// Number of entries per value
        int _num_coords;
        /// Padded number of entries per value
        int _stride;
        /// Value array
        CoordType* _vertices;
        /// Name identifier
        String _identifier;

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
         *
         * \param[in] identifier
         * Name identifier, defaults to the empty String.
         */
        explicit MeshAttribute(
          Index num_vertices,
          int num_coords,
          int stride = 0,
          String identifier_ = "")
          :
            _num_vertices(num_vertices),
            _num_coords(num_coords),
            _stride(stride > 0 ? stride : num_coords),
            _vertices(nullptr),
            _identifier(identifier_)
            {
              CONTEXT(name() + "::MeshAttribute()");
              ASSERT_(_num_coords >= 0);
              ASSERT_(_stride >= _num_coords);
              if((_num_vertices > 0) && (_num_coords > 0))
              {
                _vertices = new CoordType[_num_vertices * Index(_stride)];

                for(Index i(0); i < _num_vertices*Index(_stride); ++i)
                  _vertices[i] = CoordType(0);

              }
            }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other
         * The vertex set that is to be copied.
         */
        //template<typename Coord2_>
        //MeshAttribute(const MeshAttribute<Coord2_>& other) :
        MeshAttribute(const MeshAttribute& other) :
          _num_vertices(other.get_num_vertices()),
          _num_coords(other.get_num_coords()),
          _stride(other.get_stride()),
          _vertices(nullptr),
          _identifier(other.get_identifier())
          {
            if((_num_vertices > 0) && (_num_coords > 0))
            {
              _vertices = new CoordType[_num_vertices * Index(_stride)];
              for(Index i(0); i < _num_vertices; ++i)
              {
                for(Index j(0); j < Index(_num_coords); ++j)
                {
                  _vertices[i * Index(_stride) + j] = CoordType(other[i][j]);
                }
                for(Index j = Index(_num_coords); j < Index(_stride); ++j)
                {
                  _vertices[i * Index(_stride) + j] = CoordType(0);
                }
              }
            }
          }

        /**
         * \brief Move Constructor
         *
         * \param[in] other
         * The vertex set to be moved
         */
        MeshAttribute(MeshAttribute&& other) :
          _num_vertices(other.get_num_vertices()),
          _num_coords(other.get_num_coords()),
          _stride(other.get_stride()),
          _vertices(nullptr),
          _identifier(other.get_identifier())
          {
            if((_num_vertices > 0) && (_num_coords > 0))
              _vertices = other._vertices;

            other._num_vertices = 0;
            other._num_coords = 0;
            other._stride = 0;
            other._vertices = nullptr;
            other._identifier = "";
          }

        /// virtual destructor
        virtual ~MeshAttribute()
        {
          CONTEXT(name() + "::~MeshAttribute()");
          if(_vertices != nullptr)
          {
            delete [] _vertices;
          }
        }

        /**
         * \brief Copy assignment operator
         *
         * This is needed for replacing one MeshAttribute by another when inserting into a \ref MeshAttributeHolder.
         */
        MeshAttribute& operator=(const MeshAttribute& other)
        {
          _num_vertices = other.get_num_vertices();
          _num_coords = other.get_num_coords();
          _stride = other.get_stride();
          _identifier = other.get_identifier();

          if(_vertices != nullptr)
            delete[] _vertices;

          if(other._vertices != nullptr)
          {
            if((_num_vertices > 0) && (_num_coords > 0))
            {
              _vertices = new CoordType[_num_vertices * Index(_stride)];
              for(Index i(0); i < _num_vertices; ++i)
              {
                for(Index j(0); j < Index(_num_coords); ++j)
                  _vertices[i * Index(_stride) + j] = other[i][j];

                for(Index j = Index(_num_coords); j < Index(_stride); ++j)
                  _vertices[i * Index(_stride) + j] = CoordType(0);
              }
            }
          }
          return *this;
        }

        /**
         * \brief Returns a copy of the name identifier
         */
        String get_identifier() const
        {
          return _identifier;
        }

        /**
         * \brief Returns the number of coordinates per vertex.
         *
         * \warning This function may return 0.
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
          return &_vertices[Index(_stride) * i];
        }

        /** \copydoc operator[]() */
        ConstVertexReference operator[](Index i) const
        {
          CONTEXT(name() + "::operator[]() [const]");
          ASSERT_(_vertices != nullptr);
          ASSERT_(i < _num_vertices);
          return &_vertices[Index(_stride) * i];
        }

        /// Returns the name of the class.
        static String name()
        {
          return "MeshAttribute<...>";
        }
    }; // class MeshAttribute<...>
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_MESH_ATTRIBUTE

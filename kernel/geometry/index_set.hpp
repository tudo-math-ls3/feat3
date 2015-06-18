#pragma once
#ifndef KERNEL_GEOMETRY_INDEX_SET_HPP
#define KERNEL_GEOMETRY_INDEX_SET_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
{
  /**
   * \brief Geometry namespace
   */
  namespace Geometry
  {
    /**
     * \brief Conformal Index-Set class template
     *
     * This is used for storing information about i.e. which vertices are present in each
     * - edge
     * - face
     * - cell
     *
     * and similar types of information, like edges\@faces, edges\@cells, faces\@cells, depending
     * on the dimension.
     *
     * \author Peter Zajac
     */
    template<int num_indices_>
    class IndexSet
    {
      static_assert(num_indices_ > 0, "invalid index count");

    public:
      /// number of indices per entry
      static constexpr int num_indices = num_indices_;

      /// index vector type
      typedef Index IndexVectorType[num_indices];

      /// index vector reference type
      typedef IndexVectorType& IndexVectorReference;

      /// const index vector reference type
      typedef const IndexVectorType& ConstIndexVectorReference;

      /// ImageIterator type for Adjunctor interface implementation
      typedef const Index* ImageIterator;

    protected:
      /// number of entities, i.e. number of cells
      Index _num_entities;
      /// index bound;, i.e. number of edges if the IndexSet is to hold the information edges\@cells.
      /// Necessary for Adjunctor implementation
      Index _index_bound;

      /// index vector array
      IndexVectorType* _indices;

    private:
      /// Copy assignment operator
      IndexSet& operator=(const IndexSet&);

    public:
      /**
       * \brief Constructor
       *
       * \param[in] num_entities
       * The number of entities for the index set.
       *
       * \param[in] index_bound
       * The index bound.
       */
      explicit IndexSet(
        Index num_entities = 0,
        Index index_bound = 0)
          :
        _num_entities(num_entities),
        _index_bound(index_bound),
        _indices(nullptr)
      {
        CONTEXT(name() + "::IndexSet()");
        if(num_entities > 0)
        {
          _indices = new IndexVectorType[num_entities];
        }
      }

      /**
       * \brief Copy Constructor
       *
       * \param[in] other
       * The index set that is to be copied.
       */
      IndexSet(const IndexSet& other) :
        _num_entities(other._num_entities),
        _index_bound(other._index_bound),
        _indices(nullptr)
      {
        if(_num_entities > 0)
        {
          _indices = new IndexVectorType[_num_entities];
          for(Index i(0); i < _num_entities; ++i)
          {
            for(int j(0); j < num_indices_; ++j)
            {
              _indices[i][j] = other._indices[i][j];
            }
          }
        }
      }

      /// virtual destructor
      virtual ~IndexSet()
      {
        CONTEXT(name() + "::~IndexSet()");
        if(_indices != nullptr)
        {
          delete [] _indices;
        }
      }

      /// Move assignment operator
      IndexSet& operator=(IndexSet&& other)
      {
        if(_indices != nullptr)
          delete[] _indices;

        _indices = (other._indices);

        _num_entities = other._num_entities;
        _index_bound = other._index_bound;

        other._indices = nullptr;
        other._num_entities = 0;
        other._index_bound = 0;
        return *this;
      }

      /**
       * \brief Returns the number of indices per entity.
       * \returns
       * The number of indices per entitiy.
       */
      int get_num_indices() const
      {
        return num_indices_;
      }

      /**
       * \brief Returns the number of entities.
       * \returns
       * The number of entities;
       */
      Index get_num_entities() const
      {
        CONTEXT(name() + "::get_num_entities()");
        return _num_entities;
      }

      /**
       * \brief Returns the index bound.
       * \returns
       * The index bound.
       */
      Index get_index_bound() const
      {
        CONTEXT(name() + "::get_index_bound()");
        return _index_bound;
      }

      /**
       * \brief Returns the index vector array.
       * \returns
       * The index vector array.
       */
      IndexVectorType* get_indices()
      {
        CONTEXT(name() + "::get_indices()");
        return _indices;
      }

      /** \copydoc get_indices() */
      const IndexVectorType* get_indices() const
      {
        CONTEXT(name() + "::get_indices() [const]");
        return _indices;
      }

      /**
       * \brief Returns a reference to an index vector.
       *
       * \param[in] i
       * The index of the entity whose index vector is to be returned.
       *
       * \returns
       * A (const) reference to the index vector of the entity.
       */
      IndexVectorReference operator[](Index i)
      {
        CONTEXT(name() + "::operator[]()");
        ASSERT_(_indices != nullptr);
        ASSERT_(i < _num_entities);
        return _indices[i];
      }

      /** \copydoc operator[]() */
      ConstIndexVectorReference operator[](Index i) const
      {
        CONTEXT(name() + "::operator[]() [const]");
        ASSERT_(_indices != nullptr);
        ASSERT_(i < _num_entities);
        return _indices[i];
      }

      /**
       * \brief Maps a face index.
       *
       * \param[in] i
       * The index of the entity whose face index is to be returned.
       *
       * \param[in] j
       * The index of the local face that is to be mapped.
       *
       * \returns
       * A (const) reference to the index of the local face \p j on entity \p i.
       */
      Index& operator()(Index i, Index j)
      {
        CONTEXT(name() + "::operator()()");
        ASSERT_(_indices != nullptr);
        ASSERT_(i < _num_entities);
        ASSERT_(j < Index(num_indices));
        return _indices[i][j];
      }

      // copy-paste documentation because \copydoc does not work with operators
      /**
       * \brief Maps a face index.
       *
       * \param[in] i
       * The index of the entity whose face index is to be returned.
       *
       * \param[in] j
       * The index of the local face that is to be mapped.
       *
       * \returns
       * A (const) reference to the index of the local face \p j on entity \p i.
       */
      const Index& operator()(Index i, Index j) const
      {
        CONTEXT(name() + "::operator()()");
        ASSERT_(_indices != nullptr);
        ASSERT_(i < _num_entities);
        ASSERT_(j < Index(num_indices));
        return _indices[i][j];
      }

      /**
       * \brief Sets the index bound.
       *
       * \param[in] bound
       * The new index bound.
       */
      void set_index_bound(Index bound)
      {
        CONTEXT(name() + "::set_index_bound()");
        _index_bound = bound;
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "IndexSet<" + stringify(num_indices_) + ">";
      }

      /* *************************************************************************************** */
      /*          A D J A C T O R     I N T E R F A C E     I M P L E M E N T A T I O N          */
      /* *************************************************************************************** */
      /** \copydoc Adjacency::Adjactor::get_num_nodes_domain() */
      Index get_num_nodes_domain() const
      {
        return _num_entities;
      }

      /** \copydoc Adjacency::Adjactor::get_num_nodes_image() */
      Index get_num_nodes_image() const
      {
        return _index_bound;
      }

      /** \copydoc Adjacency::Adjactor::image_begin() */
      ImageIterator image_begin(Index domain_node) const
      {
        CONTEXT(name() + "::image_begin()");
        ASSERT_(_indices != nullptr);
        ASSERT_(domain_node < _num_entities);
        return &_indices[domain_node][0];
      }

      /** \copydoc Adjacency::Adjactor::image_end() */
      ImageIterator image_end(Index domain_node) const
      {
        CONTEXT(name() + "::image_end()");
        ASSERT_(_indices != nullptr);
        ASSERT_(domain_node < _num_entities);
        return &_indices[domain_node][num_indices];
      }
    }; // class IndexSet<...>

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    /// \cond internal
    template<
      typename Shape_,
      int face_dim_ = Shape_::dimension - 1>
    class IndexSetWrapper :
      public IndexSetWrapper<Shape_, face_dim_ - 1>
    {
      static_assert(face_dim_ < Shape_::dimension, "invalid face dimension");
      static_assert(face_dim_ > 0, "invalid face dimension");

    public:
      typedef IndexSetWrapper<Shape_, face_dim_ - 1> BaseClass;
      typedef IndexSet<Shape::FaceTraits<Shape_, face_dim_>::count> IndexSetType;

    protected:
      IndexSetType _index_set;

    public:
      explicit IndexSetWrapper(const Index num_entities[]) :
        BaseClass(num_entities),
        _index_set(
          num_entities[Shape_::dimension],
          num_entities[face_dim_])
      {
        CONTEXT(name() + "::IndexSetWrapper(const Index*)");
      }

      explicit IndexSetWrapper() :
        BaseClass(),
        _index_set()
      {
        CONTEXT(name() + "::IndexSetWrapper()");
      }

      IndexSetWrapper(const IndexSetWrapper& other) :
        BaseClass(other),
        _index_set(other._index_set)
      {
        CONTEXT(name() + "::IndexSetWrapper(const IndexSetWrapper&)");
      }

      virtual ~IndexSetWrapper()
      {
        CONTEXT(name() + "::~IndexSetWrapper()");
      }

      template<int face_dim__>
      IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set()
      {
        CONTEXT(name() + "::get_index_set()");
        static_assert(face_dim__ >= 0, "invalid face dimension");
        static_assert(face_dim__ < Shape_::dimension, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      template<int face_dim__>
      const IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set() const
      {
        CONTEXT(name() + "::get_index_set() [const]");
        static_assert(face_dim__ >= 0, "invalid face dimension");
        static_assert(face_dim__ < Shape_::dimension, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      static String name()
      {
        return "IndexSetWrapper<" + Shape_::name() + "," + stringify(face_dim_) + ">";
      }
    };

    template<typename Shape_>
    class IndexSetWrapper<Shape_, 0>
    {
      static_assert(Shape_::dimension > 0, "invalid shape dimension");

    public:
      typedef IndexSet<Shape::FaceTraits<Shape_, 0>::count> IndexSetType;

    protected:
      IndexSetType _index_set;

    public:
      explicit IndexSetWrapper(const Index num_entities[]) :
        _index_set(
          num_entities[Shape_::dimension],
          num_entities[0])
      {
        CONTEXT(name() + "::IndexSetWrapper(const Index*)");
      }

      IndexSetWrapper(const IndexSetWrapper& other) :
        _index_set(other._index_set)
      {
        CONTEXT(name() + "::IndexSetWrapper(const IndexSetWrapper&)");
      }

      explicit IndexSetWrapper() :
        _index_set()
      {
        CONTEXT(name() + "::IndexSetWrapper()");
      }

      virtual ~IndexSetWrapper()
      {
        CONTEXT(name() + "::~IndexSetWrapper()");
      }

      template<int face_dim__>
      IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set()
      {
        CONTEXT(name() + "::get_index_set<" + stringify(face_dim__) + ">()");
        static_assert(face_dim__ == 0, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      template<int face_dim__>
      const IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set() const
      {
        CONTEXT(name() + "::get_index_set<" + stringify(face_dim__) + ">() [const]");
        static_assert(face_dim__ == 0, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      static String name()
      {
        return "IndexSetWrapper<" + Shape_::name() + ",0>";
      }
    };

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    /**
     * \brief Class for holding a cross product of IndexSets
     *
     * \tparam Shape_
     * The shape type of the highest dimensional cell type, i.e. Hypercube<3>
     *
     * Depending on the shape type and the dimension, the mesh consists of objects associated with every dimension,
     * i.e. dim =
     *  0: vertices
     *  1: edges
     *  2: faces (Shape_::dimension==3) or cells (Shape_::dimension==2)
     *  3: cells.
     *
     * An IndexSetHolder stores the complete adjacency structure in a matrix:
     *
     * @      |vertex|edge|face|cell
     * -----------------------------
     * vertex |        v@e  v@f  v@c
     * edge   |             e@f  e@c
     * face   |                  f@c
     * cell   |
     *
     * Each entry in the IndexSetHolder is an IndexSet containing information of the type:
     * For every cell: Which edges are present in that cell (e@c) etc.
     *
     *
     */
    template<typename Shape_>
    class IndexSetHolder :
      public IndexSetHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType>
    {
    public:
      typedef IndexSetHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType> BaseClass;
      typedef IndexSetWrapper<Shape_> IndexSetWrapperType;

      template<
        int cell_dim_,
        int face_dim_>
      struct IndexSet
      {
        /// index set type
        typedef Geometry::IndexSet
        <
          Shape::FaceTraits
          <
            typename Shape::FaceTraits<Shape_, cell_dim_> ::ShapeType,
            face_dim_
          >
          ::count
        > Type;
      }; // struct IndexSet<...>

    protected:
      IndexSetWrapperType _index_set_wrapper;

    public:
      explicit IndexSetHolder(const Index num_entities[]) :
        BaseClass(num_entities),
        _index_set_wrapper(num_entities)
      {
        CONTEXT(name() + "::IndexSetHolder(Index*)");
      }

      IndexSetHolder(const IndexSetHolder& other) :
        BaseClass(other),
        _index_set_wrapper(other._index_set_wrapper)
      {
        CONTEXT(name() + "::IndexSetHolder(IndexSetHolder&)");
      }

      explicit IndexSetHolder() :
        BaseClass(),
        _index_set_wrapper()
        {
          CONTEXT(name() + "::IndexSetHolder()");
        }

      virtual ~IndexSetHolder()
      {
        CONTEXT(name() + "::~IndexSetHolder()");
      }

      template<int shape_dim_>
      IndexSetWrapper<typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType>& get_index_set_wrapper()
      {
        CONTEXT(name() + "::get_index_set_wrapper<" + stringify(shape_dim_) + ">()");
        static_assert(shape_dim_ > 0, "invalid shape dimension");
        static_assert(shape_dim_ <= Shape_::dimension, "invalid shape dimension");
        typedef typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType ShapeType;
        return IndexSetHolder<ShapeType>::_index_set_wrapper;
      }

      template<int shape_dim_>
      const IndexSetWrapper<typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType>& get_index_set_wrapper() const
      {
        CONTEXT(name() + "::get_index_set_wrapper<" + stringify(shape_dim_) + ">() [const]");
        static_assert(shape_dim_ > 0, "invalid shape dimension");
        static_assert(shape_dim_ <= Shape_::dimension, "invalid shape dimension");
        typedef typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType ShapeType;
        return IndexSetHolder<ShapeType>::_index_set_wrapper;
      }

      template<
        int cell_dim_,
        int face_dim_>
      Geometry::IndexSet<
        Shape::FaceTraits<
          typename Shape::FaceTraits<
            Shape_,
            cell_dim_>
          ::ShapeType,
          face_dim_>
        ::count>&
      get_index_set()
      {
        CONTEXT(name() + "::get_index_set<" + stringify(cell_dim_) + "," + stringify(face_dim_) + ">");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");
        static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
        static_assert(face_dim_ >= 0, "invalid face dimension");
        return get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
      }

      template<
        int cell_dim_,
        int face_dim_>
      const Geometry::IndexSet<
        Shape::FaceTraits<
          typename Shape::FaceTraits<
            Shape_,
            cell_dim_>
          ::ShapeType,
          face_dim_>
        ::count>&
      get_index_set() const
      {
        CONTEXT(name() + "::get_index_set<" + stringify(cell_dim_) + "," + stringify(face_dim_) + "> [const]");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");
        static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
        static_assert(face_dim_ >= 0, "invalid face dimension");
        return get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
      }

      static String name()
      {
        return "IndexSetHolder<" + Shape_::name() + ">";
      }
    };

    template<>
    class IndexSetHolder<Shape::Vertex>
    {
    public:
      explicit IndexSetHolder(const Index* /*num_entities*/)
      {
        CONTEXT(name() + "::IndexSetHolder(Index*)");
      }

      IndexSetHolder(const IndexSetHolder&)
      {
        CONTEXT(name() + "::IndexSetHolder(IndexSetHolder&)");
      }

      explicit IndexSetHolder()
      {
        CONTEXT(name() + "::IndexSetHolder()");
      }

      virtual ~IndexSetHolder()
      {
        CONTEXT(name() + "::~IndexSetHolder()");
      }

      static String name()
      {
        return "IndexSetHolder<Vertex>";
      }
    };

    /**
     * \brief Helper class to update the sizes of elements of an IndexSetHolder
     *
     * \tparam shape_dim_
     * Dimension of highest dimensional shape
     *
     * \tparam co_dim_
     * Codimension wrt. the subshape that causes the size updates.
     *
     * When using meshes that contain only the vertex@shape information, the number of edges/faces is unknown so the
     * corresponding IndexSets in the IndexSetHolder are initialised with size 0. If the IndexSet to dimensions
     * (shape_dim_-co_dim_, 0) is created, all IndexSets using the number of entities of dimension shape_dim_-co_dim_
     * need to be updated.
     *
     * To make this more clear:
     * _shape_dim = 1: Nothing to do, as vertex@shape is already present.
     * _shape_dim = 2: After vert@edge is computed, update the size of edge@cell.
     * _shape_dim = 3: After vert@face is computed, update the size of face@cell.
     *                 After vert@edge is computed, update the sizes of edge@cell and edge@face.
     *
     * \warning Because of the information hierarchy, these routines assume that the vertex@shape information is
     * present and correct for the highest dimensional shape and for all shapes with lower codimension!
     *
     * \author Jordi Paul
     *
     */
    template<int shape_dim_, int co_dim_>
    struct IndexSetHolderDimensionUpdater
    {
      /**
       * \brief Updates the size of certain IndexSets in an IndexSetHolder depending on shape_dim_ and co_dim_
       *
       * \tparam IndexSetHolderType_
       * Type of the IndexSetHolder
       *
       * \param[in,out] ish
       * IndexSetHolder whose IndexSets get updated
       */
      template<typename IndexSetHolderType_>
      static void update(IndexSetHolderType_& DOXY(ish))
      {
        //dummy
      }
    };

    /**
     * \brief Specialisation for codimension 1
     *
     */
    template<int shape_dim_>
    struct IndexSetHolderDimensionUpdater<shape_dim_,1>
    {
      template<typename IndexSetHolderType_>
      static void update(IndexSetHolderType_& ish)
      {
        auto& vert_at_subshape_index_set(ish.template get_index_set<shape_dim_-1,0>());

        Index num_shapes(ish.template get_index_set<shape_dim_,0>().get_num_entities());
        // Get the type of the IndexSet to dimensions <Shape_::dimension, face_dim_>
        typename std::remove_reference<decltype( (ish.template get_index_set<shape_dim_, shape_dim_-1>()) )>::type
          // The output IndexSet has num_shapes entities and the maximum index is the number of vertices
          subshape_at_shape_index_set(num_shapes, vert_at_subshape_index_set.get_num_entities());

        // Replace the corresponding IndexSet in the IndexSetHolder by the just generated IndexSet of the right
        // size
        ish.template get_index_set<shape_dim_, shape_dim_-1>() = std::move(subshape_at_shape_index_set);
      }
    };

    /**
     * \brief Specialisation for codimension 2
     *
     */
    template<int shape_dim_>
    struct IndexSetHolderDimensionUpdater<shape_dim_,2>
    {
      template<typename IndexSetHolderType_>
      static void update(IndexSetHolderType_& ish)
      {
        // Update edge@cell IndexSet dimensions
        // vert@edge IndexSet
        auto& vert_at_subshape_index_set(ish.template get_index_set<shape_dim_-2,0>());
        // Number of cells
        Index num_shapes(ish.template get_index_set<shape_dim_,0>().get_num_entities());
        // Get the type of the IndexSet to dimensions <Shape_::dimension, Shape::dimension-2>
        typename std::remove_reference<decltype( (ish.template get_index_set<shape_dim_, shape_dim_-2>()) )>::type
          // The output IndexSet has num_shapes entities and the maximum index is the number of subshapes (=edges)
          subshape_at_shape_index_set(num_shapes, vert_at_subshape_index_set.get_num_entities());
        // Replace the corresponding IndexSet in the IndexSetHolder by the just generated IndexSet of the right
        // size
        ish.template get_index_set<shape_dim_, shape_dim_-2>() = std::move(subshape_at_shape_index_set);

        // Now do it again for the edge@face IndexSet
        // vert@face IndexSet
        auto& vert_at_intermediate_index_set(ish.template get_index_set<shape_dim_-1,0>());
        // Number of faces
        Index num_intermediate_shapes(ish.template get_index_set<shape_dim_-1,0>().get_num_entities());
        // Get the type of the IndexSet to dimensions <Shape_::dimension-1, Shape::dimension-2>
        typename std::remove_reference<decltype( (ish.template get_index_set<shape_dim_-1, shape_dim_-2>()) )>::type
          // The output IndexSet has num_intermediate_shapes entities and the maximum index is the number of edges
          subshape_at_intermediate_index_set(
            num_intermediate_shapes, vert_at_intermediate_index_set.get_num_entities());
        // Replace the corresponding IndexSet in the IndexSetHolder by the just generated IndexSet of the right
        // size
        ish.template get_index_set<shape_dim_-1, shape_dim_-2>() = std::move(subshape_at_intermediate_index_set);
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INDEX_SET_HPP

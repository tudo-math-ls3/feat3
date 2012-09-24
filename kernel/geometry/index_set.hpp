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
     * \todo detailed documentation
     *
     * \author Peter Zajac
     */
    template<int num_indices_>
    class IndexSet
    {
      static_assert(num_indices_ > 0, "invalid index count");

    public:
      /// dummy enum
      enum
      {
        /// number of indices per entry
        num_indices = num_indices_
      };

      /// index vector type
      typedef Index IndexVectorType[num_indices];

      /// index vector reference type
      typedef IndexVectorType& IndexVectorReference;

      /// const index vector reference type
      typedef const IndexVectorType& ConstIndexVectorReference;

      /// ImageIterator type for Adjunctor interface implementation
      typedef const Index* ImageIterator;

    protected:
      /// number of entities
      Index _num_entities;
      /// index bound; necessary for Adjunctor implementation
      Index _index_bound;

      /// index vector array
      IndexVectorType* _indices;

    private:
      IndexSet(const IndexSet&);
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

      /// virtual destructor
      virtual ~IndexSet()
      {
        CONTEXT(name() + "::~IndexSet()");
        if(_indices != nullptr)
        {
          delete [] _indices;
        }
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
      /** \copydoc Adjactor::get_num_nodes_domain() */
      Index get_num_nodes_domain() const
      {
        return _num_entities;
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      Index get_num_nodes_image() const
      {
        return _index_bound;
      }

      /** \copydoc Adjactor::image_begin() */
      ImageIterator image_begin(Index domain_node) const
      {
        CONTEXT(name() + "::image_begin()");
        ASSERT_(_indices != nullptr);
        ASSERT_(domain_node < _num_entities);
        return &_indices[domain_node][0];
      }

      /** \copydoc Adjactor::image_end() */
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
        CONTEXT(name() + "::IndexSetWrapper()");
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

    template<typename Shape_>
    class IndexSetHolder :
      public IndexSetHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType>
    {
    public:
      typedef IndexSetHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType> BaseClass;
      typedef IndexSetWrapper<Shape_> IndexSetWrapperType;

    protected:
      IndexSetWrapperType _index_set_wrapper;

    public:
      explicit IndexSetHolder(const Index num_entities[]) :
        BaseClass(num_entities),
        _index_set_wrapper(num_entities)
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
      IndexSet<
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
        return get_index_set_wrapper<cell_dim_>().get_index_set<face_dim_>();
      }

      template<
        int cell_dim_,
        int face_dim_>
      const IndexSet<
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
        return get_index_set_wrapper<cell_dim_>().get_index_set<face_dim_>();
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
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INDEX_SET_HPP

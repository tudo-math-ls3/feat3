// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/shape.hpp>
#include <kernel/adjacency/permutation.hpp>

// includes, system
#include <array>
#include <vector>

namespace FEAT
{
  /**
   * \brief Geometry namespace
   */
  namespace Geometry
  {
    /**
     * \brief Index Tuple class template
     *
     * This class is a simple wrapper around a fixed-size array,
     * similar to the std::array class template.
     * Its primary purpose is to avoid compiler warnings and error
     * when used as a value-type in an std::vector.
     *
     * \author Peter Zajac
     */
    template<int num_indices_>
    struct IndexTuple
    {
      static_assert(num_indices_ > 0, "invalid number of indices");

      /// number of indices per tuple
      static constexpr int num_indices = num_indices_;

      /// indices array
      Index indices[num_indices];

      /// access operator
      Index& operator[](int i)
      {
        ASSERT(i >= 0);
        ASSERT(i < num_indices);
        return indices[i];
      }

      /// access operator
      const Index& operator[](int i) const
      {
        ASSERT(i >= 0);
        ASSERT(i < num_indices);
        return indices[i];
      }

      void permute_map(const Adjacency::Permutation& inv_perm)
      {
        for(int i(0); i < num_indices; ++i)
          indices[i] = inv_perm.map(indices[i]);
      }
    }; // struct IndexTuple

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

      /// index tuple type
      typedef IndexTuple<num_indices> IndexTupleType;

      /// ImageIterator type for Adjactor interface implementation
      typedef const Index* ImageIterator;

    protected:
      /// index bound;, i.e. number of edges if the IndexSet is to hold the information edges\@cells.
      /// Necessary for Adjactor implementation
      Index _index_bound;

      /// index vector array
      std::vector<IndexTupleType> _indices;

      /// internal clone constructor
      explicit IndexSet(Index idx_bnd, const std::vector<IndexTupleType>& idx) :
        _index_bound(idx_bnd), _indices(idx)
      {
      }

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
        _index_bound(index_bound),
        _indices(std::size_t(num_entities))
      {
      }

      /// move constructor
      IndexSet(IndexSet&& other) :
        _index_bound(other._index_bound),
        _indices(std::forward<std::vector<IndexTupleType>>(other._indices))
      {
      }

      /// move-assignment operator
      IndexSet& operator=(IndexSet&& other)
      {
        // avoid self-move
        if(this == &other)
          return *this;

        _index_bound = other._index_bound;
        _indices = std::forward<std::vector<IndexTupleType>>(other._indices);

        return *this;
      }

      /// virtual destructor
      virtual ~IndexSet()
      {
      }

      /// \returns An independent clone of this target set object.
      IndexSet clone() const
      {
        return IndexSet(this->_index_bound, this->_indices);
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return _indices.size() * sizeof(IndexTupleType);
      }

      /**
       * \brief Returns the number of indices per entity.
       * \returns
       * The number of indices per entity.
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
        return Index(_indices.size());
      }

      /**
       * \brief Returns the index bound.
       * \returns
       * The index bound.
       */
      Index get_index_bound() const
      {
        return _index_bound;
      }

      /**
       * \brief Returns the index vector array.
       * \returns
       * The index vector array.
       */
      IndexTupleType* get_indices()
      {
        return _indices.data();
      }

      /** \copydoc get_indices() */
      const IndexTupleType* get_indices() const
      {
        return _indices.data();
      }

      /**
       * \brief Returns a reference to an index tuple.
       *
       * \param[in] i
       * The index of the entity whose index tuple is to be returned.
       *
       * \returns
       * A (const) reference to the index tuple of the entity.
       */
      IndexTupleType& operator[](Index i)
      {
        ASSERT(i < get_num_entities());
        return _indices[i];
      }

      /** \copydoc operator[]() */
      const IndexTupleType& operator[](Index i) const
      {
        ASSERT(i < get_num_entities());
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
      Index& operator()(Index i, int j)
      {
        ASSERT(i < get_num_entities());
        ASSERT(j < num_indices);
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
      const Index& operator()(Index i, int j) const
      {
        ASSERT(i < get_num_entities());
        ASSERT(j < num_indices);
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
        _index_bound = bound;
      }

      /**
       * \brief Permutes this index set in-situ.
       *
       * \perm[in] perm
       * A \transient reference to the permutation of the shape indices.
       *
       * \param[in] inv_perm_face
       * A \transient reference to the inverse permutation of the face indices.
       */
      void permute(const Adjacency::Permutation& perm, const Adjacency::Permutation& inv_perm_face)
      {
        if(_indices.empty())
          return;

        // first, permute the actual index tuples by the forward permutation
        if(!perm.empty())
        {
          XASSERT(perm.size() == get_num_entities());
          perm.apply(_indices.data());
        }

        // now, permute all indices by the inverse face permutation
        if(!inv_perm_face.empty())
        {
          XASSERT(inv_perm_face.size() == get_index_bound());
          const Index n = Index(_indices.size());
          for(Index i(0); i < n; ++i)
            _indices.at(i).permute_map(inv_perm_face);
        }
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
        return get_num_entities();
      }

      /** \copydoc Adjacency::Adjactor::get_num_nodes_image() */
      Index get_num_nodes_image() const
      {
        return _index_bound;
      }

      /** \copydoc Adjacency::Adjactor::image_begin() */
      ImageIterator image_begin(Index domain_node) const
      {
        ASSERT(domain_node < get_num_entities());
        return reinterpret_cast<const Index*>(&_indices.data()[domain_node]);
      }

      /** \copydoc Adjacency::Adjactor::image_end() */
      ImageIterator image_end(Index domain_node) const
      {
        ASSERT(domain_node < get_num_entities());
        return reinterpret_cast<const Index*>(&_indices.data()[domain_node+1]);
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
      IndexSetWrapper() :
        BaseClass(),
        _index_set()
      {
      }

      explicit IndexSetWrapper(const Index num_entities[]) :
        BaseClass(num_entities),
        _index_set(
          num_entities[Shape_::dimension],
          num_entities[face_dim_])
      {
      }

      template<std::size_t m_>
      explicit IndexSetWrapper(const std::array<Index, m_>& num_entities) :
        BaseClass(num_entities),
        _index_set(
          num_entities[Shape_::dimension],
          num_entities[face_dim_])
      {
        static_assert(int(m_) > Shape_::dimension, "invalid array size");
      }

      IndexSetWrapper(IndexSetWrapper&& other) :
        BaseClass(std::forward<BaseClass>(other)),
        _index_set(std::forward<IndexSetType>(other._index_set))
      {
      }

      IndexSetWrapper& operator=(IndexSetWrapper&& other)
      {
        if(this == &other)
          return *this;

        BaseClass::operator=(std::forward<BaseClass>(other));
        _index_set = std::forward<IndexSetType>(other._index_set);

        return *this;
      }

      virtual ~IndexSetWrapper()
      {
      }

      void clone(const IndexSetWrapper& other)
      {
        BaseClass::clone(other);
        this->_index_set = other._index_set.clone();
      }

      IndexSetWrapper clone() const
      {
        IndexSetWrapper isw;
        isw.clone(*this);
        return isw;
      }

      template<int face_dim__>
      IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set()
      {
        static_assert(face_dim__ >= 0, "invalid face dimension");
        static_assert(face_dim__ < Shape_::dimension, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      template<int face_dim__>
      const IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set() const
      {
        static_assert(face_dim__ >= 0, "invalid face dimension");
        static_assert(face_dim__ < Shape_::dimension, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      template<std::size_t np_>
      void permute(const Adjacency::Permutation& shape_perm,
        const std::array<Adjacency::Permutation, np_>& inv_perm)
      {
        BaseClass::permute(shape_perm, inv_perm);
        _index_set.permute(shape_perm, inv_perm.at(face_dim_));
      }

      static String name()
      {
        return "IndexSetWrapper<" + Shape_::name() + "," + stringify(face_dim_) + ">";
      }

      std::size_t bytes() const
      {
        return BaseClass::bytes() + _index_set.bytes();
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
      IndexSetWrapper() :
        _index_set()
      {
      }

      explicit IndexSetWrapper(const Index num_entities[]) :
        _index_set(
          num_entities[Shape_::dimension],
          num_entities[0])
      {
      }

      template<std::size_t m_>
      explicit IndexSetWrapper(const std::array<Index, m_>& num_entities) :
        _index_set(
          num_entities[Shape_::dimension],
          num_entities[0])
      {
        static_assert(int(m_) > Shape_::dimension, "invalid array size");
      }

      IndexSetWrapper(IndexSetWrapper&& other) :
        _index_set(std::forward<IndexSetType>(other._index_set))
      {
      }

      IndexSetWrapper& operator=(IndexSetWrapper&& other)
      {
        if(this == &other)
          return *this;

        _index_set = std::forward<IndexSetType>(other._index_set);

        return *this;
      }

      virtual ~IndexSetWrapper()
      {
      }

      void clone(const IndexSetWrapper& other)
      {
        this->_index_set = other._index_set.clone();
      }

      IndexSetWrapper clone() const
      {
        IndexSetWrapper isw;
        isw.clone(*this);
        return isw;
      }

      template<int face_dim__>
      IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set()
      {
        static_assert(face_dim__ == 0, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      template<int face_dim__>
      const IndexSet<Shape::FaceTraits<Shape_, face_dim__>::count>& get_index_set() const
      {
        static_assert(face_dim__ == 0, "invalid face dimension");
        return IndexSetWrapper<Shape_, face_dim__>::_index_set;
      }

      template<std::size_t np_>
      void permute(const Adjacency::Permutation& shape_perm,
        const std::array<Adjacency::Permutation, np_>& inv_perm)
      {
        _index_set.permute(shape_perm, inv_perm.at(0));
      }

      static String name()
      {
        return "IndexSetWrapper<" + Shape_::name() + ",0>";
      }

      std::size_t bytes() const
      {
        return _index_set.bytes();
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
     * - 0: vertices
     * - 1: edges
     * - 2: faces
     * - 3: cells.
     *
     * An IndexSetHolder stores the complete adjacency structure in a matrix:
     *
       \verbatim
       @      |vertex|edge|face|cell
       -----------------------------
       vertex |        v@e  v@f  v@c
       edge   |             e@f  e@c
       face   |                  f@c
       cell   |
       \endverbatim
     *
     * Each entry in the IndexSetHolder is an IndexSet containing information of the type:
     * For every cell: Which edges are present in that cell (e\@c) etc.
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
      IndexSetHolder() :
        BaseClass(),
        _index_set_wrapper()
      {
      }

      explicit IndexSetHolder(const Index num_entities[]) :
        BaseClass(num_entities),
        _index_set_wrapper(num_entities)
      {
      }

      template<std::size_t m_>
      explicit IndexSetHolder(const std::array<Index, m_>& num_entities) :
        BaseClass(num_entities),
        _index_set_wrapper(num_entities)
      {
        static_assert(int(m_) > Shape_::dimension, "invalid array size");
      }

      IndexSetHolder(IndexSetHolder&& other) :
        BaseClass(std::forward<BaseClass>(other)),
        _index_set_wrapper(std::forward<IndexSetWrapperType>(other._index_set_wrapper))
      {
      }

      IndexSetHolder& operator=(IndexSetHolder&& other)
      {
        if(this == &other)
          return *this;

        BaseClass::operator=(std::forward<BaseClass>(other));
        _index_set_wrapper = std::forward<IndexSetWrapperType>(other._index_set_wrapper);

        return *this;
      }

      virtual ~IndexSetHolder()
      {
      }

      void clone(const IndexSetHolder& other)
      {
        BaseClass::clone(other);
        this->_index_set_wrapper.clone(other._index_set_wrapper);
      }

      IndexSetHolder clone() const
      {
        IndexSetHolder ish;
        ish.clone(*this);
        return ish;
      }

      template<int shape_dim_>
      IndexSetWrapper<typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType>& get_index_set_wrapper()
      {
        static_assert(shape_dim_ > 0, "invalid shape dimension");
        static_assert(shape_dim_ <= Shape_::dimension, "invalid shape dimension");
        typedef typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType ShapeType;
        return IndexSetHolder<ShapeType>::_index_set_wrapper;
      }

      template<int shape_dim_>
      const IndexSetWrapper<typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType>& get_index_set_wrapper() const
      {
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
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");
        static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
        static_assert(face_dim_ >= 0, "invalid face dimension");
        return get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
      }

      template<std::size_t np_>
      void permute(const std::array<Adjacency::Permutation, np_>& perms,
        const std::array<Adjacency::Permutation, np_>& inv_perms)
      {
        BaseClass::permute(perms, inv_perms);
        _index_set_wrapper.permute(perms.at(Shape_::dimension), inv_perms);
      }

      static String name()
      {
        return "IndexSetHolder<" + Shape_::name() + ">";
      }

      std::size_t bytes() const
      {
        return BaseClass::bytes() + _index_set_wrapper.bytes();
      }
    };

    template<>
    class IndexSetHolder<Shape::Vertex>
    {
    public:
      IndexSetHolder()
      {
      }

      explicit IndexSetHolder(const Index* /*num_entities*/)
      {
      }

      template<std::size_t m_>
      explicit IndexSetHolder(const std::array<Index, m_>&)
      {
      }

      IndexSetHolder(IndexSetHolder&&)
      {
      }

      IndexSetHolder& operator=(IndexSetHolder&&)
      {
        return *this;
      }

      virtual ~IndexSetHolder()
      {
      }

      void clone(const IndexSetHolder&)
      {
      }

      IndexSetHolder clone() const
      {
        return IndexSetHolder();
      }

      template<std::size_t np_>
      void permute(const std::array<Adjacency::Permutation, np_>&,
        const std::array<Adjacency::Permutation, np_>&)
      {
      }

      static String name()
      {
        return "IndexSetHolder<Vertex>";
      }

      std::size_t bytes() const
      {
        return std::size_t(0);
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
     * When using meshes that contain only the vertex\@shape information, the number of edges/faces is unknown so the
     * corresponding IndexSets in the IndexSetHolder are initialized with size 0. If the IndexSet to dimensions
     * (shape_dim_-co_dim_, 0) is created, all IndexSets using the number of entities of dimension shape_dim_-co_dim_
     * need to be updated.
     *
     * To make this more clear:
     * _shape_dim = 1: Nothing to do, as vertex\@shape is already present.
     * _shape_dim = 2: After vert\@edge is computed, update the size of edge\@cell.
     * _shape_dim = 3: After vert\@face is computed, update the size of face\@cell.
     *                 After vert\@edge is computed, update the sizes of edge\@cell and edge\@face.
     *
     * \warning Because of the information hierarchy, these routines assume that the vertex\@shape information is
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
     * \brief Specialization for codimension 1
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
          // The output IndexSet has num_shapes entities and the maximum index is the number of subshapes
          subshape_at_shape_index_set(num_shapes, vert_at_subshape_index_set.get_num_entities());

        // Replace the corresponding IndexSet in the IndexSetHolder by the just generated IndexSet of the right
        // size
        ish.template get_index_set<shape_dim_, shape_dim_-1>() = std::move(subshape_at_shape_index_set);
      }
    };

    /**
     * \brief Specialization for codimension 2
     *
     */
    template<int shape_dim_>
    struct IndexSetHolderDimensionUpdater<shape_dim_,2>
    {
      template<typename IndexSetHolderType_>
      static void update(IndexSetHolderType_& ish)
      {
        // Update edge\@cell IndexSet dimensions
        // vert\@edge IndexSet
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

        // Now do it again for the edge\@face IndexSet
        // Number of faces
        Index num_intermediate_shapes(ish.template get_index_set<shape_dim_-1,0>().get_num_entities());
        // Get the type of the IndexSet to dimensions <Shape_::dimension-1, Shape::dimension-2>
        typename std::remove_reference<decltype( (ish.template get_index_set<shape_dim_-1, shape_dim_-2>()) )>::type
          // The output IndexSet has num_intermediate_shapes entities and the maximum index is the number of subshapes
          subshape_at_intermediate_index_set(
            num_intermediate_shapes, vert_at_subshape_index_set.get_num_entities());
        // Replace the corresponding IndexSet in the IndexSetHolder by the just generated IndexSet of the right
        // size
        ish.template get_index_set<shape_dim_-1, shape_dim_-2>() = std::move(subshape_at_intermediate_index_set);
      }
    };

    /**
     * \brief Helper class for extracting num_entities information from an IndexSetHolder
     *
     * \tparam dim_
     * Maximum dimension of the shape to extract said information for.
     *
     * If the mesh file read by the MeshStreamer used by the MeshStreamerFactory does not provide complete
     * information about edge\@cell, edge\@face and/or face\@cell, the number of edges/faces is not correct in
     * the _mesh_data object. After fill_index_sets() is called, the size of these index set in the corresponding
     * IndexSetHolder gives the right number of edges/faces.
     *
     * Because the get_index_set routine uses the dimensions as template parameters, the loop has to be realized
     * by template recursion and as the end of that recursion needs partial specialization, a helper class.
     *
     * \author Jordi Paul
     */
    template<int dim_>
    struct NumEntitiesExtractor
    {
      /**
       * \brief Computes certain entries of num_entities from an IndexSetHolder
       *
       * \tparam IndexSetHolderType_
       * Type of the IndexSetHolder
       *
       * \param[in] ish
       * IndexSetHolder to set num_entities from
       *
       * \param[out] num_entities
       * num_entities[d] = Number of objects of dimension d
       *
       */
      template<typename IndexSetHolderType_>
      static void set_num_entities(const IndexSetHolderType_& ish, Index* num_entities)
      {
        // The number of entities of dimension dim_ is the length of the vertex\@shape[dim_] IndexSet
        num_entities[dim_] = ish.template get_index_set<dim_,0>().get_num_entities();
        // Recurse down
        NumEntitiesExtractor<dim_-1>::set_num_entities(ish, num_entities);
      }

      template<typename IndexSetHolderType_>
      static void set_num_entities_with_numverts(const IndexSetHolderType_& ish, Index* num_entities)
      {
        // The number of entities of dimension dim_ is the length of the vertex\@shape[dim_] IndexSet
        num_entities[dim_] = ish.template get_index_set<dim_,0>().get_num_entities();
        // Recurse down
        NumEntitiesExtractor<dim_-1>::set_num_entities_with_numverts(ish, num_entities);
      }
    };

    /**
     * \brief Full specialization of NumEntitiesExtractor as end of the template recursion
     *
     * As num_entities[0] = number of vertices and this information must be present and correct, this stops at 1.
     *
     */
    template<>
    struct NumEntitiesExtractor<1>
    {
      template<typename IndexSetHolderType_>
      static void set_num_entities(const IndexSetHolderType_& ish, Index* num_entities)
      {
        num_entities[1] = ish.template get_index_set<1,0>().get_num_entities();
      }

      template<typename IndexSetHolderType_>
      static void set_num_entities_with_numverts(const IndexSetHolderType_& ish, Index* num_entities)
      {
        num_entities[0] = ish.template get_index_set<1,0>().get_index_bound();
        num_entities[1] = ish.template get_index_set<1,0>().get_num_entities();
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

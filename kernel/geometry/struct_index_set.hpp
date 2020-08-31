// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_STRUCT_INDEX_SET_HPP
#define KERNEL_GEOMETRY_STRUCT_INDEX_SET_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>
#include <kernel/geometry/intern/struct_num_entities.hpp>
#include <kernel/geometry/intern/struct_index_mapping.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Structured Index-Set class template
     *
     * This class provides the functionality of an IndexSet for the StructuredMesh class template.
     *
     * \author Peter Zajac
     */
    template<
      int shape_dim_,
      int cell_dim_,
      int face_dim_>
    class StructIndexSet
    {
      static_assert(shape_dim_ >= cell_dim_, "invalid shape dimension");
      static_assert(cell_dim_ > face_dim_, "invalid cell dimension");
      static_assert(face_dim_ >= 0, "invalid face dimension");

    public:
      /// shape type
      typedef Shape::Hypercube<shape_dim_> ShapeType;
      /// cell type
      typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
      /// number of indices
      static constexpr int num_indices = Shape::FaceTraits<CellType, face_dim_>::count;

      /**
       * \brief ImageIterator implementation for Adjactor interface
       */
      class ImageIterator
      {
      private:
        /// pointer to num_slices array of the index-set
        const Index* _num_slices;
        /// cell index
        Index _i;
        /// face index
        Index _j;

      public:
        /// default constructor
        ImageIterator() :
          _num_slices(nullptr),
          _i(0),
          _j(0)
        {
        }

        /// constructor
        explicit ImageIterator(const Index* num_slices, Index i, Index j) :
          _num_slices(num_slices),
          _i(i),
          _j(j)
        {
        }

        /** \copydoc Adjacency::Adjactor::ImageIterator::operator*() */
        Index operator*() const
        {
          return Intern::StructIndexMapping<shape_dim_, cell_dim_, face_dim_>::compute(_i, _j, _num_slices);
        }

        /** \copydoc Adjacency::Adjactor::ImageIterator::operator++() */
        ImageIterator& operator++()
        {
          ++_j;
          return *this;
        }

        /** \copydoc Adjacency::Adjactor::ImageIterator::operator!=() */
        bool operator!=(const ImageIterator& other) const
        {
          return (_i != other._i) || (_j != other._j);
        }
      }; // class ImageIterator

      /**
       * \brief Index tuple class for structured meshes
       */
      class IndexTupleType
      {
      public:
        /// number of indices per tuple
        static constexpr int num_indices = StructIndexSet::num_indices;

      private:
        /// pointer to num_slices array of the index-set
        const Index* _num_slices;
        /// cell index
        Index _i;

      public:
        /// default constructor
        IndexTupleType() :
          _num_slices(nullptr),
          _i(0)
        {
        }

        /// constructor
        explicit IndexTupleType(const Index* num_slices, Index i) :
          _num_slices(num_slices),
          _i(i)
        {
        }

        /// access operator
        Index operator[](int j) const
        {
          return Intern::StructIndexMapping<shape_dim_, cell_dim_, face_dim_>::compute(_i, Index(j), _num_slices);
        }
      }; // class IndexTupleType

    private:
      /// number of slices
      Index _num_slices[shape_dim_];
      /// number of entities
      Index _num_entities;
      /// index bound
      Index _index_bound;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_slices
       * Specifies the number of slices in each direction.
       */
      explicit StructIndexSet(const Index num_slices[]) :
        _num_entities(Intern::StructNumEntities<shape_dim_, cell_dim_>::compute(num_slices)),
        _index_bound(Intern::StructNumEntities<shape_dim_, face_dim_>::compute(num_slices))
      {
        for(int i(0); i < shape_dim_; ++i)
        {
          _num_slices[i] = num_slices[i];
        }
      }

      /// move constructor
      StructIndexSet(StructIndexSet&& other) :
        _num_entities(other._num_entities),
        _index_bound(other._index_bound)
      {
        for(int i(0); i < shape_dim_; ++i)
        {
          _num_slices[i] = other._num_slices[i];
        }
      }

      /// move-assignment operator
      StructIndexSet& operator=(StructIndexSet&& other)
      {
        // avoid self-move
        if(this == &other)
          return *this;

        _num_entities = other._num_entities;
        _index_bound = other._index_bound;
        for(int i(0); i < shape_dim_; ++i)
        {
          _num_slices[i] = other._num_slices[i];
        }
        return *this;
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return std::size_t(0);
      }

      /**
       * \brief Returns the number of indices per entity.
       * \returns
       * The number of indices per entity.
       */
      int get_num_indices() const
      {
        return num_indices;
      }

      /**
       * \brief Returns the number of entities.
       * \returns
       * The number of entities;
       */
      Index get_num_entities() const
      {
        return _num_entities;
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
       * \brief Returns an index tuple.
       *
       * \param[in] i
       * The index of the entity whose index tuple is to be returned.
       *
       * \returns
       * The index tuple of the entity.
       */
      IndexTupleType operator[](Index i) const
      {
        ASSERT(i < Index(_num_entities));
        return IndexTupleType(_num_slices, i);
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
       * The index of the local face \p j on entity \p i.
       */
      Index operator()(Index i, int j) const
      {
        ASSERT(i < Index(_num_entities));
        ASSERT(j < num_indices);
        return Intern::StructIndexMapping<shape_dim_, cell_dim_, face_dim_>::compute(i, Index(j), _num_slices);
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
        ASSERT(domain_node < Index(_num_entities));
        return ImageIterator(_num_slices, domain_node, 0);
      }

      /** \copydoc Adjacency::Adjactor::image_end() */
      ImageIterator image_end(Index domain_node) const
      {
        ASSERT(domain_node < Index(_num_entities));
        return ImageIterator(_num_slices, domain_node, num_indices);
      }
    }; // class StructIndexSet<...>


    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    /// \cond internal
    template<
      int shape_dim_,
      int cell_dim_ = shape_dim_,
      int face_dim_ = cell_dim_ - 1>
    class StructIndexSetWrapper :
      public StructIndexSetWrapper<shape_dim_, cell_dim_, face_dim_ - 1>
    {
    protected:
      typedef StructIndexSetWrapper<shape_dim_, cell_dim_, face_dim_ - 1> BaseClass;

      StructIndexSet<shape_dim_, cell_dim_, face_dim_> _index_set;

    public:
      explicit StructIndexSetWrapper(const Index num_slices[]) :
        BaseClass(num_slices),
        _index_set(num_slices)
      {
      }

      StructIndexSetWrapper(StructIndexSetWrapper&& other) :
        BaseClass(std::forward<BaseClass>(other)),
        _index_set(std::forward<StructIndexSet<shape_dim_, cell_dim_, face_dim_>>(other._index_set))
      {
      }

      StructIndexSetWrapper& operator=(StructIndexSetWrapper&& other)
      {
        if(this == &other)
          return *this;

        BaseClass::operator=(std::forward<BaseClass>(other));
        _index_set = std::forward<StructIndexSet<shape_dim_, cell_dim_, face_dim_>>(other._index_set);

        return *this;
      }

      template<int face_dim__>
      const StructIndexSet<shape_dim_, cell_dim_, face_dim__>& get_index_set() const
      {
        static_assert((face_dim__ >= 0) && (face_dim__ <= face_dim_), "invalid face dimension");
        return StructIndexSetWrapper<shape_dim_, cell_dim_, face_dim__>::_index_set;
      }

      std::size_t bytes() const
      {
        return BaseClass::bytes() + _index_set.bytes();
      }
    };

    template<
      int shape_dim_,
      int cell_dim_>
    class StructIndexSetWrapper<shape_dim_, cell_dim_, 0>
    {
    protected:
      StructIndexSet<shape_dim_, cell_dim_, 0> _index_set;

    public:
      explicit StructIndexSetWrapper(const Index num_slices[]) :
        _index_set(num_slices)
      {
      }

      StructIndexSetWrapper(StructIndexSetWrapper&& other) :
        _index_set(std::forward<StructIndexSet<shape_dim_, cell_dim_, 0>>(other._index_set))
      {
      }

      StructIndexSetWrapper& operator=(StructIndexSetWrapper&& other)
      {
        if(this == &other)
          return *this;

        _index_set = std::forward<StructIndexSet<shape_dim_, cell_dim_, 0>>(other._index_set);

        return *this;
      }

      template<int face_dim__>
      const StructIndexSet<shape_dim_, cell_dim_, face_dim__>& get_index_set() const
      {
        static_assert(face_dim__ == 0, "invalid face dimension");
        return _index_set;
      }

      std::size_t bytes() const
      {
        return _index_set.bytes();
      }
    };

    template<
      int shape_dim_,
      int cell_dim_ = shape_dim_>
    class StructIndexSetHolder :
      public StructIndexSetHolder<shape_dim_, cell_dim_ - 1>
    {
    protected:
      typedef StructIndexSetHolder<shape_dim_, cell_dim_ - 1> BaseClass;

      StructIndexSetWrapper<shape_dim_, cell_dim_> _index_set_wrapper;

    public:
      explicit StructIndexSetHolder(const Index num_slices[]) :
        BaseClass(num_slices),
        _index_set_wrapper(num_slices)
      {
      }

      StructIndexSetHolder(StructIndexSetHolder&& other) :
        BaseClass(std::forward<BaseClass>(other)),
        _index_set_wrapper(std::forward<StructIndexSetWrapper<shape_dim_, cell_dim_>>(other._index_set_wrapper))
      {
      }

      StructIndexSetHolder& operator=(StructIndexSetHolder&& other)
      {
        if(this == &other)
          return *this;

        BaseClass::operator=(std::forward<BaseClass>(other));
        _index_set_wrapper = std::forward<StructIndexSetWrapper<shape_dim_, cell_dim_>>(other._index_set_wrapper);

        return *this;
      }

      template<int cell_dim__>
      const StructIndexSetWrapper<shape_dim_, cell_dim__>& get_index_set_wrapper() const
      {
        return StructIndexSetHolder<shape_dim_, cell_dim__>::_index_set_wrapper;
      }

      template<int cell_dim__, int face_dim__>
      const StructIndexSet<shape_dim_, cell_dim__, face_dim__>& get_index_set() const
      {
        return get_index_set_wrapper<cell_dim__>().template get_index_set<face_dim__>();
      }

      std::size_t bytes() const
      {
        return BaseClass::bytes() + _index_set_wrapper.bytes();
      }
    };

    template<int shape_dim_>
    class StructIndexSetHolder<shape_dim_, 0>
    {
    public:
      explicit StructIndexSetHolder(const Index /*num_slices*/[])
      {
      }

      StructIndexSetHolder(StructIndexSetHolder&&)
      {
      }

      StructIndexSetHolder& operator=(StructIndexSetHolder&&)
      {
        return *this;
      }

      std::size_t bytes() const
      {
        return std::size_t(0);
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_STRUCT_INDEX_SET_HPP

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_TARGET_SET_HPP
#define KERNEL_GEOMETRY_TARGET_SET_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Target set class
     *
     * A MeshPart refers to its parent mesh through several TargetSets.
     *
     * \author Peter Zajac
     */
    class TargetSet
    {
    protected:
      /// Index vector. _indices[i] = j means that entity i represents entity j in the parent
      std::vector<Index> _indices;

    public:
      /// standard constructor
      TargetSet() :
        _indices()
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] num_entities
       * The number of entities that are to be indexed.
       */
      explicit TargetSet(Index num_entities) :
        _indices(std::size_t(num_entities))
      {
      }

      /// move constructor
      TargetSet(TargetSet&& other) :
        _indices(std::forward<std::vector<Index>>(other._indices))
      {
      }

      /// move-assignment operator
      TargetSet& operator=(TargetSet&& other)
      {
        // avoid self-move
        if(this == &other)
          return *this;

        _indices = std::forward<std::vector<Index>>(other._indices);

        return *this;
      }

      /// virtual destructor
      virtual ~TargetSet()
      {
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return _indices.size() * sizeof(Index);
      }

      /// Returns the number of entities.
      Index get_num_entities() const
      {
        return Index(_indices.size());
      }

      /**
       * \returns A pointer to the target index array.
       */
      Index* get_indices()
      {
        return _indices.data();
      }

      /** \copydoc get_indices() */
      const Index* get_indices() const
      {
        return _indices.data();
      }

      /**
       * \brief Returns a target index.
       *
       * \param[in] i
       * The index of the entity whose target index is to be returned.
       *
       * \returns
       * A reference to the target index of entity \p i.
       */
      Index& operator[](Index i)
      {
        ASSERT(i < get_num_entities());
        return _indices[i];
      }

      /** \copydoc operator[]() */
      const Index& operator[](Index i) const
      {
        ASSERT(i < get_num_entities());
        return _indices[i];
      }
    }; // class TargetSet

    /// \cond internal
    /**
     * \brief Template recursive array of TargetSets
     *
     * A MeshPart referring to a mesh of Shape_ can have Shape_::dimension+1 TargetSets and this class provides the
     * means of accessing them. It inherits from the TargetSetHolder class wrt. the same shape of one dimension less.
     *
     * \tparam Shape_
     * Shape type this class refers to.
     *
     */
    template<typename Shape_>
    class TargetSetHolder :
      public TargetSetHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType>
    {
    public:
      typedef Shape_ ShapeType;
      static constexpr int shape_dim = ShapeType::dimension;

    protected:
      typedef TargetSetHolder<typename Shape::FaceTraits<ShapeType, shape_dim - 1>::ShapeType> BaseClass;

      TargetSet _target_set;

    public:
      TargetSetHolder() :
        BaseClass(),
        _target_set()
      {
      }

      explicit TargetSetHolder(const Index num_entities[]) :
        BaseClass(num_entities),
        _target_set(num_entities[shape_dim])
      {
      }

      TargetSetHolder(TargetSetHolder&& other) :
        BaseClass(std::forward<BaseClass>(other)),
        _target_set(std::forward<TargetSet>(other._target_set))
      {
      }

      TargetSetHolder& operator=(TargetSetHolder&& other)
      {
        if(this == &other)
          return *this;
        BaseClass::operator=(std::forward<BaseClass>(other));
        _target_set = std::forward<TargetSet>(other._target_set);
        return *this;
      }

      virtual ~TargetSetHolder()
      {
      }

      std::size_t bytes() const
      {
        return BaseClass::bytes() + _target_set.bytes();
      }

      template<int dim_>
      TargetSet& get_target_set()
      {
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return TargetSetHolder<CellType>::_target_set;
      }

      template<int dim_>
      const TargetSet& get_target_set() const
      {
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return TargetSetHolder<CellType>::_target_set;
      }

      Index get_num_entities(int dim) const
      {
        XASSERTM(dim <= shape_dim, "invalid dimension parameter");
        if(dim == shape_dim)
        {
          return _target_set.get_num_entities();
        }
        return BaseClass::get_num_entities(dim);
      }

      static String name()
      {
        return "TargetSetHolder<" + Shape_::name() + ">";
      }
    };

    /**
     * \brief TargetSetHolder for lowest possible dimension as end of the recursive inheritance
     */
    template<>
    class TargetSetHolder<Shape::Vertex>
    {
    public:
      typedef Shape::Vertex ShapeType;
      static constexpr int shape_dim = ShapeType::dimension;

    protected:
      TargetSet _target_set;

    public:
      TargetSetHolder() :
        _target_set()
      {
      }

      explicit TargetSetHolder(const Index num_entities[]) :
        _target_set(num_entities[0])
      {
      }

      TargetSetHolder(TargetSetHolder&& other) :
        _target_set(std::forward<TargetSet>(other._target_set))
      {
      }

      TargetSetHolder& operator=(TargetSetHolder&& other)
      {
        if(this == &other)
          return *this;
        _target_set = std::forward<TargetSet>(other._target_set);
        return *this;
      }

      virtual ~TargetSetHolder()
      {
      }

      std::size_t bytes() const
      {
        return _target_set.bytes();
      }

      template<int dim_>
      TargetSet& get_target_set()
      {
        static_assert(dim_ == 0, "invalid dimension");
        return _target_set;
      }

      template<int dim_>
      const TargetSet& get_target_set() const
      {
        static_assert(dim_ == 0, "invalid dimension");
        return _target_set;
      }

      Index get_num_entities(int dim) const
      {
        XASSERTM(dim == 0, "invalid dimension parameter");
        return _target_set.get_num_entities();
      }

      static String name()
      {
        return "TargetSetHolder<Vertex>";
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TARGET_SET_HPP

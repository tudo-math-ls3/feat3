#pragma once
#ifndef KERNEL_GEOMETRY_TARGET_SET_HPP
#define KERNEL_GEOMETRY_TARGET_SET_HPP 1

// includes, FEAST
#include <kernel/geometry/shape.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Target set class
     *
     * \author Peter Zajac
     */
    class TargetSet
    {
    protected:
      /// number of entities
      Index _num_entities;

      /// index array
      Index* _indices;

    private:
      TargetSet(const TargetSet&);
      TargetSet& operator=(const TargetSet&);

    public:
      /**
       * \brief Constructor
       *
       * \param[in] num_entities
       * The number of entities that are to be indexed.
       */
      explicit TargetSet(Index num_entities = 0) :
        _num_entities(num_entities),
        _indices(nullptr)
      {
        CONTEXT("TargetSet::TargetSet(Index)");
        if(num_entities > 0)
        {
          _indices = new Index[num_entities];
        }
      }

      /// virtual destructor
      virtual ~TargetSet()
      {
        CONTEXT("TargetSet::~TargetSet()");
        if(_indices != nullptr)
        {
          delete [] _indices;
        }
      }

      /// Returns the number of entities.
      Index get_num_entities() const
      {
        return _num_entities;
      }

      /// Returns the target index array.
      Index* get_indices()
      {
        return _indices;
      }

      /** \copydoc get_indices() */
      const Index* get_indices() const
      {
        return _indices;
      }
      /*
      Index& get_index(Index i)
      {
        return _indices[i];
      }

      const Index& get_index(Index i) const
      {
        return _indices[i];
      }*/

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
        return _indices[i];
      }

      /** \copydoc operator[]() */
      const Index& operator[](Index i) const
      {
        return _indices[i];
      }
    }; // class TargetSet

    /* ***************************************************************************************** */
    /// \cond internal
    template<typename Shape_>
    class TargetSetHolder :
      public TargetSetHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType>
    {
    public:
      typedef Shape_ ShapeType;
      enum
      {
        shape_dim = ShapeType::dimension
      };

    protected:
      typedef TargetSetHolder<typename Shape::FaceTraits<ShapeType, shape_dim - 1>::ShapeType> BaseClass;

      TargetSet _target_set;

    public:
      explicit TargetSetHolder(const Index num_entities[]) :
        BaseClass(num_entities),
        _target_set(num_entities[shape_dim])
      {
        CONTEXT(name() + "::TargetSetHolder(const Index[])");
      }

      virtual ~TargetSetHolder()
      {
        CONTEXT(name() + "::~TargetSetHolder()");
      }

      template<int dim_>
      TargetSet& get_target_set()
      {
        CONTEXT(name() + "::get_target_set()");
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return TargetSetHolder<CellType>::_target_set;
      }

      template<int dim_>
      const TargetSet& get_target_set() const
      {
        CONTEXT(name() + "::get_target_set() [const]");
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return TargetSetHolder<CellType>::_target_set;
      }

      Index get_num_entities(int dim) const
      {
        CONTEXT(name() + "::get_num_entities()");
        ASSERT(dim <= shape_dim, "invalid dimension parameter");
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

    template<>
    class TargetSetHolder<Shape::Vertex>
    {
    public:
      typedef Shape::Vertex ShapeType;
      enum
      {
        shape_dim = ShapeType::dimension
      };

    protected:
      TargetSet _target_set;

    public:
      TargetSetHolder() :
        _target_set()
      {
        CONTEXT(name() + "::TargetSetHolder()");
      }

      explicit TargetSetHolder(const Index num_entities[]) :
        _target_set(num_entities[0])
      {
        CONTEXT(name() + "::TargetSetHolder(const Index[])");
      }

      virtual ~TargetSetHolder()
      {
        CONTEXT(name() + "::~TargetSetHolder()");
      }

      template<int dim_>
      TargetSet& get_target_set()
      {
        CONTEXT(name() + "::get_target_set()");
        static_assert(dim_ == 0, "invalid dimension");
        return _target_set;
      }

      template<int dim_>
      const TargetSet& get_target_set() const
      {
        CONTEXT(name() + "::get_target_set() [const]");
        static_assert(dim_ == 0, "invalid dimension");
        return _target_set;
      }

      Index get_num_entities(int dim) const
      {
        CONTEXT(name() + "::get_num_entities()");
        ASSERT(dim == 0, "invalid dimension parameter");
        return _target_set.get_num_entities();
      }

      static String name()
      {
        return "TargetSetHolder<Vertex>";
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TARGET_SET_HPP

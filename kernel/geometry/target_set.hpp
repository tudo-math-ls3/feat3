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
    template<int dim_>
    class TargetSetHolder :
      public TargetSetHolder<dim_ - 1>
    {
    protected:
      typedef TargetSetHolder<dim_ - 1> BaseClass;

      TargetSet _target_set;

    public:
      explicit TargetSetHolder(const Index num_entities[]) :
        BaseClass(num_entities),
        _target_set(num_entities[dim_])
      {
        CONTEXT(name() + "::TargetSetHolder(const Index[])");
      }

      virtual ~TargetSetHolder()
      {
        CONTEXT(name() + "::~TargetSetHolder()");
      }

      template<int dim__>
      TargetSet& get_target_set()
      {
        CONTEXT(name() + "::get_target_set()");
        static_assert(dim__ >= 0, "invalid dimension");
        static_assert(dim__ <= dim_, "invalid dimension");
        return TargetSetHolder<dim__>::_target_set;
      }

      template<int dim__>
      const TargetSet& get_target_set() const
      {
        CONTEXT(name() + "::get_target_set() [const]");
        static_assert(dim__ >= 0, "invalid dimension");
        static_assert(dim__ <= dim_, "invalid dimension");
        return TargetSetHolder<dim__>::_target_set;
      }

      Index get_num_entities(int dim) const
      {
        CONTEXT(name() + "::get_num_entities()");
        ASSERT(dim <= dim_, "invalid dimension parameter");
        if(dim == dim_)
        {
          return _target_set.get_num_entities();
        }
        return BaseClass::get_num_entities(dim);
      }

      static String name()
      {
        return "TargetSetHolder<" + stringify(dim_) + ">";
      }
    };

    template<>
    class TargetSetHolder<0>
    {
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

      template<int dim__>
      TargetSet& get_target_set()
      {
        CONTEXT(name() + "::get_target_set()");
        static_assert(dim__ == 0, "invalid dimension");
        return _target_set;
      }

      template<int dim__>
      const TargetSet& get_target_set() const
      {
        CONTEXT(name() + "::get_target_set() [const]");
        static_assert(dim__ == 0, "invalid dimension");
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
        return "TargetSetHolder<0>";
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TARGET_SET_HPP

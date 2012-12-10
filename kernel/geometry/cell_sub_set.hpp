#pragma once
#ifndef KERNEL_GEOMETRY_CELL_SUB_SET_HPP
#define KERNEL_GEOMETRY_CELL_SUB_SET_HPP 1

// includes, FEAST
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/intern/standard_subset_refiner.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Cell Sub-Set class template
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class CellSubSet
    {
    public:
      /// shape type
      typedef Shape_ ShapeType;
      /// base class typedef
      typedef TargetSetHolder<ShapeType> BaseClass;
      /// target set holder type
      typedef TargetSetHolder<ShapeType> TargetSetHolderType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = ShapeType::dimension
      };

      /**
       * \brief Target set type class template.
       *
       * This nested class template is used to define the return type of the CellSubSet::get_target_set()
       * function template.
       *
       * \tparam cell_dim_
       * The cell dimesnion parameter as passed to the CellSubSet::get_target_set() function template.
       */
      template<int cell_dim_>
      struct TargetSet
      {
        /// target set type
        typedef FEAST::Geometry::TargetSet Type;
      }; // struct TargetSet<...>

    protected:
      /// the target sets of the cell subset
      TargetSetHolderType _target_set_holder;

    private:
      CellSubSet(const CellSubSet&);
      CellSubSet& operator=(const CellSubSet&);

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_entities
       * An array of length at least #shape_dim + 1 holding the number of entities for each shape dimension.
       * Must not be \c nullptr.
       */
      explicit CellSubSet(const Index num_entities[]) :
        _target_set_holder(num_entities)
      {
        CONTEXT(name() + "::CellSubSet()");
      }

      /**
       * \brief Factory constructor
       *
       * \param[in] factory
       * The factory that is to be used to create the cell subset.
       */
      explicit CellSubSet(Factory<CellSubSet>& factory) :
        _target_set_holder(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities)
      {
        CONTEXT(name() + "::CellSubSet() [factory]");
        factory.fill_target_sets(_target_set_holder);
      }

      /// virtual destructor
      virtual ~CellSubSet()
      {
        CONTEXT(name() + "::~CellSubSet()");
      }

      /**
       * \brief Returns the number of entities.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= #shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      Index get_num_entities(int dim) const
      {
        CONTEXT(name() + "::get_num_entities()");
        ASSERT_(dim >= 0);
        ASSERT_(dim <= shape_dim);
        return _target_set_holder.get_num_entities(dim);
      }

      /**
       * \brief Return the reference to a target set.
       *
       * \tparam cell_dim_
       * The dimension fo the entity whose target set is to be returned.
       *
       * \returns
       * A reference to the target set.
       */
      template<int cell_dim_>
      typename TargetSet<cell_dim_>::Type& get_target_set()
      {
        CONTEXT(name() + "::get_target_set<" + stringify(cell_dim_) + ">()");
        static_assert(cell_dim_ >= 0, "invalid cell dimension");
        static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
        return _target_set_holder.template get_target_set<cell_dim_>();
      }

      /** \copydoc get_target_set() */
      template<int cell_dim_>
      const typename TargetSet<cell_dim_>::Type& get_target_set() const
      {
        CONTEXT(name() + "::get_target_set<" + stringify(cell_dim_) + ">() [const]");
        static_assert(cell_dim_ >= 0, "invalid cell dimension");
        static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
        return _target_set_holder.template get_target_set<cell_dim_>();
      }

      /// \cond internal
      TargetSetHolderType& get_target_set_holder()
      {
        CONTEXT(name() + "::get_target_set_holder()");
        return _target_set_holder;
      }

      const TargetSetHolderType& get_target_set_holder() const
      {
        CONTEXT(name() + "::get_target_set_holder() [const]");
        return _target_set_holder;
      }
      /// \endcond

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "CellSubSet<" + ShapeType::name() + ">";
      }
    }; // class CellSubSet<...>

    /* ************************************************************************************************************* */

    /**
     * \brief Factory specialisation for CellSubSet class template.
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class Factory< CellSubSet<Shape_> >
    {
    public:
      /// mesh type
      typedef CellSubSet<Shape_> MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    public:
      /// virtual destructor
      virtual ~Factory()
      {
      }

      /**
       * \brief Returns the number of entities.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      virtual Index get_num_entities(int dim) = 0;

      /**
       * \brief Fills the target sets.
       *
       * \param[in,out] target_set_holder
       * The target set holder whose target sets are to be filled.
       */
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) = 0;
    }; // class Factory<CellSubSet<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief StandardRefinery implementation for CellSubSet
     *
     * \author Peter Zajac
     */
    template<typename Shape_, typename Parent_>
    class StandardRefinery<CellSubSet<Shape_>, Parent_> :
      public Factory< CellSubSet<Shape_> >
    {
    public:
      /// shape type
      typedef Shape_ ShapeType;
      /// mesh type
      typedef CellSubSet<Shape_> MeshType;
      /// parent type
      typedef Parent_ ParentType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = ShapeType::dimension
      };

    protected:
      /// coarse cell subset reference
      const MeshType& _coarse_mesh;
      /// coarse parent reference
      const ParentType& _parent;
      /// number of entities in refined cell subset
      Index _num_entities_fine[shape_dim + 1];
      /// number of entities in coarse parent
      Index _num_entities_parent[shape_dim + 1];

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] coarse_mesh
       * A reference to the coarse cell subset that is to be refined.
       *
       * \param[in] parent
       * A reference to the coarse parent.
       */
      explicit StandardRefinery(const MeshType& coarse_mesh, const ParentType& parent) :
        _coarse_mesh(coarse_mesh),
        _parent(parent)
      {
        // get number of entities in coarse mesh
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities_fine[i] = coarse_mesh.get_num_entities(i);
          _num_entities_parent[i] = parent.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        Intern::EntityCountWrapper<ShapeType>::query(_num_entities_fine);
      }

      /// virtual destructor
      virtual ~StandardRefinery()
      {
      }

      /**
       * \brief Returns the number of entities of the refined cellset.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= #shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      virtual Index get_num_entities(int dim)
      {
        return _num_entities_fine[dim];
      }

      /**
       * \brief Fills the target sets of the refined cell subset.
       *
       * \param[in,out] target_set_holder
       * The target set holder whose target sets are to be filled.
       */
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder)
      {
        // refine subset target indices
        Intern::SubSetRefineWrapper<ShapeType>
          ::refine(target_set_holder, _num_entities_parent, _coarse_mesh.get_target_set_holder());
      }
    };
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CELL_SUB_SET_HPP

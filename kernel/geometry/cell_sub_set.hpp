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
    class CellSubSet :
      public TargetSetHolder<Shape_>
    {
      // friends
      template<typename Mesh_, typename Parent_>
      friend class StandardRefinery;

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

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_entities
       * An array of length at least #shape_dim + 1 holding the number of entities for each shape dimension.
       * Must not be \c nullptr.
       */
      explicit CellSubSet(const Index num_entities[]) :
        BaseClass(num_entities)
      {
        CONTEXT(name() + "::CellSubSet()");
      }

      explicit CellSubSet(const Factory<CellSubSet>& factory) :
        BaseClass(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities)
      {
        factory.fill_target_sets(*this);
      }

      /// virtual destructor
      virtual ~CellSubSet()
      {
        CONTEXT(name() + "::~CellSubSet()");
      }

      /**
       * \brief Refines the cell subset.
       *
       * This function applies the standard refinement algorithm onto the cell subset and returns the
       * refined cell subset.
       *
       * \param[in] parent
       * A reference to the (coarse) parent mesh or cell subset that this cell subset refers to.
       *
       * \returns
       * A pointer to the refined cell subset.
       */
      template<typename Parent_>
      CellSubSet* refine(const Parent_& parent) const
      {
        CONTEXT(name() + "::refine()");

        return new CellSubSet(StandardRefinery<CellSubSet, Parent_>(*this, parent));
      }

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
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= #shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      virtual Index get_num_entities(int dim) const = 0;

      /**
       * \brief Fills the target sets.
       *
       * \param[in,out] target_set_holder
       * The target set holder whose target sets are to be filled.
       */
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) const = 0;
    }; // class Factory<CellSubSet<...>>

    /* ************************************************************************************************************* */

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
      virtual Index get_num_entities(int dim) const
      {
        return _num_entities_fine[dim];
      }

      /**
       * \brief Fills the target sets of the refined cell subset.
       *
       * \param[in,out] target_set_holder
       * The target set holder whose target sets are to be filled.
       */
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) const
      {
        // refine subset target indices
        Intern::SubSetRefineWrapper<ShapeType>::refine(target_set_holder, _num_entities_parent, _coarse_mesh);
      }
    };
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CELL_SUB_SET_HPP

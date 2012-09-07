#pragma once
#ifndef KERNEL_GEOMETRY_CELL_SUB_SET_HPP
#define KERNEL_GEOMETRY_CELL_SUB_SET_HPP 1

// includes, FEAST
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
    public:
      /// shape type
      typedef Shape_ ShapeType;
      /// base class typedef
      typedef TargetSetHolder<ShapeType> BaseClass;

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

        // get number of entities in coarse mesh
        Index num_entities_fine[BaseClass::shape_dim + 1];
        Index num_entities_parent[BaseClass::shape_dim + 1];
        for(int i(0); i <= BaseClass::shape_dim; ++i)
        {
          num_entities_fine[i] = BaseClass::get_num_entities(i);
          num_entities_parent[i] = parent.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        Intern::EntityCountWrapper<ShapeType>::query(num_entities_fine);

        // allocate a fine subset
        CellSubSet* fine_set = new CellSubSet(num_entities_fine);

        // refine subset target indices
        Intern::SubSetRefineWrapper<ShapeType>::refine(*fine_set, num_entities_parent, *this);

        // return fine subset
        return fine_set;
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
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CELL_SUB_SET_HPP

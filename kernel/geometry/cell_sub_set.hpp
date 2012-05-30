#pragma once
#ifndef KERNEL_GEOMETRY_CELL_SUB_SET_HPP
#define KERNEL_GEOMETRY_CELL_SUB_SET_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal/standard_refinement/subset_refiner.hpp>

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
      typedef Shape_ ShapeType;
      typedef TargetSetHolder<ShapeType> BaseClass;

    public:
      explicit CellSubSet(const Index num_entities[]) :
        BaseClass(num_entities)
      {
        CONTEXT(name() + "::CellSubSet()");
      }

      virtual ~CellSubSet()
      {
        CONTEXT(name() + "::~CellSubSet()");
      }

      template<typename ParentMesh_>
      CellSubSet* refine(const ParentMesh_& parent_mesh) const
      {
        CONTEXT(name() + "::refine()");
        using namespace Conformal::StandardRefinement;

        // get number of entities in coarse mesh
        Index num_entities_fine[BaseClass::shape_dim + 1];
        Index num_entities_parent[BaseClass::shape_dim + 1];
        for(int i(0); i <= BaseClass::shape_dim; ++i)
        {
          num_entities_fine[i] = BaseClass::get_num_entities(i);
          num_entities_parent[i] = parent_mesh.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        EntityCountWrapper<ShapeType>::query(num_entities_fine);

        // allocate a fine subset
        CellSubSet* fine_set = new CellSubSet(num_entities_fine);

        // refine subset target indices
        SubSetRefineWrapper<ShapeType>::refine(*fine_set, num_entities_parent, *this);

        // return fine subset
        return fine_set;
      }

      static String name()
      {
        return "CellSubSet<" + ShapeType::name() + ">";
      }
    }; // class CellSubSet<...>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CELL_SUB_SET_HPP

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_SIMPLE_TARGET_REFINER_HPP
#define KERNEL_GEOMETRY_INTERN_SIMPLE_TARGET_REFINER_HPP 1

// includes, FEAT
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/intern/entity_counter.hpp>
#include <kernel/geometry/intern/standard_refinement_traits.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        int cell_dim_>
      struct SimpleTargetRefiner
      {
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index offset,
          const Index index_offsets[],
          const TargetSetType& target_set_in)
        {
          // fetch number of cells
          Index num_cells = target_set_in.get_num_entities();
          static constexpr int shape_dim = ShapeType::dimension;
          static constexpr int num_children = StandardRefinementTraits<ShapeType,cell_dim_>::count;

          for(Index i(0); i < num_cells; ++i)
          {
            for(int j(0); j < num_children; ++j)
            {
              target_set_out[offset + i*Index(num_children) + Index(j)] =
                index_offsets[shape_dim] + target_set_in[i]*Index(num_children) + Index(j);
            }
          }

          return num_cells*Index(num_children);
        }
      };

      template<
        typename Shape_,
        int cell_dim_,
        int shape_dim_ = Shape_::dimension>
      struct SimpleTargetRefineShapeWrapper
      {
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index index_offsets[],
          const TargetSetHolderType& target_set_holder_in)
        {
          // recursive call of SubSetRefineShapeWrapper
          typedef typename Shape::FaceTraits<ShapeType, shape_dim_ - 1>::ShapeType FacetType;
          Index offset = SimpleTargetRefineShapeWrapper<FacetType, cell_dim_>::refine(
            target_set_out,
            index_offsets,
            target_set_holder_in);

          // get input target set
          const TargetSetType& target_set_in = target_set_holder_in.template get_target_set<shape_dim_>();

          // call subset refiner
          SimpleTargetRefiner<ShapeType, cell_dim_>::refine(
            target_set_out,
            offset,
            index_offsets,
            target_set_in);

          // return new offset
          return offset + Intern::StandardRefinementTraits<ShapeType, cell_dim_>::count *
            target_set_holder_in.get_num_entities(shape_dim_);
        }
      };

      template<
        typename Shape_,
        int cell_dim_>
      struct SimpleTargetRefineShapeWrapper<Shape_, cell_dim_, cell_dim_>
      {
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
        typedef TargetSetHolder<CellType> TargetSetHolderType;

        static Index refine(
          TargetSetType& target_set_out,
          const Index index_offsets[],
          const TargetSetHolderType& target_set_holder_in)
        {
          // get input target set
          const TargetSetType& target_set_in = target_set_holder_in.template get_target_set<cell_dim_>();

          // call subset refiner
          SimpleTargetRefiner<ShapeType, cell_dim_>::refine(
            target_set_out,
            0,
            index_offsets,
            target_set_in);

          // return new offset
          return StandardRefinementTraits<ShapeType, cell_dim_>::count *
            target_set_holder_in.get_num_entities(cell_dim_);
        }
      };

      template<
        typename Shape_,
        int cell_dim_ = Shape_::dimension>
      struct SimpleTargetRefineWrapper
      {
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
        typedef TargetSetHolder<CellType> TargetSetHolderTypeOut;
        typedef TargetSetHolder<ShapeType> TargetSetHolderTypeIn;

        static void refine(
          TargetSetHolderTypeOut& target_set_holder_out,
          const Index num_entities_trg[],
          const TargetSetHolderTypeIn& target_set_holder_in)
        {
          // recursive call of SubSetRefineWrapper
          SimpleTargetRefineWrapper<ShapeType, cell_dim_ - 1>::refine(
            target_set_holder_out,
            num_entities_trg,
            target_set_holder_in);

          // calculate index offsets
          Index index_offsets[Shape_::dimension + 1];
          EntityCounter<StandardRefinementTraits, ShapeType, cell_dim_>::offset(index_offsets, num_entities_trg);

          // get output target set
          TargetSet& target_set_out = target_set_holder_out.template get_target_set<cell_dim_>();

          // call shape wrapper
          SimpleTargetRefineShapeWrapper<ShapeType, cell_dim_>::refine(
            target_set_out,
            index_offsets,
            target_set_holder_in);
        }
      };

      template<typename Shape_>
      struct SimpleTargetRefineWrapper<Shape_, 0>
      {
        typedef Shape_ ShapeType;
        typedef TargetSet TargetSetType;
        typedef typename Shape::FaceTraits<Shape_, 0>::ShapeType CellType;
        typedef TargetSetHolder<CellType> TargetSetHolderTypeOut;
        typedef TargetSetHolder<ShapeType> TargetSetHolderTypeIn;

        static void refine(
          TargetSetHolderTypeOut& target_set_holder_out,
          const Index num_entities_trg[],
          const TargetSetHolderTypeIn& target_set_holder_in)
        {
          // calculate index offsets
          Index index_offsets[Shape_::dimension + 1];
          EntityCounter<StandardRefinementTraits, ShapeType, 0>::offset(index_offsets, num_entities_trg);

          // get output target set
          TargetSet& target_set_out = target_set_holder_out.template get_target_set<0>();

          // call shape wrapper
          SimpleTargetRefineShapeWrapper<ShapeType, 0>::refine(
            target_set_out,
            index_offsets,
            target_set_holder_in);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_SIMPLE_TARGET_REFINER_HPP

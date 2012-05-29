#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_WRAPPERS_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_WRAPPERS_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal/standard_refinement/entity_counter.hpp>
#include <kernel/geometry/conformal/standard_refinement/target_refiner-hypercube.hpp>
#include <kernel/geometry/conformal/standard_refinement/target_refiner-simplex.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Conformal
    {
      /// \cond internal
      namespace StandardRefinement
      {
        template<
          typename Shape_,
          int cell_dim_,
          int shape_dim_ = Shape_::dimension>
        class TargetRefineShapeWrapper
        {
        public:
          typedef Shape_ ShapeType;
          typedef TargetSet TargetSetType;
          typedef TargetSetHolder<ShapeType> TargetSetHolderType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "TargetRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_)
              + "," + stringify(shape_dim_) + ">";
          }

          static Index refine(
            TargetSetType& target_set_out,
            const Index index_offsets[],
            const TargetSetHolderType& target_set_holder_in,
            const IndexSetHolderType& index_set_holder_src,
            const IndexSetHolderType& index_set_holder_trg)
          {
            CONTEXT(name() + "::refine()");

            // recursive call of TargetRefineShapeWrapper
            typedef typename Shape::FaceTraits<ShapeType, shape_dim_ - 1>::ShapeType FacetType;
            Index offset = TargetRefineShapeWrapper<FacetType, cell_dim_>::refine(
              target_set_out,
              index_offsets,
              target_set_holder_in,
              index_set_holder_src,
              index_set_holder_trg);

            // call target refiner
            TargetRefiner<ShapeType, cell_dim_>::refine(
              target_set_out,
              offset,
              index_offsets,
              target_set_holder_in,
              index_set_holder_src,
              index_set_holder_trg);

            // return new offset
            return offset + StdRefTraits<ShapeType, cell_dim_>::count *
              target_set_holder_in.get_num_entities(shape_dim_);
          }
        };

        template<
          typename Shape_,
          int cell_dim_>
        class TargetRefineShapeWrapper<Shape_, cell_dim_, cell_dim_>
        {
        public:
          typedef Shape_ ShapeType;
          typedef TargetSet TargetSetType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef TargetSetHolder<CellType> TargetSetHolderType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "TargetRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_)
              + "," + stringify(cell_dim_) + ">";
          }

          static Index refine(
            TargetSetType& target_set_out,
            const Index index_offsets[],
            const TargetSetHolderType& target_set_holder_in,
            const IndexSetHolderType& index_set_holder_src,
            const IndexSetHolderType& index_set_holder_trg)
          {
            CONTEXT(name() + "::refine()");

            // call target refiner
            TargetRefiner<ShapeType, cell_dim_>::refine(
              target_set_out,
              0,
              index_offsets,
              target_set_holder_in,
              index_set_holder_src,
              index_set_holder_trg);

            // return new offset
            return StdRefTraits<ShapeType, cell_dim_>::count * target_set_holder_in.get_num_entities(cell_dim_);
          }
        };

        template<int cell_dim_>
        class TargetRefineShapeWrapper<Shape::Vertex, cell_dim_, cell_dim_>
        {
        public:
          typedef Shape::Vertex ShapeType;
          typedef TargetSet TargetSetType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef TargetSetHolder<CellType> TargetSetHolderType;
          typedef IndexSetHolder<Shape::Vertex> IndexSetHolderType;

          static String name()
          {
            return "TargetRefineShapeWrapper<Vertex," + stringify(cell_dim_) + "," + stringify(cell_dim_) + ">";
          }

          static Index refine(
            TargetSetType& target_set_out,
            const Index index_offsets[],
            const TargetSetHolderType& target_set_holder_in,
            const IndexSetHolderType& /*index_set_holder_src*/,
            const IndexSetHolderType& /*index_set_holder_trg*/)
          {
            // call target refiner
            TargetRefiner<ShapeType, cell_dim_>::refine(target_set_out, 0, index_offsets, target_set_holder_in);

            // return new offset
            return StdRefTraits<ShapeType, cell_dim_>::count * target_set_holder_in.get_num_entities(cell_dim_);
          }
        };

        /* ********************************************************************************************************* */

        template<
          typename Shape_,
          int cell_dim_ = Shape_::dimension>
        class TargetRefineWrapper
        {
        public:
          typedef Shape_ ShapeType;
          typedef TargetSet TargetSetType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef TargetSetHolder<CellType> TargetSetHolderTypeOut;
          typedef TargetSetHolder<ShapeType> TargetSetHolderTypeIn;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "TargetRefineWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ">";
          }

          static void refine(
            TargetSetHolderTypeOut& target_set_holder_out,
            const Index num_entities_trg[],
            const TargetSetHolderTypeIn& target_set_holder_in,
            const IndexSetHolderType& index_set_holder_src,
            const IndexSetHolderType& index_set_holder_trg)
          {
            CONTEXT(name() + "::refine()");

            // recursive call of TargetRefineWrapper
            TargetRefineWrapper<ShapeType, cell_dim_ - 1>::refine(
              target_set_holder_out,
              num_entities_trg,
              target_set_holder_in,
              index_set_holder_src,
              index_set_holder_trg);

            // calculate index offsets
            Index index_offsets[Shape_::dimension+1];
            EntityCounter<ShapeType, cell_dim_>::offset(index_offsets, num_entities_trg);

            // call shape wrapper
            TargetRefineShapeWrapper<ShapeType, cell_dim_>::refine(
              target_set_holder_out.template get_target_set<cell_dim_>(),
              index_offsets,
              target_set_holder_in,
              index_set_holder_src,
              index_set_holder_trg);
          }
        };

        template<typename Shape_>
        class TargetRefineWrapper<Shape_, 0>
        {
        public:
          typedef Shape_ ShapeType;
          typedef TargetSet TargetSetType;
          typedef TargetSetHolder<Shape::Vertex> TargetSetHolderTypeOut;
          typedef TargetSetHolder<ShapeType> TargetSetHolderTypeIn;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "TargetRefineWrapper<" + Shape_::name() + ",0>";
          }

          static void refine(
            TargetSetHolderTypeOut& target_set_holder_out,
            const Index num_entities_trg[],
            const TargetSetHolderTypeIn& target_set_holder_in,
            const IndexSetHolderType& index_set_holder_src,
            const IndexSetHolderType& index_set_holder_trg)
          {
            CONTEXT(name() + "::refine()");

            // calculate index offsets
            Index index_offsets[Shape_::dimension+1];
            EntityCounter<ShapeType, 0>::offset(index_offsets, num_entities_trg);

            // call shape wrapper
            TargetRefineShapeWrapper<ShapeType, 0>::refine(
              target_set_holder_out.template get_target_set<0>(),
              index_offsets,
              target_set_holder_in,
              index_set_holder_src,
              index_set_holder_trg);
          }
        };
      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_TRG_REF_WRAPPERS_HPP

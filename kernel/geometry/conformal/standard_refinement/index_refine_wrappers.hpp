#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_WRAPPERS_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_WRAPPERS_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal/standard_refinement/entity_counter.hpp>
#include <kernel/geometry/conformal/standard_refinement/index_refiner-hypercube.hpp>
#include <kernel/geometry/conformal/standard_refinement/index_refiner-simplex.hpp>

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
          int face_dim_,
          // The following "dummy" argument is necessary for partial specialisation;
          // shape_dim_ *must* coincide with Shape_::dimension !!!
          int shape_dim_ = Shape_::dimension>
        struct IndexRefineShapeWrapper
        {
          static_assert(shape_dim_ == Shape_::dimension, "invalid shape dimension");
          // the case shape_dim_ = cell_dim_ is specialised below
          static_assert(shape_dim_ > cell_dim_, "invalid cell dimension");
          static_assert(cell_dim_ > face_dim_, "invalid face dimension");
          static_assert(face_dim_ >= 0, "invalid face_dimension");

          typedef Shape_ ShapeType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef IndexSet<Shape::FaceTraits<CellType, face_dim_>::count> IndexSetType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "IndexRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
              + stringify(face_dim_) + "," + stringify(shape_dim_) + ">";
          }

          static Index refine(
            IndexSetType& index_set_out,
            const Index index_offsets[],
            const IndexSetHolderType& index_set_holder_in)
          {
            typedef typename Shape::FaceTraits<ShapeType, shape_dim_ - 1>::ShapeType FacetType;

            // refine facets
            Index offset = IndexRefineShapeWrapper<FacetType, cell_dim_, face_dim_>
              ::refine(index_set_out, index_offsets, index_set_holder_in);

            // call index refiner and return new offset
            Index num_faces = IndexRefiner<Shape_, cell_dim_, face_dim_>
              ::refine(index_set_out, offset, index_offsets, index_set_holder_in);

#ifdef DEBUG
            // validate number of created indices
            Index num_shapes = index_set_holder_in.template get_index_set<shape_dim_, 0>().get_num_entities();
            ASSERT(num_faces == (StdRefTraits<ShapeType, cell_dim_>::count * num_shapes),
              "IndexRefiner output does not match StdRefTraits prediction");
#endif // DEBUG

            // calculate new offset and return
            return offset + num_faces;
          }
        };

        template<
          typename Shape_,
          int cell_dim_,
          int face_dim_>
        struct IndexRefineShapeWrapper<Shape_, cell_dim_, face_dim_, cell_dim_>
        {
          static_assert(cell_dim_ == Shape_::dimension, "invalid shape dimension");
          static_assert(cell_dim_ > face_dim_, "invalid face dimension");
          static_assert(face_dim_ >= 0, "invalid face_dimension");

          typedef Shape_ ShapeType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef IndexSet<Shape::FaceTraits<CellType, face_dim_>::count> IndexSetType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "IndexRefineShapeWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
              + stringify(face_dim_) + "," + stringify(cell_dim_) + ">";
          }

          static Index refine(
            IndexSetType& index_set_out,
            const Index index_offsets[],
            const IndexSetHolderType& index_set_holder_in)
          {
            // call index refiner
            Index num_faces = IndexRefiner<Shape_, cell_dim_, face_dim_>
              ::refine(index_set_out, 0, index_offsets, index_set_holder_in);

#ifdef DEBUG
            // validate number of created indices
            Index num_shapes = index_set_holder_in.template get_index_set<cell_dim_, 0>().get_num_entities();
            ASSERT(num_faces == (StdRefTraits<ShapeType, cell_dim_>::count * num_shapes),
              "IndexRefiner output does not match StdRefTraits prediction");
#endif // DEBUG

            // return offset
            return num_faces;
          }
        };

        /* ********************************************************************************************************* */

        template<
          typename Shape_,
          int cell_dim_,
          int face_dim_ = cell_dim_ - 1>
        struct IndexRefineFaceWrapper
        {
          // the case face_dim_ = 0 is handled by the partial specialisation below
          static_assert(face_dim_ > 0, "invalid face dimension");
          static_assert(cell_dim_ > face_dim_, "invalid cell dimension");
          static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

          typedef Shape_ ShapeType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef IndexSet<Shape::FaceTraits<CellType, face_dim_>::count> IndexSetType;
          typedef IndexSetWrapper<CellType, face_dim_> IndexSetWrapperType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "IndexRefineFaceWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ","
              + stringify(face_dim_) + ">";
          }

          static void refine(
            IndexSetWrapperType& index_set_wrapper_out,
            const Index num_entities[],
            const IndexSetHolderType& index_set_holder_in)
          {
            CONTEXT(name() + "::refine()");

            // recursive call of IndexRefineFaceWrapper
            IndexRefineFaceWrapper<Shape_, cell_dim_, face_dim_ - 1>
              ::refine(index_set_wrapper_out, num_entities, index_set_holder_in);

            // fetch output index set
            IndexSetType& index_set_out = index_set_wrapper_out.template get_index_set<face_dim_>();

            // calculate index offsets
            Index index_offsets[Shape_::dimension+1];
            EntityCounter<Shape_, face_dim_>::offset(index_offsets, num_entities);

            // call refinement shape wrapper
            IndexRefineShapeWrapper<Shape_, cell_dim_, face_dim_>
              ::refine(index_set_out, index_offsets, index_set_holder_in);
          }
        };

        template<
          typename Shape_,
          int cell_dim_>
        struct IndexRefineFaceWrapper<Shape_, cell_dim_, 0>
        {
          static_assert(cell_dim_ > 0, "invalid cell dimension");
          static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

          typedef Shape_ ShapeType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef IndexSet<Shape::FaceTraits<CellType, 0>::count> IndexSetType;
          typedef IndexSetWrapper<CellType, 0> IndexSetWrapperType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "IndexRefineFaceWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ",0>";
          }

          static void refine(
            IndexSetWrapperType& index_set_wrapper_out,
            const Index num_entities[],
            const IndexSetHolderType& index_set_holder_in)
          {
            CONTEXT(name() + "::refine()");

            // fetch output index set
            IndexSetType& index_set_out = index_set_wrapper_out.template get_index_set<0>();

            // calculate index offsets
            Index index_offsets[Shape_::dimension+1];
            EntityCounter<Shape_, 0>::offset(index_offsets, num_entities);

            // call refinement shape wrapper
            IndexRefineShapeWrapper<Shape_, cell_dim_, 0>
              ::refine(index_set_out, index_offsets, index_set_holder_in);
          }
        };

        /* ********************************************************************************************************* */

        template<
          typename Shape_,
          int cell_dim_ = Shape_::dimension>
        struct IndexRefineWrapper
        {
          // the case cell_dim_ = 1 is handled by the partial specialisation below
          static_assert(cell_dim_ > 1, "invalid cell dimension");
          static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

          typedef Shape_ ShapeType;
          typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;
          typedef IndexSetWrapper<CellType> IndexSetWrapperType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "IndexRefineWrapper<" + Shape_::name() + "," + stringify(cell_dim_) + ">";
          }

          static void refine(
            IndexSetHolderType& index_set_holder_out,
            const Index num_entities[],
            const IndexSetHolderType& index_set_holder_in)
          {
            CONTEXT(name() + "::refine()");

            // recursive call of IndexRefineWrapper
            IndexRefineWrapper<Shape_, cell_dim_ - 1>
              ::refine(index_set_holder_out, num_entities, index_set_holder_in);

            // fetch output index set wrapper
            IndexSetWrapperType& index_set_wrapper_out = index_set_holder_out
              .template get_index_set_wrapper<cell_dim_>();

            // call face wrapper
            IndexRefineFaceWrapper<Shape_, cell_dim_>
              ::refine(index_set_wrapper_out, num_entities, index_set_holder_in);
          }
        };

        template<typename Shape_>
        struct IndexRefineWrapper<Shape_, 1>
        {
          static_assert(1 <= Shape_::dimension, "invalid shape dimension");

          typedef Shape_ ShapeType;
          typedef typename Shape::FaceTraits<ShapeType, 1>::ShapeType CellType;
          typedef IndexSetWrapper<CellType> IndexSetWrapperType;
          typedef IndexSetHolder<ShapeType> IndexSetHolderType;

          static String name()
          {
            return "IndexRefineWrapper<" + Shape_::name() + ",1>";
          }

          static void refine(
            IndexSetHolderType& index_set_holder_out,
            const Index num_entities[],
            const IndexSetHolderType& index_set_holder_in)
          {
            CONTEXT(name() + "::refine()");

            // fetch output index set wrapper
            IndexSetWrapperType& index_set_wrapper_out = index_set_holder_out.template get_index_set_wrapper<1>();

            // call face wrapper
            IndexRefineFaceWrapper<Shape_, 1>::refine(index_set_wrapper_out, num_entities, index_set_holder_in);
          }
        };
      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_IDX_REF_WRAPPERS_HPP

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_MACRO_INDEX_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_MACRO_INDEX_MAPPING_HPP 1

// includes, FEAST
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        int face_dim_ = Shape_::dimension>
      struct DynamicNumFaces
      {
        static int value(int dim)
        {
          ASSERT_((dim <= face_dim_) && (dim >= 0));
          if(dim == face_dim_)
            return Shape::FaceTraits<Shape_, face_dim_>::count;
          else
            return DynamicNumFaces<Shape_, face_dim_ - 1>::value(dim);
        }
      };

      template<typename Shape_>
      struct DynamicNumFaces<Shape_, 0>
      {
        static int value(int /*dim*/)
        {
          return Shape::FaceTraits<Shape_, 0>::count;
        }
      };

      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_,
        int shape_dim_ = Shape_::dimension>
      struct MacroIndexBuilder
      {
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
        enum
        {
          num_cells = Shape::FaceTraits<Shape_, cell_dim_>::count,
          num_indices = Shape::FaceTraits<CellType, face_dim_>::count
        };

        static void build(IndexSet<num_indices>& idx)
        {
          typedef FaceIndexMapping<Shape_, cell_dim_, face_dim_> FimType;
          for(int i(0); i < num_cells; ++i)
          {
            for(int j(0); j < num_indices; ++j)
            {
              idx(i, j) = FimType::map(i, j);
            }
          }
        }
      };

      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_>
      struct MacroIndexBuilder<Shape_, cell_dim_, face_dim_, cell_dim_>
      {
        enum
        {
          num_indices = Shape::FaceTraits<Shape_, face_dim_>::count
        };

        static void build(IndexSet<num_indices>& idx)
        {
          for(int j(0); j < num_indices; ++j)
          {
            idx(0, j) = Index(j);
          }
        }
      };

      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_ = cell_dim_ - 1>
      struct MacroIndexHelper
      {
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
        static void build(IndexSetWrapper<CellType, face_dim_>& idx)
        {
          // recurse down
          MacroIndexHelper<Shape_, cell_dim_, face_dim_ - 1>::build(idx);

          // call builder
          MacroIndexBuilder<Shape_, cell_dim_, face_dim_>::build(idx.template get_index_set<face_dim_>());
        }
      };

      template<
        typename Shape_,
        int cell_dim_>
      struct MacroIndexHelper<Shape_, cell_dim_, 0>
      {
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
        static void build(IndexSetWrapper<CellType, 0>& idx)
        {
          // call builder
          MacroIndexBuilder<Shape_, cell_dim_, 0>::build(idx.template get_index_set<0>());
        }
      };

      template<
        typename Shape_,
        int cell_dim_ = Shape_::dimension>
      struct MacroIndexWrapper
      {
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
        static void build(IndexSetHolder<CellType>& idx)
        {
          // recurse down
          MacroIndexWrapper<Shape_, cell_dim_ - 1>::build(idx);

          // call helper
          MacroIndexHelper<Shape_, cell_dim_>::build(idx.template get_index_set_wrapper<cell_dim_>());
        }
      };

      template<typename Shape_>
      struct MacroIndexWrapper<Shape_, 0>
      {
        static void build(IndexSetHolder<Shape::Vertex>&)
        {
          // do nothing
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INTERN_MACRO_INDEX_MAPPING_HPP

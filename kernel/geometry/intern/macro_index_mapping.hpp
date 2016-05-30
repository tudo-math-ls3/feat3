#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_MACRO_INDEX_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_MACRO_INDEX_MAPPING_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>

namespace FEAT
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
        static constexpr int num_cells = Shape::FaceTraits<Shape_, cell_dim_>::count;
        static constexpr int num_indices = Shape::FaceTraits<CellType, face_dim_>::count;

        static void build(IndexSet<num_indices>& idx)
        {
          typedef FaceIndexMapping<Shape_, cell_dim_, face_dim_> FimType;
          for(Index i(0); i < Index(num_cells); ++i)
          {
            for(Index j(0); j < Index(num_indices); ++j)
            {
              idx(i, j) = Index(FimType::map(int(i), int(j)));
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
        static constexpr int num_indices = Shape::FaceTraits<Shape_, face_dim_>::count;

        static void build(IndexSet<num_indices>& idx)
        {
          for(Index j(0); j < Index(num_indices); ++j)
          {
            idx(0, j) = j;
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

      template<typename Shape_, int face_dim_ = Shape_::dimension-1>
      struct MacroTargetWrapper
      {
        static void build(TargetSetHolder<Shape_>& trg, const IndexSetHolder<Shape_>& ish, const Index cell_idx)
        {
          // get our target set
          TargetSet& target = trg.template get_target_set<face_dim_>();

          // get our number of local faces
          static constexpr Index num_faces = Index(Shape::FaceTraits<Shape_, face_dim_>::count);

          // validate size
          ASSERT(target.get_num_entities() == num_faces, "invalid target set size");

          // get the index set
          const IndexSet<int(num_faces)>& idx = ish.template get_index_set<Shape_::dimension, face_dim_>();

          // loop over all faces
          for(Index i(0); i < num_faces; ++i)
          {
            target[i] = idx(cell_idx, i);
          }

          // recurse down
          MacroTargetWrapper<Shape_, face_dim_-1>::build(trg, ish, cell_idx);
        }
      };

      template<typename Shape_>
      struct MacroTargetWrapper<Shape_, 0>
      {
        static void build(TargetSetHolder<Shape_>& trg, const IndexSetHolder<Shape_>& ish, const Index cell_idx)
        {
          // get our target set
          TargetSet& target = trg.template get_target_set<0>();

          // get our number of local faces
          static constexpr Index num_faces = Index(Shape::FaceTraits<Shape_,0>::count);

          // validate size
          ASSERT(target.get_num_entities() == num_faces, "invalid target set size");

          // get the index set
          const IndexSet<int(num_faces)>& idx = ish.template get_index_set<Shape_::dimension, 0>();

          // loop over all faces
          for(Index i(0); i < num_faces; ++i)
          {
            target[i] = idx(cell_idx, i);
          }
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_MACRO_INDEX_MAPPING_HPP

#pragma once
#ifndef KERNEL_GEOMETRY_CONGRUENCY_SUB_INDEX_MAPPING_HPP
#define KERNEL_GEOMETRY_CONGRUENCY_SUB_INDEX_MAPPING_HPP 1

// includes, FEAST
#include <kernel/geometry/local_index_mapping.hpp>
#include <kernel/geometry/congruency/hypercube.hpp>
#include <kernel/geometry/congruency/simplex.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Congruency
    {
      /**
       * \brief Sub-Index Mapping class template
       *
       * \todo detailed description
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_>
      class SubIndexMapping
      {
        static_assert(cell_dim_ < Shape_::dimension, "invalid cell dimension");
        static_assert(face_dim_ < cell_dim_, "invalid cell/face dimension");
        static_assert(face_dim_ >= 0, "invalid face dimension");

      protected:
        /// \cond internal
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;

        template<typename Outer_>
        class CompIndexMap
        {
        protected:
          const Outer_& _outer;
          int _idx;

        public:
          explicit CompIndexMap(
            const Outer_& outer,
            int idx)
              :
            _outer(outer),
            _idx(idx)
          {
          }

          Index operator[](int i) const
          {
            typedef LocalIndexMapping<Shape_, cell_dim_, 0> LimType;
            return _outer[LimType::map(_idx,i)];
          }
        }; // class CompIndexMap
        /// \endcond

        /// dummy enum
        enum
        {
          /// number of cells
          num_cells = Shape::FaceTraits<Shape_, cell_dim_>::count
        };

        /// cell orientation codes
        int _cell_orient[num_cells];

      public:
        template<
          typename ShapeVerts_,
          typename ShapeCells_,
          typename CellVerts_>
        explicit SubIndexMapping(
          const ShapeVerts_& shape_verts,
          const ShapeCells_& shape_cells,
          const CellVerts_& cell_verts)
        {
          typedef Sampler<CellType> SamplerType;

          // calculate orientations
          for(int cell(0); cell < num_cells; ++cell)
          {
            CompIndexMap<ShapeVerts_> comp_index_map(shape_verts, cell);
            _cell_orient[cell] = SamplerType::compare(comp_index_map, cell_verts[shape_cells[cell]]);
          }
        }

        int map(int cell, int face) const
        {
          typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;

          return IndexMapping<CellType, face_dim_>::map(_cell_orient[cell], face);
        }
      }; // class SubIndexMapping
    } // namespace Congruency
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONGRUENCY_SUB_INDEX_MAPPING_HPP

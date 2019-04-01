// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STRUCTURED_TARGET_REFINER_HPP
#define KERNEL_GEOMETRY_INTERN_STRUCTURED_TARGET_REFINER_HPP 1

// includes, FEAT
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/intern/entity_counter.hpp>
#include <kernel/geometry/intern/standard_refinement_traits.hpp>
#include <kernel/geometry/intern/struct_index_coding.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /// \todo implement for 3D
      template<int shape_dim_, int cell_dim_, int child_dim_>
      struct StructuredTargetRefiner;
      /*{
        static Index refine_simple(TargetSet&,const TargetSet&,const Index,const Index[],const Index[])
        {
          std::cout << "SPECIALISE ME: StructuredTargetRefiner<" << shape_dim_
            << "," << cell_dim_ << "," << child_dim_ << ">" << std::endl;
          return Index(0);
        }
      };*/

      // vertex -> vertex
      template<int shape_dim_>
      struct StructuredTargetRefiner<shape_dim_, 0, 0>
      {
        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSet& target_set_in,
          const Index offset,
          const Index num_slices_c[],
          const Index num_slices_f[])
        {
          const Index num_verts = target_set_in.get_num_entities();
          Index v[shape_dim_];

          for(Index i(0); i < num_verts; ++i)
          {
            // decode coarse vertex index
            StructIndexCoding<shape_dim_, 0>::decode(v, target_set_in[i], num_slices_c);

            // transform vertex index-coords
            for(int j(0); j < shape_dim_; ++j)
              v[j] *= 2;

            // encode fine vertex index
            target_set_out[offset+i] = StructIndexCoding<shape_dim_, 0>::encode(v, num_slices_f);
          }

          return num_verts;
        }
      };

      // edge -> vertex
      template<int shape_dim_>
      struct StructuredTargetRefiner<shape_dim_, 1, 0>
      {
        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSet& target_set_in,
          const Index offset,
          const Index num_slices_c[],
          const Index num_slices_f[])
        {
          const Index num_edges = target_set_in.get_num_entities();
          Index v[shape_dim_+1];

          for(Index i(0); i < num_edges; ++i)
          {
            // decode coarse edge index
            StructIndexCoding<shape_dim_, 1>::decode(v, target_set_in[i], num_slices_c);

            // transform edge index-coords
            for(int j(0); j < shape_dim_; ++j)
              v[j] *= 2;

            // funny trick :)
            ++v[v[shape_dim_]];

            // encode fine vertex index
            target_set_out[offset+i] = StructIndexCoding<shape_dim_, 0>::encode(v, num_slices_f);
          }

          return num_edges;
        }
      };

      // edge -> edge
      template<int shape_dim_>
      struct StructuredTargetRefiner<shape_dim_, 1, 1>
      {
        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSet& target_set_in,
          const Index offset,
          const Index num_slices_c[],
          const Index num_slices_f[])
        {
          const Index num_edges = target_set_in.get_num_entities();
          Index v[shape_dim_+1];

          for(Index i(0); i < num_edges; ++i)
          {
            // decode coarse edge index
            StructIndexCoding<shape_dim_, 1>::decode(v, target_set_in[i], num_slices_c);

            // transform edge index-coords
            for(int j(0); j < shape_dim_; ++j)
              v[j] *= 2;

            // add first child edge
            target_set_out[offset+2*i+0] = StructIndexCoding<shape_dim_, 1>::encode(v, num_slices_f);

            // funny trick :)
            ++v[v[shape_dim_]];

            // encode fine vertex index
            target_set_out[offset+2*i+1] = StructIndexCoding<shape_dim_, 1>::encode(v, num_slices_f);
          }

          return 2*num_edges;
        }
      };

      // 2D quad -> vertex
      template<int shape_dim_>
      struct StructuredTargetRefiner<shape_dim_, 2, 0>
      {
        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSet& target_set_in,
          const Index offset,
          const Index num_slices_c[],
          const Index num_slices_f[])
        {
          const Index num_quads = target_set_in.get_num_entities();
          Index v[shape_dim_+1];

          for(Index i(0); i < num_quads; ++i)
          {
            // decode coarse quad index
            StructIndexCoding<shape_dim_, 2>::decode(v, target_set_in[i], num_slices_c);

            // transform quad index-coords
            for(int j(0); j < shape_dim_; ++j)
              (v[j] *= 2) += 1;

            // encode fine vertex index
            target_set_out[offset+i] = StructIndexCoding<shape_dim_, 0>::encode(v, num_slices_f);
          }

          return num_quads;
        }
      };

      // 2D quad -> edge
      template</*int shape_dim_*/>
      struct StructuredTargetRefiner</*shape_dim_*/2, 2, 1>
      {
        static constexpr int shape_dim_ = 2;

        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSet& target_set_in,
          const Index offset,
          const Index num_slices_c[],
          const Index num_slices_f[])
        {
          const Index num_quads = target_set_in.get_num_entities();
          Index v[shape_dim_+1];

          for(Index i(0); i < num_quads; ++i)
          {
            // decode coarse quad index
            StructIndexCoding<shape_dim_, 2>::decode(v, target_set_in[i], num_slices_c);

            // transform quad index-coords
            for(int j(0); j < shape_dim_; ++j)
              v[j] *= 2;

            // add horizontal edges
            v[2] = 0;
            // left edge
            ++v[1];
            target_set_out[offset+4*i+0] = StructIndexCoding<shape_dim_, 1>::encode(v, num_slices_f);
            // right edge
            ++v[0];
            target_set_out[offset+4*i+1] = StructIndexCoding<shape_dim_, 1>::encode(v, num_slices_f);

            // add vertical edges
            v[2] = 1;
            // bottom edge
            --v[1];
            target_set_out[offset+4*i+2] = StructIndexCoding<shape_dim_, 1>::encode(v, num_slices_f);
            // top edge
            ++v[1];
            target_set_out[offset+4*i+3] = StructIndexCoding<shape_dim_, 1>::encode(v, num_slices_f);
          }

          return Index(4)*num_quads;
        }
      };

      // 2D quad -> quad
      template</*int shape_dim_*/>
      struct StructuredTargetRefiner</*shape_dim_*/ 2, 2, 2>
      {
        static constexpr int shape_dim_ = 2;

        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSet& target_set_in,
          const Index offset,
          const Index num_slices_c[],
          const Index num_slices_f[])
        {
          const Index num_quads = target_set_in.get_num_entities();
          Index v[shape_dim_+1];

          for(Index i(0); i < num_quads; ++i)
          {
            // decode coarse quad index
            StructIndexCoding<shape_dim_, 2>::decode(v, target_set_in[i], num_slices_c);

            // transform quad index-coords
            for(int j(0); j < shape_dim_; ++j)
              v[j] *= 2;

            // encode fine quad indices:
            // lower left child
            target_set_out[offset+4*i+0] = StructIndexCoding<shape_dim_, 2>::encode(v, num_slices_f);

            // lower right child
            ++v[0];
            target_set_out[offset+4*i+1] = StructIndexCoding<shape_dim_, 2>::encode(v, num_slices_f);

            // upper left child
            --v[0];
            ++v[1];
            target_set_out[offset+4*i+2] = StructIndexCoding<shape_dim_, 2>::encode(v, num_slices_f);

            // upper right child
            ++v[0];
            target_set_out[offset+4*i+3] = StructIndexCoding<shape_dim_, 2>::encode(v, num_slices_f);
          }

          return Index(4) * num_quads;
        }
      };

      template<int shape_dim_, int child_dim_, int cell_dim_ = shape_dim_>
      struct StructuredTargetRefineShapeWrapper
      {
        typedef Shape::Hypercube<shape_dim_> ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;

        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSetHolder<CellType>& target_set_holder_in,
          const Index num_slices_c[], const Index num_slices_f[])
        {
          // recursive call of SubSetRefineShapeWrapper
          Index offset = StructuredTargetRefineShapeWrapper<shape_dim_, child_dim_, cell_dim_-1>::refine_simple(
            target_set_out, target_set_holder_in, num_slices_c, num_slices_f);

          // call target refiner
          Index count = StructuredTargetRefiner<shape_dim_, cell_dim_, child_dim_>::refine_simple(
            target_set_out, target_set_holder_in.template get_target_set<cell_dim_>(),
            offset, num_slices_c, num_slices_f);

          // return new offset
          return offset + count;
        }
      };

      template<int shape_dim_, int cell_dim_>
      struct StructuredTargetRefineShapeWrapper<shape_dim_, cell_dim_, cell_dim_>
      {
        typedef Shape::Hypercube<shape_dim_> ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;

        static Index refine_simple(
          TargetSet& target_set_out,
          const TargetSetHolder<CellType>& target_set_holder_in,
          const Index num_slices_c[], const Index num_slices_f[])
        {
          // call target refiner
          Index count = StructuredTargetRefiner<shape_dim_, cell_dim_, cell_dim_>::refine_simple(
            target_set_out, target_set_holder_in.template get_target_set<cell_dim_>(),
            Index(0), num_slices_c, num_slices_f);

          // return offset
          return count;
        }
      };

      template<int shape_dim_, int cell_dim_ = shape_dim_>
      struct StructuredTargetRefineWrapper
      {
        typedef Shape::Hypercube<shape_dim_> ShapeType;
        typedef Shape::Hypercube<cell_dim_> CellType;

        static void refine_simple(
          TargetSetHolder<CellType>& target_set_holder_out,
          const TargetSetHolder<ShapeType>& target_set_holder_in,
          const Index num_slices_c[], const Index num_slices_f[])
        {
          // recursive call
          StructuredTargetRefineWrapper<shape_dim_, cell_dim_-1>::refine_simple(
            target_set_holder_out, target_set_holder_in, num_slices_c, num_slices_f);

          // call shape wrapper
          StructuredTargetRefineShapeWrapper<shape_dim_, cell_dim_>::refine_simple(
            target_set_holder_out.template get_target_set<cell_dim_>(),
            target_set_holder_in, num_slices_c, num_slices_f);
        }
      };

      template<int shape_dim_>
      struct StructuredTargetRefineWrapper<shape_dim_, 0>
      {
        typedef Shape::Hypercube<shape_dim_> ShapeType;
        //typedef typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType CellType;

        static void refine_simple(
          TargetSetHolder<Shape::Vertex>& target_set_holder_out,
          const TargetSetHolder<ShapeType>& target_set_holder_in,
          const Index num_slices_c[], const Index num_slices_f[])
        {
          // call shape wrapper
          StructuredTargetRefineShapeWrapper<shape_dim_, 0>::refine_simple(
            target_set_holder_out.template get_target_set<0>(),
            target_set_holder_in, num_slices_c, num_slices_f);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_STRUCTURED_TARGET_REFINER_HPP

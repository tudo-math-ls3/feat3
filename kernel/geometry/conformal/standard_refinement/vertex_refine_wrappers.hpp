#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_VTX_REF_WRAPPERS_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_VTX_REF_WRAPPERS_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal/standard_refinement/entity_counter.hpp>
#include <kernel/geometry/conformal/standard_refinement/vertex_refiner.hpp>

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
          typename VertexSet_>
        struct VertexRefineWrapper
        {
          typedef VertexSet_ VertexSetType;
          typedef IndexSetHolder<Shape_> IndexSetHolderType;
          typedef IndexSet<Shape::FaceTraits<Shape_, 0>::count> IndexSetType;

          static Index refine(
            VertexSetType& vertex_set_out,
            const VertexSetType& vertex_set_in,
            const IndexSetHolderType& index_set_holder_in)
          {
            CONTEXT("VertexRefineWrapper<...>::refine()");
            typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;

            // recursive call of VertexRefineWrapper
            Index offset = VertexRefineWrapper<FacetType, VertexSet_>
              ::refine(vertex_set_out, vertex_set_in, index_set_holder_in);

            // get index set of current shape
            const IndexSetType& index_set_in = index_set_holder_in.
              template get_index_set_wrapper<Shape_::dimension>().template get_index_set<0>();

            // call VertexRefiner
            Index num_verts = VertexRefiner<Shape_, VertexSet_>
              ::refine(offset, vertex_set_out, vertex_set_in, index_set_in);

            // validate number of created vertices
            ASSERT(num_verts == (StdRefTraits<Shape_, 0>::count * index_set_in.get_num_entities()),
              "VertexRefiner output does not match StdRefTraits prediction");

            // return new offset
            return offset + num_verts;
          }
        }; // struct VertexRefineWrapper<...>

        template<typename VertexSet_>
        struct VertexRefineWrapper<Shape::Vertex, VertexSet_>
        {
          typedef VertexSet_ VertexSetType;
          typedef IndexSetHolder<Shape::Vertex> IndexSetHolderType;

          static Index refine(
            VertexSetType& vertex_set_out,
            const VertexSetType& vertex_set_in,
            const IndexSetHolderType& /*index_set_holder_in*/)
          {
            CONTEXT("VertexRefineWrapper<Vertex,...>::refine()");

            // call VertexRefiner
            Index num_verts =  VertexRefiner<Shape::Vertex,VertexSet_>
              ::refine(0, vertex_set_out, vertex_set_in);

            // validate number of created vertices
            ASSERT(num_verts == (StdRefTraits<Shape::Vertex, 0>::count * vertex_set_in.get_num_vertices()),
              "VertexRefiner output does not match StdRefTraits prediction");

            // return new offset
            return num_verts;
          }
        }; // struct VertexRefineWrapper<Vertex,...>
      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_VTX_REF_WRAPPERS_HPP

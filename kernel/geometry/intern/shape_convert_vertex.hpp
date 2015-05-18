#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_SHAPE_CONVERT_VERTEX_HPP
#define KERNEL_GEOMETRY_INTERN_SHAPE_CONVERT_VERTEX_HPP 1

// includes, FEAST
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/intern/vertex_abacus.hpp>
#include <kernel/geometry/intern/shape_convert_traits.hpp>
#include <kernel/geometry/intern/entity_counter.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        typename VertexSet_>
      struct ShapeConvertVertex
      {
        // This generic implementation works for both Simplex<n> and Hypercube<m> shapes,
        // where 1<=n<=3 and 2<=m<=3. The case Hypercube<1> is specialised below.
        static_assert(Shape_::dimension > 0, "invalid shape");
        static_assert(Shape_::dimension <= 3, "invalid shape");

        typedef Shape_ ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index offset,
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetType& index_set_in)
        {
          CONTEXT("ShapeConvertVertex<...>::refine()");

          typedef typename IndexSetType::ConstIndexVectorReference ConstIndexVectorReference;
          typedef typename VertexSetType::VertexReference VertexReference;
          typedef typename VertexSetType::CoordType CoordType;

          // scaling factor
          static const CoordType scale = CoordType(1) / CoordType(IndexSetType::num_indices);

          // get number of cells
          Index num_cells = index_set_in.get_num_entities();

          // create a vertex abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // loop over all cells
          for(Index i(0); i < num_cells; ++i)
          {
            // get input index vector
            ConstIndexVectorReference idx_in = index_set_in[i];

            // get output vertex
            VertexReference vtx_out = vertex_set_out[offset + i];

            // clear output vertex
            abacus.clear(vtx_out);

            // add all other vertices onto it
            for(Index k(0); k < IndexSetType::num_indices; ++k)
            {
              abacus.add(vtx_out, vertex_set_in[idx_in[k]]);
            }

            // scale the output vertex
            abacus.scale(vtx_out, scale);
          }

          // return number of created vertices
          return num_cells;
        }
      };

      /**
       * \brief ShapeConvertVertex implementation for Vertex shape
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      struct ShapeConvertVertex<Shape::Vertex, VertexSet_>
      {
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index offset,
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in)
        {
          CONTEXT("ShapeConvertVertex<Vertex,...>::refine()");

          // get the number of coarse mesh vertices
          Index num_verts = vertex_set_in.get_num_vertices();
          ASSERT_(vertex_set_out.get_num_vertices() >= num_verts);

          // create a vertex-abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // loop over all vertices
          for(Index i(0); i < num_verts; ++i)
          {
            // copy source vertex
            abacus.copy(vertex_set_out[offset + i], vertex_set_in[i]);
          }

          // return number of created vertices
          return num_verts;
        }
      }; // class ShapeConvertVertex<Vertex,...>

      /**
       * \brief ShapeConvertVertex implementation for Simplex<1> shape
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      struct ShapeConvertVertex<Shape::Hypercube<1>, VertexSet_>
      {
        typedef Shape::Hypercube<1> ShapeType;
        typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
        typedef VertexSet_ VertexSetType;

        static Index refine(
          Index /*offset*/,
          VertexSetType& /*vertex_set_out*/,
          const VertexSetType& /*vertex_set_in*/,
          const IndexSetType& /*index_set_in*/)
        {
          CONTEXT("ShapeConvertVertex<Hypercube<1>,...>::refine()");

          // Hypercube edges do not generate new vertices
          return 0;
        }
      }; // class ShapeConvertVertex<Shape::Hypercube<1>,...>

      /* *************************************************************************************** */
      /* *************************************************************************************** */
      /* *************************************************************************************** */

      /**
       * \brief ShapeConvertVertexWrapper class template
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        typename VertexSet_>
      struct ShapeConvertVertexWrapper
      {
        typedef VertexSet_ VertexSetType;
        typedef IndexSetHolder<Shape_> IndexSetHolderType;
        typedef IndexSet<Shape::FaceTraits<Shape_, 0>::count> IndexSetType;

        static Index refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetHolderType& index_set_holder_in)
        {
          CONTEXT("ShapeConvertVertexWrapper<...>::refine()");
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;

          // recursive call of ShapeConvertVertexWrapper
          Index offset = ShapeConvertVertexWrapper<FacetType, VertexSet_>
            ::refine(vertex_set_out, vertex_set_in, index_set_holder_in);

          // get index set of current shape
          const IndexSetType& index_set_in = index_set_holder_in.
            template get_index_set_wrapper<Shape_::dimension>().template get_index_set<0>();

          // call VertexRefiner
          Index num_verts = ShapeConvertVertex<Shape_, VertexSet_>
            ::refine(offset, vertex_set_out, vertex_set_in, index_set_in);

          // validate number of created vertices
          ASSERT(num_verts == (ShapeConvertTraits<Shape_, 0>::count * index_set_in.get_num_entities()),
            "ShapeConvertVertex output does not match ShapeConvertTraits prediction");

          // return new offset
          return offset + num_verts;
        }
      }; // struct StandardVertexRefineWrapper<...>

      template<typename VertexSet_>
      struct ShapeConvertVertexWrapper<Shape::Vertex, VertexSet_>
      {
        typedef VertexSet_ VertexSetType;
        typedef IndexSetHolder<Shape::Vertex> IndexSetHolderType;

        static Index refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const IndexSetHolderType& /*index_set_holder_in*/)
        {
          CONTEXT("ShapeConvertVertexWrapper<Vertex,...>::refine()");

          // call VertexRefiner
          Index num_verts =  ShapeConvertVertex<Shape::Vertex,VertexSet_>
            ::refine(0, vertex_set_out, vertex_set_in);

          // validate number of created vertices
          ASSERT(num_verts == (ShapeConvertTraits<Shape::Vertex, 0>::count * vertex_set_in.get_num_vertices()),
            "ShapeConvertVertex output does not match ShapeConvertTraits prediction");

          // return new offset
          return num_verts;
        }
      }; // struct StandardVertexRefineWrapper<Vertex,...>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INTERN_SHAPE_CONVERT_VERTEX_HPP

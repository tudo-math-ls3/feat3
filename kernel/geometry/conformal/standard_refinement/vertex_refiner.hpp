#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_VERTEX_REFINER_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_VERTEX_REFINER_HPP 1

// includes, FEAST
#include <kernel/geometry/std_ref_traits.hpp>
#include <kernel/geometry/vertex_set.hpp>
#include <kernel/geometry/vertex_abacus.hpp>
#include <kernel/geometry/conformal/index_set.hpp>

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
        struct VertexRefiner;

        template<typename VertexSet_>
        struct VertexRefiner<Shape::Vertex, VertexSet_>
        {
          typedef VertexSet_ VertexSetType;

          static Index refine(
            Index offset,
            VertexSetType& vertex_set_out,
            const VertexSetType& vertex_set_in)
          {
            CONTEXT("VertexRefiner<Vertex,...>::refine()");

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
        }; // class VertexRefiner<Vertex,...>

        template<
          int cell_dim_,
          typename VertexSet_>
        struct VertexRefiner<Shape::Hypercube<cell_dim_>, VertexSet_>
        {
          typedef Shape::Hypercube<cell_dim_> ShapeType;
          typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
          typedef VertexSet_ VertexSetType;

          static Index refine(
            Index offset,
            VertexSetType& vertex_set_out,
            const VertexSetType& vertex_set_in,
            const IndexSetType& index_set_in)
          {
            CONTEXT("VertexRefiner<...,Hypercube<...>>::refine()");

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
              for(int k(0); k < IndexSetType::num_indices; ++k)
              {
                abacus.add(vtx_out, vertex_set_in[idx_in[k]]);
              }

              // scale the output vertex
              abacus.scale(vtx_out, scale);
            }

            // return number of created vertices
            return num_cells;
          }
        }; // struct VertexRefiner<Hypercube<...>,...>

        //////////////////////////////////////////////////

        // specialisation for "simplex edges"
        template<
          typename VertexSet_>
        struct VertexRefiner<Shape::Simplex<1>, VertexSet_>
        {
          typedef Shape::Simplex<1> ShapeType;
          typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
          typedef VertexSet_ VertexSetType;

          static Index refine(
            Index offset,
            VertexSetType& vertex_set_out,
            const VertexSetType& vertex_set_in,
            const IndexSetType& index_set_in)
          {
            CONTEXT("VertexRefiner<...,Simplex<...>>::refine()");

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
              for(int k(0); k < IndexSetType::num_indices; ++k)
              {
                abacus.add(vtx_out, vertex_set_in[idx_in[k]]);
              }

              // scale the output vertex
              abacus.scale(vtx_out, scale);
            }

            // return number of created vertices
            return num_cells;
          }
        }; // struct VertexRefiner<Simplex<...>,...>

        // specialisation for triangles
        template<
          typename VertexSet_>
        struct VertexRefiner<Shape::Simplex<2>, VertexSet_>
        {
          typedef Shape::Simplex<2> ShapeType;
          typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
          typedef VertexSet_ VertexSetType;

          static Index refine(
            Index /*offset*/,
            VertexSetType& /*vertex_set_out*/,
            const VertexSetType& /*vertex_set_in*/,
            const IndexSetType& /*index_set_in*/)
          {
            CONTEXT("VertexRefiner<...,Simplex<...>>::refine()");

            // return number of created vertices
            return 0;
          }
        }; // struct VertexRefiner<Simplex<...>,...>

        // specialisation for tetraeder

        template<
          typename VertexSet_>
        struct VertexRefiner<Shape::Simplex<3>, VertexSet_>
        {
          typedef Shape::Simplex<3> ShapeType;
          typedef IndexSet<Shape::FaceTraits<ShapeType, 0>::count> IndexSetType;
          typedef VertexSet_ VertexSetType;

          static Index refine(
            Index offset,
            VertexSetType& vertex_set_out,
            const VertexSetType& vertex_set_in,
            const IndexSetType& index_set_in)
          {
            CONTEXT("VertexRefiner<...,Simplex<...>>::refine()");

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
              for(int k(0); k < IndexSetType::num_indices; ++k)
              {
                abacus.add(vtx_out, vertex_set_in[idx_in[k]]);
              }

              // scale the output vertex
              abacus.scale(vtx_out, scale);
            }

            // return number of created vertices
            return num_cells;
          }
        }; // struct VertexRefiner<Simplex<...>,...>

        ////////////////////////////////////////////////////// by CC

      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_VERTEX_REFINER_HPP

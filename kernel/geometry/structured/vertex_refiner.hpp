#pragma once
#ifndef KERNEL_GEOMETRY_STRUCTURED_VERTEX_REFINER_HPP
#define KERNEL_GEOMETRY_STRUCTURED_VERTEX_REFINER_HPP 1

// includes, FEAST
#include <kernel/geometry/vertex_set.hpp>
#include <kernel/geometry/vertex_abacus.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Structured namespace
     *
     * This namespace encapsulates classes relates to the implementation of the structured mesh.
     */
    namespace Structured
    {
      /// \cond internal
      template<
        typename Shape_,
        typename VertexSet_>
      class VertexRefiner;

      /**
       * \brief Structured VertexRefiner implementation for Hypercube<1>
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      class VertexRefiner<Shape::Hypercube<1>, VertexSet_>
      {
      public:
        typedef VertexSet_ VertexSetType;

        static void refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const Index num_slices_coarse[])
        {
          CONTEXT("VertexRefiner<Hypercube<1>,...>::refine");

          // create a vertex-abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // get number of X-slices
          const Index n = num_slices_coarse[0];

          // copy first vertex
          abacus.copy(vertex_set_out[0], vertex_set_in[0]);

          // loop over all slices
          for(Index i(0); i < n; ++i)
          {
            // calculate edge midpoint
            abacus.mean(vertex_set_out[2*i + 1], vertex_set_in[i], vertex_set_in[i + 1]);

            // copy right vertex
            abacus.copy(vertex_set_out[2*i + 2], vertex_set_in[i + 1]);
          }
        }
      }; // class VertexRefiner<Hypercube<1>,...>


      /**
       * \brief Structured VertexRefiner implementation for Hypercube<2>
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      class VertexRefiner<Shape::Hypercube<2>, VertexSet_>
      {
      public:
        typedef VertexSet_ VertexSetType;

        static void refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const Index num_slices_coarse[])
        {
          CONTEXT("VertexRefiner<Hypercube<2>,...>::refine");

          // create a vertex-abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // get number of X- and Y-slices
          const Index m = num_slices_coarse[0];
          const Index n = num_slices_coarse[1];

          // copy lower left vertex
          abacus.copy(vertex_set_out[0], vertex_set_in[0]);

          // loop over all edges in bottom row
          for(Index i(0); i < m; ++i)
          {
            // calculate edge midpoint
            abacus.mean(vertex_set_out[2*i + 1], vertex_set_in[i], vertex_set_in[i + 1]);

            // copy right vertex
            abacus.copy(vertex_set_out[2*i + 2], vertex_set_in[i + 1]);
          }

          // loop over all rows
          for(Index j(0); j < n; ++j)
          {
            // calculate offsets
            Index k0 =      j *(m + 1);     // lower-left coarse mesh vertex
            Index k1 = (1 + j)*(m + 1);     // upper-left coarse mesh vertex
            Index l0 = (2*j + 1)*(2*m + 1); // left fine mesh vertex; generated from coarse mesh edge midpoint
            Index l1 = 2*(j + 1)*(2*m + 1); // upper-left fine mesh vertex

            // calculate left edge midpoint
            abacus.mean(vertex_set_out[l0], vertex_set_in[k0], vertex_set_in[k1]);

            // copy upper left vertex
            abacus.copy(vertex_set_out[l1], vertex_set_in[k1]);

            // loop over all cells in current row
            for(Index i(0); i < m; ++i)
            {
              // calculate quad midpoint
              abacus.mean(vertex_set_out[l0 + 2*i + 1],
                vertex_set_in[k0 + i], vertex_set_in[k0 + i + 1],
                vertex_set_in[k1 + i], vertex_set_in[k1 + i + 1]);

              // calculate right edge midpoint
              abacus.mean(vertex_set_out[l0 + 2*i + 2], vertex_set_in[k0 + i + 1], vertex_set_in[k1 + i + 1]);

              // calculate upper edge midpoint
              abacus.mean(vertex_set_out[l1 + 2*i + 1], vertex_set_in[k1 + i], vertex_set_in[k1 + i + 1]);

              // copy upper right vertex
              abacus.copy(vertex_set_out[l1 + 2*i + 2], vertex_set_in[k1 + i + 1]);
            }
          }
        }
      }; // class VertexRefiner<Hypercube<2>,...>
      /// \endcond
    } // namespace Structured
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_STRUCTURED_VERTEX_REFINER_HPP

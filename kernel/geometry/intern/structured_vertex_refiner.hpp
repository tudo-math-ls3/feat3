#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STRUCTURED_VERTEX_REFINER_HPP
#define KERNEL_GEOMETRY_INTERN_STRUCTURED_VERTEX_REFINER_HPP 1

// includes, FEAT
#include <kernel/geometry/intern/vertex_abacus.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Shape_,
        typename VertexSet_>
      class StructuredVertexRefiner;

      /**
       * \brief Structured StructuredVertexRefiner implementation for Hypercube<1>
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      class StructuredVertexRefiner<Shape::Hypercube<1>, VertexSet_>
      {
      public:
        typedef VertexSet_ VertexSetType;

        static void refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const Index num_slices_coarse[])
        {
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
      }; // class StructuredVertexRefiner<Hypercube<1>,...>


      /**
       * \brief Structured StructuredVertexRefiner implementation for Hypercube<2>
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      class StructuredVertexRefiner<Shape::Hypercube<2>, VertexSet_>
      {
      public:
        typedef VertexSet_ VertexSetType;

        static void refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const Index num_slices_coarse[])
        {
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
      }; // class StructuredVertexRefiner<Hypercube<2>,...>

       /**
       * \brief Structured StructuredVertexRefiner implementation for Hypercube<3>
       */
      template<typename VertexSet_>
      class StructuredVertexRefiner<Shape::Hypercube<3>, VertexSet_>
      {
      public:
        typedef VertexSet_ VertexSetType;

        static void refine(
          VertexSetType& vertex_set_out,
          const VertexSetType& vertex_set_in,
          const Index num_slices_coarse[])
        {
          // create a vertex-abacus object
          VertexAbacus<VertexSetType> abacus(vertex_set_in);

          // get number of X-, Y- and Z-slices
          const Index m = num_slices_coarse[0];
          const Index n = num_slices_coarse[1];
          const Index l = num_slices_coarse[2];

          // copy lower left front vertex
          abacus.copy(vertex_set_out[0], vertex_set_in[0]);

          // indices of the bottom layer (low,up -> Y-axis, left,right -> X-axis)

          // indices of the lower edge
          for(Index i(0); i < m; ++i)
          {
            // calculate edge midpoint
            abacus.mean(vertex_set_out[2*i + 1], vertex_set_in[i], vertex_set_in[i + 1]);

            // copy right vertex
            abacus.copy(vertex_set_out[2*i + 2], vertex_set_in[i + 1]);
          }

          // loop over all Y-slices
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

            // loop over all cells in the current Y-slice
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
          // indices of the bottom face done

          // loop over all Z-slices of the coarse mesh
          for(Index k(0); k < l; ++k)
          {
            // coarse and fine mesh index of the left front vertex of the current slice (= offsets)
            Index oc = (k+1)*(m+1)*(n+1);
            Index of = 2*(k+1)*(2*m+1)*(2*n+1);

            // copy lower left front vertex
            abacus.copy(vertex_set_out[of], vertex_set_in[oc]);

            // indices of the front edge
            for(Index i(0); i < m; ++i)
            {
              // calculate edge midpoint
              abacus.mean(vertex_set_out[of + 2*i + 1], vertex_set_in[oc + i], vertex_set_in[oc + i + 1]);

              // copy right vertex
              abacus.copy(vertex_set_out[of + 2*i + 2], vertex_set_in[oc + i + 1]);
            }

            // loop over all Y-slices
            for(Index j(0); j < n; ++j)
            {
              // calculate (local) offsets
              Index k0 =      j *(m + 1);     // lower-left coarse mesh vertex
              Index k1 = (1 + j)*(m + 1);     // upper-left coarse mesh vertex
              Index l0 = (2*j + 1)*(2*m + 1); // left fine mesh vertex; generated from coarse mesh edge midpoint
              Index l1 = 2*(j + 1)*(2*m + 1); // upper-left fine mesh vertex

              // calculate left edge midpoint
              abacus.mean(vertex_set_out[l0 + of], vertex_set_in[k0 + oc], vertex_set_in[k1 + oc]);

              // copy upper left vertex
              abacus.copy(vertex_set_out[l1 + of], vertex_set_in[k1+ oc]);

              // loop over all cells in the current Y-slice
              for(Index i(0); i < m; ++i)
              {
                // calculate quad midpoint
                abacus.mean(vertex_set_out[l0 + 2*i + 1 + of],
                vertex_set_in[k0 + i + oc], vertex_set_in[k0 + i + 1 + oc],
                vertex_set_in[k1 + i + oc], vertex_set_in[k1 + i + 1 + oc]);

                // calculate right edge midpoint
                abacus.mean(vertex_set_out[l0 + 2*i + 2 + of], vertex_set_in[k0 + i + 1 + oc], vertex_set_in[k1 + i + 1 + oc]);

               // calculate upper edge midpoint
                abacus.mean(vertex_set_out[l1 + 2*i + 1 + of], vertex_set_in[k1 + i + oc], vertex_set_in[k1 + i + 1 + oc]);

                // copy upper right vertex
                abacus.copy(vertex_set_out[l1 + 2*i + 2 + of], vertex_set_in[k1 + i + 1 + oc]);
              } // i-loop (cells)
            } // j-loop (Y-slices)

            // calculate indices of the layer between the original coarse mesh layers

            // loop over all Y-slices
            for(Index j(0); j < 2*n+1; ++j)
            {
              // loop over all cells
              for(Index i(0); i < 2*m+1; ++i)
              {
                // calculate midpoint
                abacus.mean(vertex_set_out[ i + j*(2*m+1) + of - (2*m+1)*(2*n+1)],
                vertex_set_out[i + j*(2*m+1) + of], vertex_set_out[i + j*(2*m+1) + of - 2*(2*m+1)*(2*n+1)]);
              } // i-loop
            } // j-loop
          } //k-loop
        }
      }; // class StructuredVertexRefiner<Hypercube<3>,...>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_STRUCTURED_VERTEX_REFINER_HPP

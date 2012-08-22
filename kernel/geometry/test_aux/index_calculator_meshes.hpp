#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_INDEX_CALCULATOR_MESHES_HPP
#define KERNEL_GEOMETRY_TEST_AUX_INDEX_CALCULATOR_MESHES_HPP 1

// includes, FEAST
#include <kernel/geometry/index_calculator.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {

      typedef IndexTree< Shape::Hypercube<2> > QuadIndexTree;
      typedef IndexTree< Shape::Hypercube<1> > HEdgeIndexTree;
      typedef IndexTree< Shape::Simplex<2> > TriaIndexTree;
      typedef IndexTree< Shape::Simplex<1> > SEdgeIndexTree;

      /**
       * \brief Checks, if the vertex-at-edge index set of the tetris-hexa mesh was
       * parsed correctly.
       *
       * \author Constantin Christof
       */
      void validate_hypercube_edge_index_set(const HEdgeIndexTree& tree);

      /**
       * \brief Checks, if the vertex-at-quad index set of the tetris-hexa mesh was
       * parsed correctly.
       *
       * \author Constantin Christof
       */
      void validate_hypercube_quad_index_set(const QuadIndexTree& tree);

      /**
       * \brief Checks, if the vertex-at-edge index set of the big-tetra mesh was
       * parsed correctly.
       *
       * \author Constantin Christof
       */
      void validate_simplex_edge_index_set(const SEdgeIndexTree& tree);

      /**
       * \brief Checks, if the vertex-at-triangle index set of the big-tetra mesh was
       * parsed correctly.
       *
       * \author Constantin Christof
       */

      void validate_simplex_triangle_index_set(const TriaIndexTree& tree);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_INDEX_CALCULATOR_MESHES_HPP

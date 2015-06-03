#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_TETRIS_QUAD_HPP
#define KERNEL_GEOMETRY_TEST_AUX_TETRIS_QUAD_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      typedef ConformalMesh<Shape::Quadrilateral> QuadMesh;
      typedef MeshPart<QuadMesh> QuadSubMesh;
      typedef MeshPart<QuadMesh> QuadCellSubSet;

      /**
       * \brief Creates the 2D tetris mesh.
       *
       * \verbatim
         7----L----8----M----9
         |         |         |
         I   Q_2   J   Q_3   K
         |         |         |
         3----F----4----G----5----H----6
                   |         |         |
                   C   Q_0   D   Q_1   E
                   |         |         |
                   0----A----1----B----2 \endverbatim
       * \note
       * - All horizontal edges are pointing from left to right.
       * - All vertical edges are pointing from down to up.
       * - For all quads
       *   - the lower left vertex is the first local vertex
       *   - the lower right vertex is the second local vertex
       *
       * \author Peter Zajac
       */
      QuadMesh* create_tetris_mesh_2d();

      /**
       * \brief Validates the refined 2d tetris mesh.
       *
       * \verbatim
           7----W---21----X----8----Y---22----Z----9
           |         |         |         |         |
           R   Q_10  J'  Q_11  T   Q_14  N'  Q_15  V
           |         |         |         |         |
          18---K'---25---L'---19---O'---26---P'---20
           |         |         |         |         |
           Q   Q_8   I'  Q_9   S   Q_12  M'  Q_13  U
           |         |         |         |         |
           3----K---15----L----4----M---16----N----5----O---17----P----6
                               |         |         |         |         |
                               F   Q_2   B'  Q_3   H   Q_6   F'  Q_7   J
                               |         |         |         |         |
                              12---C'---23---D'---13---G'---24---H'---14
                               |         |         |         |         |
                               E   Q_0   A'  Q_1   G   Q_4   E'  Q_5   I
                               |         |         |         |         |
                               0----A---10----B----1----C---11----D----2 \endverbatim
       *
       * \author Peter Zajac
       */
      void validate_refined_tetris_mesh_2d(const QuadMesh& mesh);

      QuadSubMesh* create_tetris_edge_submesh_2d();

      void validate_refined_tetris_edge_submesh_2d(const QuadSubMesh& mesh);

      QuadSubMesh* create_tetris_quad_submesh_2d();

      void validate_refined_tetris_quad_submesh_2d(const QuadSubMesh& mesh);

      QuadSubMesh* create_tetris_quad_edge_submesh_2d();

      void validate_refined_tetris_quad_edge_submesh_2d(const QuadSubMesh& mesh);

      QuadCellSubSet* create_tetris_quad_cellsubset_2d();

      void validate_refined_tetris_quad_cellsubset_2d(const QuadCellSubSet& subset);

      QuadCellSubSet* create_tetris_quad_edge_cellsubset_2d();

      void validate_refined_tetris_quad_edge_cellsubset_2d(const QuadCellSubSet& subset);

      void validate_tetris_quad_boundary_cellsubset_2d(const QuadCellSubSet& subset);
      void validate_refined_tetris_quad_boundary_cellsubset_2d(const QuadCellSubSet& subset);
    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_TETRIS_QUAD_HPP

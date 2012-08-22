#include <kernel/geometry/test_aux/standard_quad.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace TestAux
    {
      QuadMesh* create_quad_mesh_2d(int orientation)
      {
        Index num_entities[] =
        {
          4, // vertices
          4, // edges
          1   // quads
        };

        // create mesh
        QuadMesh* mesh = new QuadMesh(num_entities);

        //// first possibility ////

        //             e_1
        //    v_2_______>_______v_3
        //     |                 |
        //     |     -->---      |
        //     |      \          |
        //     |       ^         |
        // e_2 ^        \        ^ e_3
        //     |         \       |
        //     |     -->---      |
        //     |________>________|
        //    v_0      e_0      v_1
        //

        // set up vertex coordinates array
        static const Real vtx0[4*2] =
        {
          0.0, 0.0,
          1.0, 0.0,
          0.0, 1.0,
          1.0, 1.0
        };

        // set up vertices-at-edge array
        static const Index v_e0[4*2] =
        {
          0, 1,
          2, 3,
          0, 2,
          1, 3
        };

        // set up vertices-at-quad array
        static const Index v_q0[1*4] =
        {
          0, 1, 2, 3
        };

        // set up edges-at-quad array
        static const Index e_q0[1*4] =
        {
          0,  1,  2,  3
        };

        ////second possibility ////

        //             e_0
        //    v_1_______<_______v_2
        //     |                 |
        //     |     |     /|    |
        //     |     |    / |    |
        //     |     |   /  |    |
        // e_2 v     ^  v   ^    v e_1
        //     |     | /    |    |
        //     |     |/     |    |
        //     |________<________|
        //    v_3      e_3      v_0
        //

        // set up vertex coordinates array
        static const Real vtx1[4*2] =
        {
          1.0, 0.0,
          0.0, 1.0,
          1.0, 1.0,
          0.0, 0.0
        };

        // set up vertices-at-edge array
        static const Index v_e1[4*2] =
        {
          2, 1,
          2, 0,
          1, 3,
          0, 3
        };

        // set up vertices-at-quad array
        static const Index v_q1[1*4] =
        {
          0, 2, 3, 1
        };

        // set up edges-at-quad array
        static const Index e_q1[1*4] =
        {
          1,  2,  3,  0
        };

        ////third possibility ////

        //             e_2
        //    v_2_______>_______v_1
        //     |                 |
        //     |     --<---      |
        //     |      \          |
        //     |       v         |
        // e_3 ^        \        ^ e_0
        //     |         \       |
        //     |     --<---      |
        //     |________>________|
        //    v_0      e_1      v_3
        //

        // set up vertex coordinates array
        static const Real vtx2[4*2] =
        {
          0.0, 0.0,
          1.0, 1.0,
          0.0, 1.0,
          1.0, 0.0
        };

        // set up vertices-at-edge array
        static const Index v_e2[4*2] =
        {
          3, 1,
          0, 3,
          2, 1,
          0, 2
        };

        // set up vertices-at-quad array
        static const Index v_q2[1*4] =
        {
          1, 2, 3, 0
        };

        // set up edges-at-quad array
        static const Index e_q2[1*4] =
        {
          2,  1,  0,  3
        };

        ////fourth possibility ////

        //             e_2
        //    v_1_______>_______v_2
        //     |                 |
        //     |     |     /|    |
        //     |     |    / |    |
        //     |     |   /  |    |
        // e_1 ^     v  ^   v    v e_0
        //     |     | /    |    |
        //     |     |/     |    |
        //     |________<________|
        //    v_0      e_3      v_3
        //

        // set up vertex coordinates array
        static const Real vtx3[4*2] =
        {
          0.0, 0.0,
          0.0, 1.0,
          1.0, 1.0,
          1.0, 0.0
        };

        // set up vertices-at-edge array
        static const Index v_e3[4*2] =
        {
          2, 3,
          0, 1,
          1, 2,
          3, 0
        };

        // set up vertices-at-quad array
        static const Index v_q3[1*4] =
        {
          1, 0, 2, 3
        };

        // set up edges-at-quad array
        static const Index e_q3[1*4] =
        {
          1,  0,  2,  3
        };

        switch(orientation)
        {
          case 0:
            copy_vtx(mesh->get_vertex_set(), vtx0);
            copy_idx(mesh->get_index_set<1,0>(), v_e0);
            copy_idx(mesh->get_index_set<2,0>(), v_q0);
            copy_idx(mesh->get_index_set<2,1>(), e_q0);
            break;
          case 1:
            copy_vtx(mesh->get_vertex_set(), vtx1);
            copy_idx(mesh->get_index_set<1,0>(), v_e1);
            copy_idx(mesh->get_index_set<2,0>(), v_q1);
            copy_idx(mesh->get_index_set<2,1>(), e_q1);
            break;
          case 2:
            copy_vtx(mesh->get_vertex_set(), vtx2);
            copy_idx(mesh->get_index_set<1,0>(), v_e2);
            copy_idx(mesh->get_index_set<2,0>(), v_q2);
            copy_idx(mesh->get_index_set<2,1>(), e_q2);
            break;
          case 3:
            copy_vtx(mesh->get_vertex_set(), vtx3);
            copy_idx(mesh->get_index_set<1,0>(), v_e3);
            copy_idx(mesh->get_index_set<2,0>(), v_q3);
            copy_idx(mesh->get_index_set<2,1>(), e_q3);
            break;
        }
        // okay
        return mesh;
      } // create_quad_mesh_2d

      void validate_refined_quad_mesh_2d(const QuadMesh& mesh, int orientation)
      {
        // validate sizes
        if(mesh.get_num_entities(0) != 9)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 12)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 4)
          throw String("Quad count mismatch");

        //// first possibility ////

        // vertex coordinates array
        static const Real vtx0[] =
        {
          0.0, 0.0, // coarse mesh vertices (0,...,3)
          1.0, 0.0,
          0.0, 1.0,
          1.0, 1.0,
          0.5, 0.0,
          0.5, 1.0,
          0.0, 0.5,
          1.0, 0.5,
          0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e0[] =
        {
          0, 4,
          4, 1,
          2, 5,
          5, 3,
          0, 6,
          6, 2,
          1, 7,
          7, 3,
          4, 8,
          8, 5,
          6, 8,
          8, 7
        };

        // vertices-at-quad
        static const Index v_q0[] =
        {
          0, 4, 6, 8,
          4, 1, 8, 7,
          6, 8, 2, 5,
          8, 7, 5, 3
        };

        // edges-at-quad
        static const Index e_q0[] =
        {
          0, 10, 4, 8,
          1, 11, 8, 6,
          10, 2, 5, 9,
          11, 3, 9, 7
        };

        //// second possibility ////

        // vertex coordinates array
        static const Real vtx1[] =
        {
          1.0, 0.0,
          0.0, 1.0,
          1.0, 1.0,
          0.0, 0.0,
          0.5, 1.0,
          1.0, 0.5,
          0.0, 0.5,
          0.5, 0.0,
          0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e1[] =
        {
          2, 4,
          4, 1,
          2, 5,
          5, 0,
          1, 6,
          6, 3,
          0, 7,
          7, 3,
          5, 8,
          8, 6,
          7, 8,
          8, 4
        };

        // vertices-at-quad
        static const Index v_q1[] =
        {
          0, 5, 7, 8,
          5, 2, 8, 4,
          7, 8, 3, 6,
          8, 4, 6, 1
        };

        // edges-at-quad
        static const Index e_q1[] =
        {
          3, 10, 6, 8,
          2, 11, 8, 0,
          10, 5, 7, 9,
          11, 4, 9, 1
        };

        //// third possibility ////

        // vertex coordinates array
        static const Real vtx2[] =
        {
          0.0, 0.0,
          1.0, 1.0,
          0.0, 1.0,
          1.0, 0.0,
          1.0, 0.5,
          0.5, 0.0,
          0.5, 1.0,
          0.0, 0.5,
          0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e2[] =
        {
          3, 4,
          4, 1,
          0, 5,
          5, 3,
          2, 6,
          6, 1,
          0, 7,
          7, 2,
          6, 8,
          8, 5,
          4, 8,
          8, 7
        };

        // vertices-at-quad
        static const Index v_q2[] =
        {
          1, 6, 4, 8,
          6, 2, 8, 7,
          4, 8, 3, 5,
          8, 7, 5, 0
        };

        // edges-at-quad
        static const Index e_q2[] =
        {
          5, 10, 1, 8,
          4, 11, 8, 7,
          10, 3, 0, 9,
          11, 2, 9, 6
        };

        //// fourth possibility ////

        // vertex coordinates array
        static const Real vtx3[] =
        {
          0.0, 0.0,
          0.0, 1.0,
          1.0, 1.0,
          1.0, 0.0,
          1.0, 0.5,
          0.0, 0.5,
          0.5, 1.0,
          0.5, 0.0,
          0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e3[] =
        {
          2, 4,
          4, 3,
          0, 5,
          5, 1,
          1, 6,
          6, 2,
          3, 7,
          7, 0,
          5, 8,
          8, 4,
          6, 8,
          8, 7
        };

        // vertices-at-quad
        static const Index v_q3[] =
        {
          1, 5, 6, 8,
          5, 0, 8, 7,
          6, 8, 2, 4,
          8, 7, 4, 3
        };

        // edges-at-quad
        static const Index e_q3[] =
        {
          3, 10, 4, 8,
          2, 11, 8, 7,
          10, 0, 5, 9,
          11, 1, 9, 6
        };


        switch(orientation)
        {
          case 0:
            // check vertex coordinates array
            if(!comp_vtx(mesh.get_vertex_set(), vtx0))
              throw String("Vertex coordinate refinement failure");

            // check vertices-at-edge array
            if(!comp_idx(mesh.get_index_set<1,0>(), v_e0))
              throw String("Vertex-At-Edge index set refinement failure");

            // check vertices-at-quad
            if(!comp_idx(mesh.get_index_set<2,0>(), v_q0))
              throw String("Vertex-At-Quad index set refinement failure");

            // check edges-at-quad
            if(!comp_idx(mesh.get_index_set<2,1>(), e_q0))
              throw String("Edge-At-Quad index set refinement failure");
            break;

          case 1:
            // check vertex coordinates array
            if(!comp_vtx(mesh.get_vertex_set(), vtx1))
              throw String("Vertex coordinate refinement failure");

            // check vertices-at-edge array
            if(!comp_idx(mesh.get_index_set<1,0>(), v_e1))
              throw String("Vertex-At-Edge index set refinement failure");

            // check vertices-at-quad
            if(!comp_idx(mesh.get_index_set<2,0>(), v_q1))
              throw String("Vertex-At-Quad index set refinement failure");

            // check edges-at-quad
            if(!comp_idx(mesh.get_index_set<2,1>(), e_q1))
              throw String("Edge-At-Quad index set refinement failure");
            break;

          case 2:
            // check vertex coordinates array
            if(!comp_vtx(mesh.get_vertex_set(), vtx2))
              throw String("Vertex coordinate refinement failure");

            // check vertices-at-edge array
            if(!comp_idx(mesh.get_index_set<1,0>(), v_e2))
              throw String("Vertex-At-Edge index set refinement failure");

            // check vertices-at-quad
            if(!comp_idx(mesh.get_index_set<2,0>(), v_q2))
              throw String("Vertex-At-Quad index set refinement failure");

            // check edges-at-quad
            if(!comp_idx(mesh.get_index_set<2,1>(), e_q2))
              throw String("Edge-At-Quad index set refinement failure");
            break;

          case 3:
            // check vertex coordinates array
            if(!comp_vtx(mesh.get_vertex_set(), vtx3))
              throw String("Vertex coordinate refinement failure");

            // check vertices-at-edge array
            if(!comp_idx(mesh.get_index_set<1,0>(), v_e3))
              throw String("Vertex-At-Edge index set refinement failure");

            // check vertices-at-quad
            if(!comp_idx(mesh.get_index_set<2,0>(), v_q3))
              throw String("Vertex-At-Quad index set refinement failure");

            // check edges-at-quad
            if(!comp_idx(mesh.get_index_set<2,1>(), e_q3))
              throw String("Edge-At-Quad index set refinement failure");
            break;
        } //switch
      } // validate_refined_quad_mesh_2d

    } // namespace TestAux
  } // namespace Geometry
} // namespace FEAST

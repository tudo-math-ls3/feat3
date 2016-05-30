#include <kernel/geometry/test_aux/standard_hexa.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

namespace FEAT
{
  namespace Geometry
  {
    namespace TestAux
    {
      HexaMesh* create_hexa_mesh_3d(int orientation)
      {
        Index num_entities[] =
        {
          8, // vertices
          12, // edges
          6,   // quads
          1   // cubes
        };

        // create mesh
        HexaMesh* mesh = new HexaMesh(num_entities);

        //// first possibility ////

        //         v_6__________e_3___>_____v_7
        //          /|                     /|
        //      e_6^ |        q_1      e_7^ |
        //        /__|________>__________/  |
        //    v_4|   |       e_2      v_5|  |
        //       |   |                   |  |
        //       |   |e_10      q_3      |q |e_11
        //       | q ^         x         |5 ^
        //    e_8| 4 |       q_2      e_9|  |
        //       ^   |                   ^  |
        //       |   |v_2_______e_1__>___|__|v_3
        //       |  /                    |  /
        //       | ^e_4       q_0        | ^e_5
        //       |/___________>__________|/
        //      v_0          e_0        v_1
        //
        // orientation: v_0-v_1-v_2-v_3-v_4-v_5-v_6-v_7
        // quad orientation:
        // q_0 : v_0-v_1-v_2-v_3
        // q_1 : v_4-v_5-v_6-v_7
        // q_2 : v_0-v_1-v_4-v_5
        // q_3 : v_2-v_3-v_6-v_7
        // q_4 : v_0-v_2-v_4-v_6
        // q_5 : v_1-v_3-v_5-v_7

        // set up vertex coordinates array
        static const Real vtx0[3*8] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          0.0, 1.0, 1.0,
          1.0, 1.0, 1.0
        };

        // set up vertices-at-edge array
        static const Index v_e0[12*2] =
        {
          0, 1,
          2, 3,
          4, 5,
          6, 7,
          0, 2,
          1, 3,
          4, 6,
          5, 7,
          0, 4,
          1, 5,
          2, 6,
          3, 7
        };

        // set up vertices-at-quad array
        static const Index v_q0[6*4] =
        {
          0, 1, 2, 3,
          4, 5, 6, 7,
          0, 1, 4, 5,
          2, 3, 6, 7,
          0, 2, 4, 6,
          1, 3, 5, 7
        };

        // set up vertices-at-cube array
        static const Index v_c0[1*8] =
        {
          0, 1, 2, 3, 4, 5, 6, 7
        };

        // set up edges-at-quad array
        static const Index e_q0[6*4] =
        {
          0, 1, 4, 5,
          2, 3, 6, 7,
          0, 2, 8, 9,
          1, 3, 10, 11,
          4, 6, 8, 10,
          5, 7, 9, 11
        };

        // set up edges-at-cube array
        static const Index e_c0[1*12] =
        {
          0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
        };

        // set up quad-at-cube array
        static const Index q_c0[1*6] =
        {
          0, 1, 2, 3, 4, 5
        };

        //// second possibility ////

        //         v_5__________e_7__>______v_2
        //          /|                     /|
        //      e_6^ |        q_2      e_4v |
        //        /__|________<__________/  |
        //    v_4|   |       e_5      v_1|  |
        //       |   |                   |  |
        //       |   |e_10      q_5      |q |e_9
        //       | q v         x         |1 |
        //    e_8| 3 |       q_4     e_11|  |
        //       v   |                   v  |
        //       |   |v_3_______e_2__<___|__|v_7
        //       |  /                    |  /
        //       | ve_3       q_0        | ^e_1
        //       |/___________>__________|/
        //      v_0          e_0        v_6
        //
        // cube orientation: v_1-v_2-v_6-v_7-v_4-v_5-v_0-v_3
        // quad orientation:
        // q_0 : v_7-v_3-v_6-v_0
        // q_1 : v_7-v_6-v_2-v_1
        // q_2 : v_2-v_5-v_1-v_4
        // q_3 : v_0-v_3-v_4-v_5
        // q_4 : v_6-v_1-v_0-v_4
        // q_5 : v_5-v_3-v_2-v_7

        // set up vertex coordinates array
        static const Real vtx1[3*8] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 1.0,
          1.0, 1.0, 1.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          0.0, 1.0, 1.0,
          1.0, 0.0, 0.0,
          1.0, 1.0, 0.0,
        };

        // set up vertices-at-edge array
        static const Index v_e1[12*2] =
        {
          0, 6,
          6, 7,
          7, 3,
          3, 0,

          2, 1,
          1, 4,
          4, 5,
          5, 2,

          4, 0,
          2, 7,
          5, 3,
          1, 6
        };

        // set up vertices-at-quad array
        static const Index v_q1[6*4] =
        {
          7, 3, 6, 0,
          7, 6, 2, 1,
          2, 5, 1, 4,
          0, 3, 4, 5,
          6, 1, 0, 4,
          5, 3, 2, 7
        };

        // set up vertices-at-cube array
        static const Index v_c1[1*8] =
        {
          1, 2, 6, 7, 4, 5, 0, 3
        };

        // set up edges-at-quad array
        static const Index e_q1[6*4] =
        {
          2, 0, 1, 3,
          1, 4, 9, 11,
          7, 5, 4, 6,
          3, 6, 8, 10,
          11, 8, 0, 5,
          10, 9, 7, 2
        };

        // set up edges-at-cube array
        static const Index e_c1[1*12] =
        {
          4, 1, 6, 3, 11, 9, 8, 10, 5, 7, 0, 2
        };

        // set up quad-at-cube array
        static const Index q_c1[1*6] =
        {
          1, 3, 2, 0, 4, 5
        };

        //// third possibility ////

        //         v_6__________e_5__>______v_5
        //          /|                     /|
        //      e_8^ |        q_5     e_10v |
        //        /__|________<__________/  |
        //    v_2|   |       e_3      v_3|  |
        //       |   |                   |  |
        //       |   |e_7       q_1      |q |e_6
        //       | q v         x         v2 |
        //    e_1| 3 |       q_0      e_2|  |
        //       v   |                   v  |
        //       |   |v_4_______e_4__<___|__|v_1
        //       |  /                    |  /
        //       | ve_11      q_4        | ^e_9
        //       |/___________>__________|/
        //      v_0          e_0        v_7
        //
        // cube orientation: v_0-v_7-v_4-v_1-v_2-v_3-v_6-v_5
        // quad orientation:
        // q_0 : v_2-v_0-v_3-v_7
        // q_1 : v_5-v_1-v_6-v_4
        // q_2 : v_7-v_3-v_1-v_5
        // q_3 : v_6-v_4-v_2-v_0
        // q_4 : v_7-v_0-v_1-v_4
        // q_5 : v_6-v_5-v_2-v_3

        // set up vertex coordinates array
        static const Real vtx2[3*8] =
        {
          0.0, 0.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          0.0, 1.0, 0.0,
          1.0, 1.0, 1.0,
          0.0, 1.0, 1.0,
          1.0, 0.0, 0.0,
        };

        // set up vertices-at-edge array
        static const Index v_e2[12*2] =
        {
          0, 7,
          2, 0,
          3, 7,
          3, 2,

          4, 1,
          5, 6,
          1, 5,
          4, 6,

          6, 2,
          1, 7,
          3, 5,
          0, 4
        };

        // set up vertices-at-quad array
        static const Index v_q2[6*4] =
        {
          2, 0, 3, 7,
          5, 1, 6, 4,
          7, 3, 1, 5,
          6, 4, 2, 0,
          7, 0, 1, 4,
          6, 5, 2, 3
        };

        // set up vertices-at-cube array
        static const Index v_c2[1*8] =
        {
          0, 7, 4, 1, 2, 3, 6, 5
        };

        // set up edges-at-quad array
        static const Index e_q2[6*4] =
        {
          1, 2, 3, 0,
          6, 7, 5, 4,
          2, 6, 9, 10,
          7, 1, 8, 11,
          0, 4, 9, 11,
          5, 3, 8, 10
        };

        // set up edges-at-cube array
        static const Index e_c2[1*12] =
        {
          0, 4, 3, 5, 11, 9, 8, 10, 1, 2, 7, 6
        };

        // set up quad-at-cube array
        static const Index q_c2[1*6] =
        {
          4, 5, 0, 1, 3, 2
        };

        switch(orientation)
        {
          case 0:
            copy_vtx(mesh->get_vertex_set(), vtx0);
            copy_idx(mesh->get_index_set<1,0>(), v_e0);
            copy_idx(mesh->get_index_set<2,0>(), v_q0);
            copy_idx(mesh->get_index_set<3,0>(), v_c0);
            copy_idx(mesh->get_index_set<2,1>(), e_q0);
            copy_idx(mesh->get_index_set<3,1>(), e_c0);
            copy_idx(mesh->get_index_set<3,2>(), q_c0);
            break;
          case 1:
            copy_vtx(mesh->get_vertex_set(), vtx1);
            copy_idx(mesh->get_index_set<1,0>(), v_e1);
            copy_idx(mesh->get_index_set<2,0>(), v_q1);
            copy_idx(mesh->get_index_set<3,0>(), v_c1);
            copy_idx(mesh->get_index_set<2,1>(), e_q1);
            copy_idx(mesh->get_index_set<3,1>(), e_c1);
            copy_idx(mesh->get_index_set<3,2>(), q_c1);
            break;
         case 2:
            copy_vtx(mesh->get_vertex_set(), vtx2);
            copy_idx(mesh->get_index_set<1,0>(), v_e2);
            copy_idx(mesh->get_index_set<2,0>(), v_q2);
            copy_idx(mesh->get_index_set<3,0>(), v_c2);
            copy_idx(mesh->get_index_set<2,1>(), e_q2);
            copy_idx(mesh->get_index_set<3,1>(), e_c2);
            copy_idx(mesh->get_index_set<3,2>(), q_c2);
            break;
        }
        // okay
        return mesh;
      } // create_hexa_mesh_3d

      void validate_refined_hexa_mesh_3d(const HexaMesh& mesh, int orientation)
      {
        // validate sizes
        if(mesh.get_num_entities(0) != 27)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 54)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 36)
          throw String("Quad count mismatch");
        if(mesh.get_num_entities(3) != 8)
          throw String("Cube count mismatch");

        //// first possibility ////

        //         v_6__________e_3___>_____v_7
        //          /|                     /|
        //      e_6^ |        q_1      e_7^ |
        //        /__|________>__________/  |
        //    v_4|   |       e_2      v_5|  |
        //       |   |                   |  |
        //       |   |e_10      q_3      |q |e_11
        //       | q ^         x         |5 ^
        //    e_8| 4 |       q_2      e_9|  |
        //       ^   |                   ^  |
        //       |   |v_2_______e_1__>___|__|v_3
        //       |  /                    |  /
        //       | ^e_4       q_0        | ^e_5
        //       |/___________>__________|/
        //      v_0          e_0        v_1
        //
        // orientation: v_0-v_1-v_2-v_3-v_4-v_5-v_6-v_7
        // quad orientation:
        // q_0 : v_0-v_1-v_2-v_3
        // q_1 : v_5-v_6-v_7-v_8
        // q_2 : v_1-v_2-v_5-v_6
        // q_3 : v_3-v_4-v_7-v_8
        // q_4 : v_1-v_3-v_5-v_7
        // q_5 : v_2-v_4-v_6-v_8


        // vertex coordinates array
        static const Real vtx0[] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          0.0, 1.0, 1.0,
          1.0, 1.0, 1.0,

          0.5, 0.0, 0.0,
          0.5, 1.0, 0.0,
          0.5, 0.0, 1.0,
          0.5, 1.0, 1.0,
          0.0, 0.5, 0.0,
          1.0, 0.5, 0.0,
          0.0, 0.5, 1.0,
          1.0, 0.5, 1.0,
          0.0, 0.0, 0.5,
          1.0, 0.0, 0.5,
          0.0, 1.0, 0.5,
          1.0, 1.0, 0.5,

          0.5, 0.5, 0.0,
          0.5, 0.5, 1.0,
          0.5, 0.0, 0.5,
          0.5, 1.0, 0.5,
          0.0, 0.5, 0.5,
          1.0, 0.5, 0.5,

          0.5, 0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e0[] =
        {
          0, 8,
          8, 1,
          2, 9,
          9, 3,
          4, 10,
          10, 5,
          6, 11,
          11, 7,
          0, 12,
          12, 2,
          1, 13,
          13, 3,
          4, 14,
          14, 6,
          5, 15,
          15, 7,
          0, 16,
          16, 4,
          1, 17,
          17, 5,
          2, 18,
          18, 6,
          3, 19,
          19, 7,

          8, 20,
          20, 9,
          12, 20,
          20, 13,
          10, 21,
          21, 11,
          14, 21,
          21, 15,
          8, 22,
          22, 10,
          16, 22,
          22, 17,
          9, 23,
          23, 11,
          18, 23,
          23, 19,
          12, 24,
          24, 14,
          16, 24,
          24, 18,
          13, 25,
          25, 15,
          17, 25,
          25, 19,

          20, 26,
          26, 21,
          22, 26,
          26, 23,
          24, 26,
          26, 25
        };

        // vertices-at-quad
        static const Index v_q0[] =
        {
          0, 8, 12, 20,
          8, 1, 20, 13,
          12, 20, 2, 9,
          20, 13, 9, 3,
          4, 10, 14, 21,
          10, 5, 21, 15,
          14, 21, 6, 11,
          21, 15, 11, 7,
          0, 8, 16, 22,
          8, 1, 22, 17,
          16, 22, 4, 10,
          22, 17, 10, 5,
          2, 9, 18, 23,
          9, 3, 23, 19,
          18, 23, 6, 11,
          23, 19, 11, 7,
          0, 12, 16, 24,
          12, 2, 24, 18,
          16, 24, 4, 14,
          24, 18, 14, 6,
          1, 13, 17, 25,
          13, 3, 25, 19,
          17, 25, 5, 15,
          25, 19, 15, 7,


          8, 20, 22, 26,
          20, 9, 26, 23,
          22, 26, 10, 21,
          26, 23, 21, 11,
          12, 20, 24, 26,
          20, 13, 26, 25,
          24, 26, 14, 21,
          26, 25, 21, 15,
          16, 22, 24, 26,
          22, 17, 26, 25,
          24, 26, 18, 23,
          26, 25, 23, 19
        };

        // vertices-at-cube
        static const Index v_c0[] =
        {
          0, 8, 12, 20, 16, 22, 24, 26,
          8, 1, 20, 13, 22, 17, 26, 25,
          12, 20, 2, 9, 24, 26, 18, 23,
          20, 13, 9, 3, 26, 25, 23, 19,
          16, 22, 24, 26, 4, 10, 14, 21,
          22, 17, 26, 25, 10, 5, 21, 15,
          24, 26, 18, 23, 14, 21, 6, 11,
          26, 25, 23, 19, 21, 15, 11, 7
        };

        // edges-at-quad
        static const Index e_q0[] =
        {
          0, 26, 8, 24,
          1, 27, 24, 10,
          26, 2, 9, 25,
          27, 3, 25, 11,

          4, 30, 12, 28,
          5, 31, 28, 14,
          30, 6, 13, 29,
          31, 7, 29, 15,

          0, 34, 16, 32,
          1, 35, 32, 18,
          34, 4, 17, 33,
          35, 5, 33, 19,

          2, 38, 20, 36,
          3, 39, 36, 22,
          38, 6, 21, 37,
          39, 7, 37, 23,

          8, 42, 16, 40,
          9, 43, 40, 20,
          42, 12, 17, 41,
          43, 13, 41, 21,

          10, 46, 18, 44,
          11, 47, 44, 22,
          46, 14, 19, 45,
          47, 15, 45, 23,

          24, 50, 32, 48,
          25, 51, 48, 36,
          50, 28, 33, 49,
          51, 29, 49, 37,

          26, 52, 40, 48,
          27, 53, 48, 44,
          52, 30, 41, 49,
          53, 31, 49, 45,

          34, 52, 42, 50,
          35, 53, 50, 46,
          52, 38, 43, 51,
          53, 39, 51, 47
        };

        // edges-at-cube
        static const Index e_c0[] =
        {
          0, 26, 34, 52, 8, 24, 42, 50, 16, 32, 40, 48,
          1, 27, 35, 53, 24, 10, 50, 46, 32, 18, 48, 44,
          26, 2, 52, 38, 9, 25, 43, 51, 40, 48, 20, 36,
          27, 3, 53, 39, 25, 11, 51, 47, 48, 44, 36, 22,
          34, 52, 4, 30, 42, 50, 12, 28, 17, 33, 41, 49,
          35, 53, 5, 31, 50, 46, 28, 14, 33, 19, 49, 45,
          52, 38, 30, 6, 43, 51, 13, 29, 41, 49, 21, 37,
          53, 39, 31, 7, 51, 47, 29, 15, 49, 45, 37, 23
        };

        // quad-at-cube
        static const Index q_c0[] =
        {
          0, 32, 8, 28, 16, 24,
          1, 33, 9, 29, 24, 20,
          2, 34, 28, 12, 17, 25,
          3, 35, 29, 13, 25, 21,
          32, 4, 10, 30, 18, 26,
          33, 5, 11, 31, 26, 22,
          34, 6, 30, 14, 19, 27,
          35, 7, 31, 15, 27, 23
        };

        //// second possibility ////

        //         v_5__________e_7__>______v_2
        //          /|                     /|
        //      e_6^ |        q_2      e_4v |
        //        /__|________<__________/  |
        //    v_4|   |       e_5      v_1|  |
        //       |   |                   |  |
        //       |   |e_10      q_5      |q |e_9
        //       | q v         x         |1 |
        //    e_8| 3 |       q_4     e_11|  |
        //       v   |                   v  |
        //       |   |v_3_______e_2__<___|__|v_7
        //       |  /                    |  /
        //       | ve_3       q_0        | ^e_1
        //       |/___________>__________|/
        //      v_0          e_0        v_6
        //
        // cube orientation: v_1-v_2-v_6-v_7-v_4-v_5-v_0-v_3
        // quad orientation:
        // q_0 : v_7-v_3-v_6-v_0
        // q_1 : v_7-v_6-v_2-v_1
        // q_2 : v_2-v_5-v_1-v_4
        // q_3 : v_0-v_3-v_4-v_5
        // q_4 : v_6-v_1-v_0-v_4
        // q_5 : v_5-v_3-v_2-v_7

        // vertex coordinates array
        static const Real vtx1[] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 1.0,
          1.0, 1.0, 1.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          0.0, 1.0, 1.0,
          1.0, 0.0, 0.0,
          1.0, 1.0, 0.0,

          0.5, 0.0, 0.0,
          1.0, 0.5, 0.0,
          0.5, 1.0, 0.0,
          0.0, 0.5, 0.0,
          1.0, 0.5, 1.0,
          0.5, 0.0, 1.0,
          0.0, 0.5, 1.0,
          0.5, 1.0, 1.0,
          0.0, 0.0, 0.5,
          1.0, 1.0, 0.5,
          0.0, 1.0, 0.5,
          1.0, 0.0, 0.5,

          0.5, 0.5, 0.0,
          1.0, 0.5, 0.5,
          0.5, 0.5, 1.0,
          0.0, 0.5, 0.5,
          0.5, 0.0, 0.5,
          0.5, 1.0, 0.5,

          0.5, 0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e1[] =
        {
          0, 8,
          8, 6,
          6, 9,
          9, 7,
          7, 10,
          10, 3,
          3, 11,
          11, 0,
          2, 12,
          12, 1,
          1, 13,
          13, 4,
          4, 14,
          14, 5,
          5, 15,
          15, 2,
          4, 16,
          16, 0,
          2, 17,
          17, 7,
          5, 18,
          18, 3,
          1, 19,
          19, 6,

          10, 20,
          20, 8,
          9, 20,
          20, 11,
          9, 21,
          21, 12,
          17, 21,
          21, 19,
          15, 22,
          22, 13,
          12, 22,
          22, 14,
          11, 23,
          23, 14,
          16, 23,
          23, 18,
          19, 24,
          24, 16,
          8, 24,
          24, 13,
          18, 25,
          25, 17,
          15, 25,
          25, 10,

          21, 26,
          26, 23,
          22, 26,
          26, 20,
          24, 26,
          26, 25
        };

        // vertices-at-quad
        static const Index v_q1[] =
        {
          7, 10, 9, 20,
          10, 3, 20, 11,
          9, 20, 6, 8,
          20, 11, 8, 0,

          7, 9, 17, 21,
          9, 6, 21, 19,
          17, 21, 2, 12,
          21, 19, 12, 1,

          2, 15, 12, 22,
          15, 5, 22, 14,
          12, 22, 1, 13,
          22, 14, 13, 4,

          0, 11, 16, 23,
          11, 3, 23, 18,
          16, 23, 4, 14,
          23, 18, 14, 5,

          6, 19, 8, 24,
          19, 1, 24, 13,
          8, 24, 0, 16,
          24, 13, 16, 4,

          5, 18, 15, 25,
          18, 3, 25, 10,
          15, 25, 2, 17,
          25, 10, 17, 7,

          12, 21, 22, 26,
          21, 9, 26, 20,
          22, 26, 14, 23,
          26, 20, 23, 11,

          19, 21, 24, 26,
          21, 17, 26, 25,
          24, 26, 16, 23,
          26, 25, 23, 18,

          13, 22, 24, 26,
          22, 15, 26, 25,
          24, 26, 8, 20,
          26, 25, 20, 10
        };

        // vertices-at-cube
        static const Index v_c1[] =
        {
          1, 12, 19, 21, 13, 22, 24, 26,
          12, 2, 21, 17, 22, 15, 26, 25,
          19, 21, 6, 9, 24, 26, 8, 20,
          21, 17, 9, 7, 26, 25, 20, 10,
          13, 22, 24, 26, 4, 14, 16, 23,
          22, 15, 26, 25, 14, 5, 23, 18,
          24, 26, 8, 20, 16, 23, 0, 11,
          26, 25, 20, 10, 23, 18, 11, 3
        };

        // edges-at-quad
        static const Index e_q1[] =
        {
          4, 26, 3, 24,
          5, 27, 24, 6,
          26, 1, 2, 25,
          27, 0, 25, 7,

          3, 30, 19, 28,
          2, 31, 28, 23,
          30, 8, 18, 29,
          31, 9, 29, 22,

          15, 34, 8, 32,
          14, 35, 32, 13,
          34, 10, 9, 33,
          35, 11, 33, 12,

          7, 38, 17, 36,
          6, 39, 36, 21,
          38, 12, 16, 37,
          39, 13, 37, 20,

          23, 42, 1, 40,
          22, 43, 40, 10,
          42, 17, 0, 41,
          43, 16, 41, 11,

          20, 46, 14, 44,
          21, 47, 44, 5,
          46, 18, 15, 45,
          47, 19, 45, 4,

          29, 50, 34, 48,
          28, 51, 48, 26,
          50, 37, 35, 49,
          51, 36, 49, 27,

          31, 52, 40, 48,
          30, 53, 48, 45,
          52, 38, 41, 49,
          53, 39, 49, 44,

          33, 52, 43, 50,
          32, 53, 50, 46,
          52, 25, 42, 51,
          53, 24, 51, 47
        };

        // edges-at-cube
        static const Index e_c1[] =
        {
          9, 31, 33, 52, 22, 29, 43, 50, 10, 34, 40, 48,
          8, 30, 32, 53, 29, 18, 50, 46, 34, 15, 48, 45,
          31, 2, 52, 25, 23, 28, 42, 51, 40, 48, 1, 26,
          30, 3, 53, 24, 28, 19, 51, 47, 48, 45, 26, 4,

          33, 52, 12, 38, 43, 50, 16, 37, 11, 35, 41, 49,
          32, 53, 13, 39, 50, 46, 37, 20, 35, 14, 49, 44,
          52, 25, 38, 7, 42, 51, 17, 36, 41, 49, 0, 27,
          53, 24, 39, 6, 51, 47, 36, 21, 49, 44, 27, 5
        };

        // quad-at-cube
        static const Index q_c1[] =
        {
          7, 32, 10, 28, 17, 24,
          6, 33, 8, 29, 24, 22,
          5, 34, 28, 2, 16, 25,
          4, 35, 29, 0, 25, 23,
          32, 14, 11, 30, 19, 26,
          33, 15, 9, 31, 26, 20,
          34, 12, 30, 3, 18, 27,
          35, 13, 31, 1, 27, 21
        };

        //// third possibility ////

        //         v_6__________e_5__>______v_5
        //          /|                     /|
        //      e_8^ |        q_5     e_10v |
        //        /__|________<__________/  |
        //    v_2|   |       e_3      v_3|  |
        //       |   |                   |  |
        //       |   |e_7       q_1      |q |e_6
        //       | q v         x         v2 |
        //    e_1| 3 |       q_0      e_2|  |
        //       v   |                   v  |
        //       |   |v_4_______e_4__<___|__|v_1
        //       |  /                    |  /
        //       | ve_11      q_4        | ^e_9
        //       |/___________>__________|/
        //      v_0          e_0        v_7
        //
        // cube orientation: v_0-v_7-v_4-v_1-v_2-v_3-v_6-v_5
        // quad orientation:
        // q_0 : v_2-v_0-v_3-v_7
        // q_1 : v_5-v_1-v_6-v_4
        // q_2 : v_7-v_3-v_1-v_5
        // q_3 : v_6-v_4-v_2-v_0
        // q_4 : v_7-v_0-v_1-v_4
        // q_5 : v_6-v_5-v_2-v_3

        // vertex coordinates array
        static const Real vtx2[] =
        {
          0.0, 0.0, 0.0,
          1.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          1.0, 0.0, 1.0,
          0.0, 1.0, 0.0,
          1.0, 1.0, 1.0,
          0.0, 1.0, 1.0,
          1.0, 0.0, 0.0,

          0.5, 0.0, 0.0,
          0.0, 0.0, 0.5,
          1.0, 0.0, 0.5,
          0.5, 0.0, 1.0,
          0.5, 1.0, 0.0,
          0.5, 1.0, 1.0,
          1.0, 1.0, 0.5,
          0.0, 1.0, 0.5,
          0.0, 0.5, 1.0,
          1.0, 0.5, 0.0,
          1.0, 0.5, 1.0,
          0.0, 0.5, 0.0,

          0.5, 0.0, 0.5,
          0.5, 1.0, 0.5,
          1.0, 0.5, 0.5,
          0.0, 0.5, 0.5,
          0.5, 0.5, 0.0,
          0.5, 0.5, 1.0,

          0.5, 0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e2[] =
        {
          0, 8,
          8, 7,
          2, 9,
          9, 0,
          3, 10,
          10, 7,
          3, 11,
          11, 2,
          4, 12,
          12, 1,
          5, 13,
          13, 6,
          1, 14,
          14, 5,
          4, 15,
          15, 6,
          6, 16,
          16, 2,
          1, 17,
          17, 7,
          3, 18,
          18, 5,
          0, 19,
          19, 4,

          9, 20,
          20, 10,
          11, 20,
          20, 8,
          14, 21,
          21, 15,
          13, 21,
          21, 12,
          10, 22,
          22, 14,
          17, 22,
          22, 18,
          15, 23,
          23, 9,
          16, 23,
          23, 19,
          8, 24,
          24, 12,
          17, 24,
          24, 19,
          13, 25,
          25, 11,
          16, 25,
          25, 18,

          24, 26,
          26, 25,
          20, 26,
          26, 21,
          23, 26,
          26, 22
        };

        // vertices-at-quad
        static const Index v_q2[] =
        {
          2, 9, 11, 20,
          9, 0, 20, 8,
          11, 20, 3, 10,
          20, 8, 10, 7,

          5, 14, 13, 21,
          14, 1, 21, 12,
          13, 21, 6, 15,
          21, 12, 15, 4,

          7, 10, 17, 22,
          10, 3, 22, 18,
          17, 22, 1, 14,
          22, 18, 14, 5,

          6, 15, 16, 23,
          15, 4, 23, 19,
          16, 23, 2, 9,
          23, 19, 9, 0,

          7, 8, 17, 24,
          8, 0, 24, 19,
          17, 24, 1, 12,
          24, 19, 12, 4,

          6, 13, 16, 25,
          13, 5, 25, 18,
          16, 25, 2, 11,
          25, 18, 11, 3,

          8, 24, 20, 26,
          24, 12, 26, 21,
          20, 26, 11, 25,
          26, 21, 25, 13,

          19, 24, 23, 26,
          24, 17, 26, 22,
          23, 26, 16, 25,
          26, 22, 25, 18,

          9, 20, 23, 26,
          20, 10, 26, 22,
          23, 26, 15, 21,
          26, 22, 21, 14
        };

        // vertices-at-cube
        static const Index v_c2[] =
        {
          0, 8, 19, 24, 9, 20, 23, 26,
          8, 7, 24, 17, 20, 10, 26, 22,
          19, 24, 4, 12, 23, 26, 15, 21,
          24, 17, 12, 1, 26, 22, 21, 14,
          9, 20, 23, 26, 2, 11, 16, 25,
          20, 10, 26, 22, 11, 3, 25, 18,
          23, 26, 15, 21, 16, 25, 6, 13,
          26, 22, 21, 14, 25, 18, 13, 5
        };

        // edges-at-quad
        static const Index e_q2[] =
        {
          2, 26, 7, 24,
          3, 27, 24, 0,
          26, 4, 6, 25,
          27, 5, 25, 1,

          13, 30, 10, 28,
          12, 31, 28, 9,
          30, 15, 11, 29,
          31, 14, 29, 8,

          5, 34, 19, 32,
          4, 35, 32, 20,
          34, 12, 18, 33,
          35, 13, 33, 21,

          15, 38, 16, 36,
          14, 39, 36, 23,
          38, 2, 17, 37,
          39, 3, 37, 22,

          1, 42, 19, 40,
          0, 43, 40, 22,
          42, 9, 18, 41,
          43, 8, 41, 23,

          11, 46, 16, 44,
          10, 47, 44, 21,
          46, 7, 17, 45,
          47, 6, 45, 20,

          40, 50, 27, 48,
          41, 51, 48, 31,
          50, 45, 26, 49,
          51, 44, 49, 30,

          43, 52, 39, 48,
          42, 53, 48, 34,
          52, 46, 38, 49,
          53, 47, 49, 35,

          24, 52, 37, 50,
          25, 53, 50, 32,
          52, 29, 36, 51,
          53, 28, 51, 33
        };

        // edges-at-cube
        static const Index e_c2[] =
        {
          0, 43, 24, 52, 22, 40, 37, 50, 3, 27, 39, 48,
          1, 42, 25, 53, 40, 19, 50, 32, 27, 5, 48, 34,
          43, 8, 52, 29, 23, 41, 36, 51, 39, 48, 14, 31,
          42, 9, 53, 28, 41, 18, 51, 33, 48, 34, 31, 12,

          24, 52, 7, 46, 37, 50, 17, 45, 2, 26, 38, 49,
          25, 53, 6, 47, 50, 32, 45, 20, 26, 4, 49, 35,
          52, 29, 46, 11, 36, 51, 16, 44, 38, 49, 15, 30,
          53, 28, 47, 10, 51, 33, 44, 21, 49, 35, 30, 13
        };

        // quad-at-cube
        static const Index q_c2[] =
        {
          17, 32, 1, 28, 15, 24,
          16, 33, 3, 29, 24, 8,
          19, 34, 28, 7, 13, 25,
          18, 35, 29, 5, 25, 10,
          32, 22, 0, 30, 14, 26,
          33, 23, 2, 31, 26, 9,
          34, 20, 30, 6, 12, 27,
          35, 21, 31, 4, 27, 11
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

            // check vertices-at-cube
            if(!comp_idx(mesh.get_index_set<3,0>(), v_c0))
              throw String("Vertex-At-Cube index set refinement failure");

            // check edges-at-quad
            if(!comp_idx(mesh.get_index_set<2,1>(), e_q0))
              throw String("Edge-At-Quad index set refinement failure");

            // check edges-at-cube
            if(!comp_idx(mesh.get_index_set<3,1>(), e_c0))
              throw String("Edge-At-Cube index set refinement failure");

            // check quad-at-cube
            if(!comp_idx(mesh.get_index_set<3,2>(), q_c0))
              throw String("Quad-At-Cube index set refinement failure");
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

            // check vertices-at-cube
            if(!comp_idx(mesh.get_index_set<3,0>(), v_c1))
              throw String("Vertex-At-Cube index set refinement failure");

            // check edges-at-quad
            if(!comp_idx(mesh.get_index_set<2,1>(), e_q1))
              throw String("Edge-At-Quad index set refinement failure");

            // check edges-at-cube
            if(!comp_idx(mesh.get_index_set<3,1>(), e_c1))
              throw String("Edge-At-Cube index set refinement failure");

            // check quad-at-cube
            if(!comp_idx(mesh.get_index_set<3,2>(), q_c1))
              throw String("Quad-At-Cube index set refinement failure");
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

            // check vertices-at-cube
            if(!comp_idx(mesh.get_index_set<3,0>(), v_c2))
              throw String("Vertex-At-Cube index set refinement failure");

            // check edges-at-quad
            if(!comp_idx(mesh.get_index_set<2,1>(), e_q2))
              throw String("Edge-At-Quad index set refinement failure");

            // check edges-at-cube
            if(!comp_idx(mesh.get_index_set<3,1>(), e_c2))
              throw String("Edge-At-Cube index set refinement failure");

            // check quad-at-cube
            if(!comp_idx(mesh.get_index_set<3,2>(), q_c2))
              throw String("Quad-At-Cube index set refinement failure");
            break;

        } //switch
      } // validate_refined_hexa_mesh_3d

    } // namespace TestAux
  } // namespace Geometry
} // namespace FEAT

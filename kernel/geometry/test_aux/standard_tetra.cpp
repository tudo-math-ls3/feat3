// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/test_aux/standard_tetra.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

namespace FEAT
{
  namespace Geometry
  {
    namespace TestAux
    {
      TetraMesh* create_tetra_mesh_3d(int orientation)
      {

        Index num_entities[] =
        {
          4, // vertices
          6, // edges
          4, // triangles
          1  // tetrahedron
        };

        // create mesh
        TetraMesh* mesh = new TetraMesh(num_entities);

        // first possibility (standard)

        // set up vertex coordinates array
        static const Real vtx0[3*4] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0
        };

        // set up vertices-at-edge array
        static const Index v_e0[6*2] =
        {
          0, 1,
          0, 2,
          0, 3,
          1, 2,
          1, 3,
          2, 3
        };

        // set up vertices-at-triangle array
        static const Index v_t0[4*3] =
        {
          1, 2, 3,
          0, 2, 3,
          0, 1, 3,
          0, 1, 2
        };

        // set up vertices-at-tetrahedron array
        static const Index v_s0[1*4] =
        {
          0, 1, 2, 3
        };

        // set up edges-at-triangle array
        static const Index e_t0[3*4] =
        {
          5, 4, 3,
          5, 2, 1,
          4, 2, 0,
          3, 1, 0
        };

        // set up edges-at-tetrahedron array
        static const Index e_s0[1*6] =
        {
          0, 1, 2, 3, 4, 5
        };

        // set up triangle-at-tetrahedron array
        static const Index t_s0[1*4] =
        {
          0, 1, 2, 3
        };

        // second possibility

        /*                        v_3
                                  /\\
                                 /  \  \
                                /    \    \
                               /      \      \
                              /        \        \
                             /          \          \
                            /            \           v e_5
                           /              \              \
                          /                \                \
                         /                  \                  \
                        /                    \                    \
                       /                  e_3 \                      \ v_2
                  e_4 v                        v                 /    /
                     /                          \           /        /
                    /                            \      /           /
                   /                              \ /              /
                  /                       e_0  /   \              /
                 /              s_0        ^        \            /
                /                    /               \          ^ e_1
               /                /                     \        /
              /            /                           \      /
             /        /                                 \    /
            /    /                                       \  /
           //____________________<________________________\/
          v_0                    e_2                      v_1

         orientation:
            t0: v_2 - v_1 - v_0
            t1: v_2 - v_3 - v_1
            t2: v_3 - v_0 - v_2
            t3: v_0 - v_3 - v_1
         tetrahedron:
                v_0 - v_1 - v_2 - v_3
        */


        // set up vertex coordinates array
        static const Real vtx1[3*4] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0
        };

        // set up vertices-at-edge array
        static const Index v_e1[6*2] =
        {
          0, 2,
          1, 2,
          1, 0,
          3, 1,
          3, 0,
          3, 2
        };

        // set up vertices-at-triangle array
        static const Index v_t1[4*3] =
        {
          2, 1, 0,
          2, 3, 1,
          3, 0, 2,
          0, 3, 1
        };

        // set up vertices-at-tetrahedron array
        static const Index v_s1[1*4] =
        {
          0, 1, 2, 3
        };

        // set up edges-at-triangle array
        static const Index e_t1[3*4] =
        {
          2, 0, 1,
          3, 1, 5,
          0, 5, 4,
          3, 2, 4
        };

        // set up edges-at-tetrahedron array
        static const Index e_s1[1*6] =
        {
          2, 0, 4, 1, 3, 5
        };

        // set up triangle-at-tetrahedron array
        static const Index t_s1[1*4] =
        {
          1, 2, 3, 0
        };

        // third possibility

        /*                        v_2
                                  /\\
                                 /  \  \
                                /    \    \
                               /      \      \
                              /        \        \
                             /          \          \
                            /            \           ^ e_5
                           /              \              \
                          /                \                \
                         /                  \                  \
                        /                    \                    \
                       /                  e_3 \                      \ v_0
                  e_0 v                        v                 /    /
                     /                          \           /        /
                    /                            \      /           /
                   /                              \ /              /
                  /                       e_4  /   \              /
                 /              s_0        ^        \            /
                /                    /               \          ^ e_2
               /                /                     \        /
              /            /                           \      /
             /        /                                 \    /
            /    /                                       \  /
           //____________________>________________________\/
          v_3                    e_1                      v_1

         orientation:
            t0: v_0 - v_3 - v_2
            t1: v_0 - v_1 - v_3
            t2: v_0 - v_1 - v_2
            t3: v_2 - v_1 - v_3
         tetrahedron:
                v_3 - v_1 - v_0 - v_2
        */

        // set up vertex coordinates array
        static const Real vtx2[3*4] =
        {
          1.0, 2.0, 0.0,
          2.0, 1.0, 0.0,
          0.0, 0.0, 2.0,
          1.0, 1.0, 0.0
        };

        // set up vertices-at-edge array
        static const Index v_e2[6*2] =
        {
          2, 3,
          3, 1,
          1, 0,
          2, 1,
          3, 0,
          0, 2
        };

        // set up vertices-at-triangle array
        static const Index v_t2[4*3] =
        {
          0, 3, 2,
          0, 1, 3,
          0, 1, 2,
          2, 1, 3
        };

        // set up vertices-at-tetrahedron array
        static const Index v_s2[1*4] =
        {
          3, 1, 0, 2
        };

        // set up edges-at-triangle array
        static const Index e_t2[3*4] =
        {
          0, 5, 4,
          1, 4, 2,
          3, 5, 2,
          1, 0, 3
        };

        // set up edges-at-tetrahedron array
        static const Index e_s2[1*6] =
        {
          1, 4, 0, 2, 3, 5
        };

        // set up triangle-at-tetrahedron array
        static const Index t_s2[1*4] =
        {
          2, 0, 3, 1
        };

        switch(orientation)
        {
          case 0:
            copy_vtx(mesh->get_vertex_set(), vtx0);
            copy_idx(mesh->get_index_set<1,0>(), v_e0);
            copy_idx(mesh->get_index_set<2,0>(), v_t0);
            copy_idx(mesh->get_index_set<2,1>(), e_t0);
            copy_idx(mesh->get_index_set<3,0>(), v_s0);
            copy_idx(mesh->get_index_set<3,1>(), e_s0);
            copy_idx(mesh->get_index_set<3,2>(), t_s0);
            break;

          case 1:
            copy_vtx(mesh->get_vertex_set(), vtx1);
            copy_idx(mesh->get_index_set<1,0>(), v_e1);
            copy_idx(mesh->get_index_set<2,0>(), v_t1);
            copy_idx(mesh->get_index_set<2,1>(), e_t1);
            copy_idx(mesh->get_index_set<3,0>(), v_s1);
            copy_idx(mesh->get_index_set<3,1>(), e_s1);
            copy_idx(mesh->get_index_set<3,2>(), t_s1);
            break;

          case 2:
            copy_vtx(mesh->get_vertex_set(), vtx2);
            copy_idx(mesh->get_index_set<1,0>(), v_e2);
            copy_idx(mesh->get_index_set<2,0>(), v_t2);
            copy_idx(mesh->get_index_set<2,1>(), e_t2);
            copy_idx(mesh->get_index_set<3,0>(), v_s2);
            copy_idx(mesh->get_index_set<3,1>(), e_s2);
            copy_idx(mesh->get_index_set<3,2>(), t_s2);
            break;

          default:
            XABORTM("Unhandled orientation "+stringify(orientation));
        }
        // okay
        return mesh;
      } // create_tetra_mesh_3d

      void validate_refined_tetra_mesh_3d(const TetraMesh& mesh, int orientation)
      {

        // validate sizes
        if(mesh.get_num_entities(0) != 11)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 30)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 32)
          throw String("Triangle count mismatch");
        if(mesh.get_num_entities(3) != 12)
          throw String("Triangle count mismatch");


        // first possibility (standard)

        // vertex coordinates array
        static const Real vtx0[] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0,
          0.5, 0.0, 0.0,
          0.0, 0.5, 0.0,
          0.0, 0.0, 0.5,
          0.5, 0.5, 0.0,
          0.5, 0.0, 0.5,
          0.0, 0.5, 0.5,
          0.25, 0.25, 0.25
        };

        // vertices-at-edge array
        static const Index v_e0[] =
        {
          //edges
          0, 4,
          4, 1,
          0, 5,
          5, 2,
          0, 6,
          6, 3,
          1, 7,
          7, 2,
          1, 8,
          8, 3,
          2, 9,
          9, 3,

          //triangles
          7, 8,
          9, 7,
          8, 9,
          5, 6,
          9, 5,
          6, 9,
          4, 6,
          8, 4,
          6, 8,
          4, 5,
          7, 4,
          5, 7,

          //tetrahedron
          4, 10,
          5, 10,
          6, 10,
          7, 10,
          8, 10,
          9, 10
        };

        // vertices-at-triangle
        static const Index v_t0[] =
        {
          //triangles
          1, 7, 8,
          7, 2, 9,
          8, 9, 3,
          9, 8, 7,
          0, 5, 6,
          5, 2, 9,
          6, 9, 3,
          9, 6, 5,
          0, 4, 6,
          4, 1, 8,
          6, 8, 3,
          8, 6, 4,
          0, 4, 5,
          4, 1, 7,
          5, 7, 2,
          7, 5, 4,
          4, 5, 6,
          4, 7, 8,
          5, 7, 9,
          6, 8, 9,
          4, 5, 10,
          4, 6, 10,
          4, 7, 10,
          4, 8, 10,
          5, 6, 10,
          5, 7, 10,
          5, 9, 10,
          6, 8, 10,
          6, 9, 10,
          7, 8, 10,
          7, 9, 10,
          8, 9, 10
        };

        // vertices-at-tetrahedron
        static const Index v_s0[] =
        {
          4, 6, 5, 0,
          4, 5, 6, 10,
          4, 7, 8, 1,
          4, 8, 7, 10,
          5, 9, 7, 2,
          5, 7, 9, 10,
          6, 8, 9, 3,
          6, 9, 8, 10,

          7, 8, 9, 10,
          5, 9, 6, 10,
          4, 6, 8, 10,
          4, 7, 5, 10

        };

        // edges-at-triangle
        static const Index e_t0[] =
        {
          // surfaces
          12, 8, 6, //0
          10, 13, 7,
          11, 9, 14,
          12, 13, 14,

          15, 4, 2,
          10, 16, 3, //5
          11, 5, 17,
          15, 16, 17,

          18, 4, 0,
          8, 19, 1,
          9, 5, 20, //10
          18, 19, 20,

          21, 2, 0,
          6, 22, 1,
          7, 3, 23,
          21, 22, 23, //15

          // inside
          15, 18, 21,
          12, 19, 22,
          13, 16, 23,
          14, 17, 20,

          25, 24, 21, //20
          26, 24, 18,
          27, 24, 22,
          28, 24, 19,

          26, 25, 15,
          27, 25, 23, //25
          29, 25, 16,

          28, 26, 20,
          29, 26, 17,

          28, 27, 12,
          29, 27, 13, //30

          29, 28, 14
        };

        // edges-at-tetrahedron
        static const Index e_s0[] =
        {
          //corners
          18, 21, 0, 15, 4, 2,
          21, 18, 24, 15, 25, 26,

          22, 19, 1, 12, 6, 8,
          19, 22, 24, 12, 28, 27,

          16, 23, 3, 13, 10, 7,
          23, 16, 25, 13, 27, 29,

          20, 17, 5, 14, 9, 11,
          17, 20, 26, 14, 29, 28,

          // surfaces
          12, 13, 27, 14, 28, 29,
          16, 15, 25, 17, 29, 26,
          18, 19, 24, 20, 26, 28,
          22, 21, 24, 23, 27, 25,
        };

        // triangle-at-tetrahedron
        static const Index t_s0[] =
        {
          4, 12, 8, 16,
          24, 21, 20, 16,

          0, 9, 13, 17,
          29, 22, 23, 17,

          1, 14, 5, 18,
          30, 26, 25, 18,

          2, 6, 10, 19,
          31, 27, 28, 19,

          31, 30, 29, 3,
          28, 24, 26, 7,
          27, 23, 21, 11,
          25, 20, 22, 15

        };

        // second possibility

        // vertex coordinates array
        static const Real vtx1[] =
        {
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0,

          0.0, 0.5, 0.0,
          0.5, 0.5, 0.0,
          0.5, 0.0, 0.0,
          0.5, 0.0, 0.5,
          0.0, 0.0, 0.5,
          0.0, 0.5, 0.5,

          0.25, 0.25, 0.25
        };

        // vertices-at-edge array
        static const Index v_e1[] =
        {
          //edges
          0, 4,
          4, 2,
          1, 5,
          5, 2,
          1, 6,
          6, 0, //5
          3, 7,
          7, 1,
          3, 8,
          8, 0,
          3, 9, //10
          9, 2,

          //triangles
          5, 4,
          6, 5,
          4, 6,

          9, 5, //15
          7, 9,
          5, 7,

          8, 9,
          4, 8,
          9, 4, //20

          8, 6,
          7, 8,
          6, 7,

          //tetrahedron
          6, 10,
          4, 10, //25
          8, 10,
          5, 10,
          7, 10,
          9, 10
        };

        // vertices-at-triangle
        static const Index v_t1[] =
        {
          //triangles
          2, 5, 4,
          5, 1, 6,
          4, 6, 0,
          6, 4, 5,

          2, 9, 5,
          9, 3, 7,//5
          5, 7, 1,
          7, 5, 9,

          3, 8, 9,
          8, 0, 4,
          9, 4, 2,//10
          4, 9, 8,

          0, 8, 6,
          8, 3, 7,
          6, 7, 1,
          7, 6, 8,//15

          //tetrahedron
          6, 4, 8,
          6, 5, 7,
          4, 5, 9,
          8, 7, 9,

          6, 4, 10, //20
          6, 8, 10,
          6, 5, 10,
          6, 7, 10,

          4, 8, 10,
          4, 5, 10,//25
          4, 9, 10,

          8, 7, 10,
          8, 9, 10,

          5, 7, 10,
          5, 9, 10,//30

          7, 9, 10
        };

        // vertices-at-tetrahedron
        static const Index v_s1[] =
        {
          6, 8, 4, 0,
          6, 4, 8, 10,
          6, 5, 7, 1,
          6, 7, 5, 10,
          4, 9, 5, 2,
          4, 5, 9, 10,
          8, 7, 9, 3,
          8, 9, 7, 10,

          5, 7, 9, 10,
          4, 9, 8, 10,
          6, 8, 7, 10,
          6, 5, 4, 10
        };

        // edges-at-triangle
        static const Index e_t1[] =
        {
          // surfaces
          12, 1, 3,
          4, 13, 2,
          5, 0, 14,
          12, 13, 14,

          15, 3, 11,
          6, 16, 10,
          7, 2, 17,
          15, 16, 17,

          18, 10, 8,
          0, 19, 9,
          1, 11, 20,
          18, 19, 20,

          21, 5, 9,
          6, 22, 8,
          7, 4, 23,
          21, 22, 23,

          // inside
          19, 21, 14,
          17, 23, 13,
          15, 20, 12,
          16, 18, 22,

          25, 24, 14,
          26, 24, 21,
          27, 24, 13,
          28, 24, 23,

          26, 25, 19,
          27, 25, 12,
          29, 25, 20,

          28, 26, 22,
          29, 26, 18,

          28, 27, 17,
          29, 27, 15,

          29, 28, 16
        };

        // edges-at-tetrahedron
        static const Index e_s1[] =
        {
          //corners
          21, 14, 5, 19, 9, 0,
          14, 21, 24, 19, 25, 26,
          13, 23, 4, 17, 2, 7,
          23, 13, 24, 17, 28, 27,
          20, 12, 1, 15, 11, 3,
          12, 20, 25, 15, 27, 29,
          22, 18, 8, 16, 6, 10,
          18, 22, 26, 16, 29, 28,

          // surfaces
          17, 15, 27, 16, 28, 29,
          20, 19, 25, 18, 29, 26,
          21, 23, 24, 22, 26, 28,
          13, 14, 24, 12, 27, 25
        };

        // triangle-at-tetrahedron
        static const Index t_s1[] =
        {
          9, 2, 12, 16,
          24, 21, 20, 16,

          6, 14, 1, 17,
          29, 22, 23, 17,

          4, 0, 10, 18,
          30, 26, 25, 18,

          5, 8, 13, 19,
          31, 27, 28, 19,

          31, 30, 29, 7,
          28, 24, 26, 11,
          27, 23, 21, 15,
          25, 20, 22, 3
        };

        // third possibility

        // vertex coordinates array
        static const Real vtx2[] =
        {
          1.0, 2.0, 0.0,
          2.0, 1.0, 0.0,
          0.0, 0.0, 2.0,
          1.0, 1.0, 0.0,

          0.5, 0.5, 1.0,
          1.5, 1.0, 0.0,
          1.5, 1.5, 0.0,
          1.0, 0.5, 1.0,
          1.0, 1.5, 0.0,
          0.5, 1.0, 1.0,

          1.0, 1.0, 0.5
        };

        // vertices-at-edge array
        static const Index v_e2[] =
        {
          //edges
          2, 4,
          4, 3,
          3, 5,
          5, 1,
          1, 6,
          6, 0, //5
          2, 7,
          7, 1,
          3, 8,
          8, 0,
          0, 9, //10
          9, 2,

          //triangles
          8, 9,
          4, 8,
          9, 4,

          6, 8, //15
          5, 6,
          8, 5,

          6, 9,
          7, 6,
          9, 7, //20

          7, 4,
          5, 7,
          4, 5,

          //tetrahedron
          5, 10,
          8, 10, //25
          4, 10,
          6, 10,
          7, 10,
          9, 10
        };

        // vertices-at-triangle
        static const Index v_t2[] =
        {
          //triangles
          0, 8, 9, //0
          8, 3, 4,
          9, 4, 2,
          4, 9, 8,

          0, 6, 8,
          6, 1, 5, //5
          8, 5, 3,
          5, 8, 6,

          0, 6, 9,
          6, 1, 7,
          9, 7, 2,//10
          7, 9, 6,

          2, 7, 4,
          7, 1, 5,
          4, 5, 3,
          5, 4, 7,//15

          //tetrahedron
          5, 8, 4,
          5, 6, 7,
          8, 6, 9,
          4, 7, 9,

          5, 8, 10, //20
          5, 4, 10,
          5, 6, 10,
          5, 7, 10,

          8, 4, 10,
          8, 6, 10,//25
          8, 9, 10,

          4, 7, 10,
          4, 9, 10,

          6, 7, 10,
          6, 9, 10,//30

          7, 9, 10
        };

        // vertices-at-tetrahedron
        static const Index v_s2[] =
        {
          5, 4, 8, 3,
          5, 8, 4, 10,
          5, 6, 7, 1,
          5, 7, 6, 10,
          8, 9, 6, 0,
          8, 6, 9, 10,
          4, 7, 9, 2,
          4, 9, 7, 10,

          6, 7, 9, 10,
          8, 9, 4, 10,
          5, 4, 7, 10,
          5, 6, 8, 10
        };

        // edges-at-triangle
        static const Index e_t2[] =
        {
          // surfaces
          12, 10, 9,
          1, 13, 8,
          0, 11, 14,
          12, 13, 14,

          15, 9, 5,
          3, 16, 4,
          2, 8, 17,
          15, 16, 17,

          18, 10, 5,
          7, 19, 4,
          6, 11, 20,
          18, 19, 20,

          21, 0, 6,
          3, 22, 7,
          2, 1, 23,
          21, 22, 23,

          // inside
          13, 23, 17,
          19, 22, 16,
          18, 12, 15,
          20, 14, 21,

          25, 24, 17,
          26, 24, 23,
          27, 24, 16,
          28, 24, 22,

          26, 25, 13,
          27, 25, 15,
          29, 25, 12,

          28, 26, 21,
          29, 26, 14,

          28, 27, 19,
          29, 27, 18,

          29, 28, 20
        };

        // edges-at-tetrahedron
        static const Index e_s2[] =
        {
          //corners
          23, 17, 2,  13, 1, 8,
          17, 23, 24, 13, 25, 26,
          16, 22, 3, 19, 4, 7,
          22, 16, 24, 19, 28, 27,
          12, 15, 9, 18, 10, 5,
          15, 12, 25, 18, 27, 29,
          21, 14, 0, 20, 6, 11,
          14, 21, 26, 20, 29, 28,

          // surfaces

          19, 18, 27, 20, 28, 29,
          12, 13, 25, 14, 29, 26,
          23, 22, 24, 21, 26, 28,
          16, 17, 24, 15, 27, 25
        };

        // triangle-at-tetrahedron
        static const Index t_s2[] =
        {
          1, 6, 14, 16,
          24, 21, 20, 16,

          9, 13, 5, 17,
          29, 22, 23, 17,

          8, 4, 0, 18,
          30, 26, 25, 18,

          10, 2, 12, 19,
          31, 27, 28, 19,

          31, 30, 29, 11,

          28, 24, 26, 3,

          27, 23, 21, 15,

          25, 20, 22, 7
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

            // check vertices-at-triangle
            if(!comp_idx(mesh.get_index_set<2,0>(), v_t0))
              throw String("Vertex-At-Triangle index set refinement failure");

            // check vertices-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,0>(), v_s0))
              throw String("Vertex-At-Tetrahedron index set refinement failure");

            // check edges-at-triangle
            if(!comp_idx(mesh.get_index_set<2,1>(), e_t0))
              throw String("Edge-At-Triangle index set refinement failure");

            // check edges-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,1>(), e_s0))
              throw String("Edge-At-Tetrahedron index set refinement failure");

            // check triangles-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,2>(), t_s0))
              throw String("Triangle-At-Tetrahedron index set refinement failure");
            break;

          case 1:
            // check vertex coordinates array
            if(!comp_vtx(mesh.get_vertex_set(), vtx1))
              throw String("Vertex coordinate refinement failure");

            // check vertices-at-edge array
            if(!comp_idx(mesh.get_index_set<1,0>(), v_e1))
              throw String("Vertex-At-Edge index set refinement failure");

            // check vertices-at-triangle
            if(!comp_idx(mesh.get_index_set<2,0>(), v_t1))
              throw String("Vertex-At-Triangle index set refinement failure");

            // check vertices-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,0>(), v_s1))
              throw String("Vertex-At-Tetrahedron index set refinement failure");

            // check edges-at-triangle
            if(!comp_idx(mesh.get_index_set<2,1>(), e_t1))
              throw String("Edge-At-Triangle index set refinement failure");

            // check edges-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,1>(), e_s1))
              throw String("Edge-At-Tetrahedron index set refinement failure");

            // check triangles-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,2>(), t_s1))
              throw String("Triangle-At-Tetrahedron index set refinement failure");
            break;

          case 2:
            // check vertex coordinates array
            if(!comp_vtx(mesh.get_vertex_set(), vtx2))
              throw String("Vertex coordinate refinement failure");

            // check vertices-at-edge array
            if(!comp_idx(mesh.get_index_set<1,0>(), v_e2))
              throw String("Vertex-At-Edge index set refinement failure");

            // check vertices-at-triangle
            if(!comp_idx(mesh.get_index_set<2,0>(), v_t2))
              throw String("Vertex-At-Triangle index set refinement failure");

            // check vertices-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,0>(), v_s2))
              throw String("Vertex-At-Tetrahedron index set refinement failure");

            // check edges-at-triangle
            if(!comp_idx(mesh.get_index_set<2,1>(), e_t2))
              throw String("Edge-At-Triangle index set refinement failure");

            // check edges-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,1>(), e_s2))
              throw String("Edge-At-Tetrahedron index set refinement failure");

            // check triangles-at-tetrahedron
            if(!comp_idx(mesh.get_index_set<3,2>(), t_s2))
              throw String("Triangle-At-Tetrahedron index set refinement failure");
            break;

          default:
            XABORTM("Unhandled orientation "+stringify(orientation));
        } //switch
      } // validate_refined_tetra_mesh_3d

      TetraMesh* create_big_tetra_mesh_3d()
      {

        Index num_entities[] =
        {
          5, // vertices
          10, // edges
          9, // triangles
          3  // tetrahedron
        };

        // create mesh
        TetraMesh* mesh = new TetraMesh(num_entities);

        // set up vertex coordinates array
        static const Real vtx0[3*5] =
        {
          -2.0, -4.0,  0.0,
          -2.0,  4.0,  0.0,
           4.0,  0.0,  0.0,
           0.0,  0.0,  42.0,
           0.0,  0.0,  0.0
        };

        // set up vertices-at-edge array
        static const Index v_e0[10*2] =
        {
          3, 4,
          0, 4,
          1, 4,
          2, 4,
          3, 0,
          3, 1,
          3, 2,
          0, 1,
          1, 2,
          2, 0
        };

        // set up vertices-at-triangle array
        static const Index v_t0[9*3] =
        {
          1, 4, 3,//0
          1, 2, 3,
          1, 0, 3,
          3, 0, 2,
          1, 2, 4,
          2, 0, 4,//5
          1, 0, 4,
          2, 3, 4,
          0, 4, 3
        };

        // set up vertices-at-tetrahedron array
        static const Index v_s0[3*4] =
        {
          0, 1, 3, 4,
          2, 3, 4, 1,
          2, 3, 0, 4
        };

        // set up edges-at-triangle array
        static const Index e_t0[3*9] =
        {
          0, 5, 2,
          6, 5, 8,
          4, 5, 7,
          9, 6, 4,
          3, 2, 8,
          1, 3, 9,
          1, 2, 7,
          0, 3, 6,
          0, 4, 1
        };

        // set up edges-at-tetrahedron array
        static const Index e_s0[3*6] =
        {
          7, 4, 1, 5, 2, 0,
          6, 3, 8, 0, 5, 2,
          6, 9, 3, 4, 0, 1
        };

        // set up triangle-at-tetrahedron array
        static const Index t_s0[3*4] =
        {
          0, 8, 6, 2,
          0, 4, 1, 7,
          8, 5, 7, 3
        };

        copy_vtx(mesh->get_vertex_set(), vtx0);
        copy_idx(mesh->get_index_set<1,0>(), v_e0);
        copy_idx(mesh->get_index_set<2,0>(), v_t0);
        copy_idx(mesh->get_index_set<2,1>(), e_t0);
        copy_idx(mesh->get_index_set<3,0>(), v_s0);
        copy_idx(mesh->get_index_set<3,1>(), e_s0);
        copy_idx(mesh->get_index_set<3,2>(), t_s0);

        // okay
        return mesh;
      } // create_big_tetra_mesh_3d

      void validate_refined_big_tetra_mesh_3d(const TetraMesh& mesh)
      {

        // validate sizes
        if(mesh.get_num_entities(0) != 18)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 65)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 84)
          throw String("Triangle count mismatch");
        if(mesh.get_num_entities(3) != 36)
          throw String("Triangle count mismatch");

        // vertex coordinates array
        static const Real vtx0[] =
        {
          -2.0, -4.0,  0.0,
          -2.0,  4.0,  0.0,
           4.0,  0.0,  0.0,
           0.0,  0.0,  42.0,
           0.0,  0.0,  0.0,

           0.0,  0.0,  21.0,
          -1.0, -2.0,  0.0,
          -1.0,  2.0,  0.0,
           2.0,  0.0,  0.0,
          -1.0, -2.0,  21.0,
          -1.0,  2.0,  21.0,
           2.0,  0.0,  21.0,
          -2.0,  0.0,  0.0,
           1.0,  2.0,  0.0,
           1.0, -2.0,  0.0,

          -1.0,  0.0,  10.5,
           0.5,  1.0,  10.5,
           0.5, -1.0,  10.5
        };

        // vertices-at-edge array
        static const Index v_e0[] =
        {
          //edges
          3, 5, //0
          5, 4,
          0, 6,
          6, 4,
          1, 7,
          7, 4, //5
          2, 8,
          8, 4,
          3, 9,
          9, 0,
          3, 10,//10
          10, 1,
          3, 11,
          11, 2,
          0, 12,
          12, 1,//15
          1, 13,
          13, 2,
          2, 14,
          14, 0,

          //triangles
          7, 10,//20
          5, 7,
          10, 5,

          13, 10,
          11, 13,
          10, 11,//25

          12, 10,
          9, 12,
          10, 9,

          9, 11,
          14, 9,//30
          11, 14,

          13, 7,
          8, 13,
          7, 8,

          14, 8,//35
          6, 14,
          8, 6,

          12, 7,
          6, 12,
          7, 6,//40

          11, 8,
          5, 11,
          8, 5,

          6, 9,
          5, 6,//45
          9, 5,

          //tetrahedron
          12, 15,
          9, 15,
          6, 15,
          10, 15,//50
          7, 15,
          5, 15,

          11, 16,
          8, 16,
          13, 16,//55
          5, 16,
          10, 16,
          7, 16,

          11, 17,
          14, 17,//60
          8, 17,
          9, 17,
          5, 17,
          6, 17 //64
        };

        // vertices-at-triangle
        static const Index v_t0[] =
        {
          //triangles
          1, 7, 10, //0
          7, 4, 5,
          10, 5, 3,
          5, 10, 7,

          1, 13, 10,
          13, 2, 11, //5
          10, 11, 3,
          11, 10, 13,

          1, 12, 10,
          12, 0, 9,
          10, 9, 3,//10
          9, 10, 12,

          3, 9, 11,
          9, 0, 14,
          11, 14, 2,
          14, 11, 9,//15

          1, 13, 7,
          13, 2, 8,
          7, 8, 4,
          8, 7, 13,

          2, 14, 8,//20
          14, 0, 6,
          8, 6, 4,
          6, 8, 14,

          1, 12, 7,
          12, 0, 6,//25
          7, 6, 4,
          6, 7, 12,

          2, 11, 8,
          11, 3, 5,
          8, 5, 4,//30
          5, 8, 11,

          0, 6, 9,
          6, 4, 5,
          9, 5, 3,
          5, 9, 6,//35

          //tetrahedron
          12, 9, 6,
          12, 10, 7,
          9, 10, 5,
          6, 7, 5,

          12, 9, 15,//40
          12, 6, 15,
          12, 10, 15,
          12, 7, 15,
          9, 6, 15,
          9, 10, 15,//45
          9, 5, 15,
          6, 7, 15,
          6, 5, 15,
          10, 7, 15,
          10, 5, 15,//50
          7, 5, 15,

          11, 8, 13,
          11, 5, 10,
          8, 5, 7,
          13, 10, 7,//55

          11, 8, 16,
          11, 13, 16,
          11, 5, 16,
          11, 10, 16,
          8, 13, 16,//60
          8, 5, 16,
          8, 7, 16,
          13, 10, 16,
          13, 7, 16,
          5, 10, 16,//65
          5, 7, 16,
          10, 7, 16,

          11, 14, 8,
          11, 9, 5,
          14, 9, 6,//70
          8, 5, 6,

          11, 14, 17,
          11, 8, 17,
          11, 9, 17,
          11, 5, 17,//75
          14, 8, 17,
          14, 9, 17,
          14, 6, 17,
          8, 5, 17,
          8, 6, 17,//80
          9, 5, 17,
          9, 6, 17,
          5, 6, 17
        };

        // vertices-at-tetrahedron
        static const Index v_s0[] =
        {
          12, 6, 9, 0,//0
          12, 9, 6, 15,
          12, 10, 7, 1,
          12, 7, 10, 15,
          9, 5, 10, 3,
          9, 10, 5, 15,//5
          6, 7, 5, 4,
          6, 5, 7, 15,

          10, 7, 5, 15,
          9, 5, 6, 15,
          12, 6, 7, 15,//10
          12, 10, 9, 15,

          11, 13, 8, 2,
          11, 8, 13, 16,
          11, 5, 10, 3,
          11, 10, 5, 16,//15
          8, 7, 5, 4,
          8, 5, 7, 16,
          13, 10, 7, 1,
          13, 7, 10, 16,

          5, 10, 7, 16,//20
          8, 7, 13, 16,
          11, 13, 10, 16,
          11, 5, 8, 16,

          11, 8, 14, 2,
          11, 14, 8, 17,//25
          11, 9, 5, 3,
          11, 5, 9, 17,
          14, 6, 9, 0,
          14, 9, 6, 17,
          8, 5, 6, 4,//30
          8, 6, 5, 17,

          9, 5, 6, 17,
          14, 6, 8, 17,
          11, 8, 5, 17,
          11, 9, 14, 17//35
        };

        // edges-at-triangle
        static const Index e_t0[] =
        {
          // surfaces
          20, 11, 4,
          1, 21, 5,
          0, 10, 22,
          20, 21, 22,

          23, 11, 16,
          13, 24, 17,
          12, 10, 25,
          23, 24, 25,

          26, 11, 15,
          9, 27, 14,
          8, 10, 28,
          26, 27, 28,

          29, 12, 8,
          19, 30, 9,
          18, 13, 31,
          29, 30, 31,

          32, 4, 16,
          6, 33, 17,
          7, 5, 34,
          32, 33, 34,

          35, 6, 18,
          2, 36, 19,
          3, 7, 37,
          35, 36, 37,

          38, 4, 15,
          2, 39, 14,
          3, 5, 40,
          38, 39, 40,

          41, 6, 13,
          0, 42, 12,
          1, 7, 43,
          41, 42, 43,

          44, 9, 2,
          1, 45, 3,
          0, 8, 46,
          44, 45, 46, //108

          // inside
          44, 39, 27,
          20, 38, 26,
          22, 46, 28,
          21, 45, 40,

          48, 47, 27,
          49, 47, 39,
          50, 47, 26,
          51, 47, 38,
          49, 48, 44,
          50, 48, 28,
          52, 48, 46,
          51, 49, 40,
          52, 49, 45,
          51, 50, 20,
          52, 50, 22,
          52, 51, 21,

          33, 24, 41,
          22, 25, 42,
          21, 34, 43,
          20, 32, 23,

          54, 53, 41,
          55, 53, 24,
          56, 53, 42,
          57, 53, 25,
          55, 54, 33,
          56, 54, 43,
          58, 54, 34,
          57, 55, 23,
          58, 55, 32,
          57, 56, 22,
          58, 56, 21,
          58, 57, 20,

          35, 41, 31,
          46, 42, 29,
          44, 36, 30,
          45, 37, 43,

          60, 59, 31,
          61, 59, 41,
          62, 59, 29,
          63, 59, 42,
          61, 60, 35,
          62, 60, 30,
          64, 60, 36,
          63, 61, 43,
          64, 61, 37,
          63, 62, 46,
          64, 62, 44,
          64, 63, 45
        };

        // edges-at-tetrahedron
        static const Index e_s0[] =
        {
          //corners
          39, 27, 14, 44, 2, 9,
          27, 39, 47, 44, 48, 49,
          26, 38, 15, 20, 11, 4,
          38, 26, 47, 20, 51, 50,
          46, 28, 8, 22, 0, 10,
          28, 46, 48, 22, 50, 52,
          40, 45, 3, 21, 5, 1,
          45, 40, 49, 21, 52, 51,

          20, 22, 50, 21, 51, 52,
          46, 44, 48, 45, 52, 49,
          39, 38, 47, 40, 49, 51,
          26, 27, 47, 28, 50, 48,

          24, 41, 13, 33, 17, 6,
          41, 24, 53, 33, 54, 55,
          42, 25, 12, 22, 0, 10,
          25, 42, 53, 22, 57, 56,
          34, 43, 7, 21, 5, 1,
          43, 34, 54, 21, 56, 58,
          23, 32, 16, 20, 11, 4,
          32, 23, 55, 20, 58, 57,


          22, 21, 56, 20, 57, 58,
          34, 33, 54, 32, 58, 55,
          24, 25, 53, 23, 55, 57,
          42, 41, 53, 43, 56, 54,

          41, 31, 13, 35, 6, 18,
          31, 41, 59, 35, 60, 61,
          29, 42, 12, 46, 8, 0,
          42, 29, 59, 46, 63, 62,
          36, 30, 19, 44, 2, 9,
          30, 36, 60, 44, 62, 64,
          43, 37, 7, 45, 1, 3,
          37, 43, 61, 45, 64, 63,

          46, 44, 62, 45, 63, 64,
          36, 35, 60, 37, 64, 61,
          41, 42, 59, 43, 61, 63,
          29, 31, 59, 30, 62, 60
        };

        // triangle-at-tetrahedron
        static const Index t_s0[] =
        {
          32, 9, 25, 36,//0
          44, 41, 40, 36,
          0, 24, 8, 37,
          49, 42, 43, 37,
          2, 10, 34, 38,
          50, 46, 45, 38,//5
          1, 33, 26, 39,
          51, 47, 48, 39,

          51, 50, 49, 3,
          48, 44, 46, 35,
          47, 43, 41, 27,//10
          45, 40, 42, 11,

          17, 28, 5, 52,
          60, 57, 56, 52,
          2, 6, 29, 53,
          65, 58, 59, 53,//15
          1, 30, 18, 54,
          66, 62, 61, 54,
          0, 16, 4, 55,
          67, 63, 64, 55,

          67, 66, 65, 3,//20
          64, 60, 62, 19,
          63, 59, 57, 7,
          61, 56, 58, 31,

          20, 14, 28, 68,
          76, 73, 72, 68,//25
          34, 29, 12, 69,
          81, 74, 75, 69,
          32, 13, 21, 70,
          82, 78, 77, 70,
          33, 22, 30, 71,//30
          83, 79, 80, 71,

          83, 82, 81, 35,
          80, 76, 78, 23,
          79, 75, 73, 31,
          77, 72, 74, 15//35
        };

        // check vertex coordinates array
        if(!comp_vtx(mesh.get_vertex_set(), vtx0))
          throw String("Vertex coordinate refinement failure");

        // check vertices-at-edge array
        if(!comp_idx(mesh.get_index_set<1,0>(), v_e0))
          throw String("Vertex-At-Edge index set refinement failure");

        // check vertices-at-triangle
        if(!comp_idx(mesh.get_index_set<2,0>(), v_t0))
          throw String("Vertex-At-Triangle index set refinement failure");

        // check vertices-at-tetrahedron
        if(!comp_idx(mesh.get_index_set<3,0>(), v_s0))
          throw String("Vertex-At-Tetrahedron index set refinement failure");

        // check edges-at-triangle
        if(!comp_idx(mesh.get_index_set<2,1>(), e_t0))
          throw String("Edge-At-Triangle index set refinement failure");

        // check edges-at-tetrahedron
        if(!comp_idx(mesh.get_index_set<3,1>(), e_s0))
          throw String("Edge-At-Tetrahedron index set refinement failure");

        // check triangles-at-tetrahedron
        if(!comp_idx(mesh.get_index_set<3,2>(), t_s0))
          throw String("Triangle-At-Tetrahedron index set refinement failure");

      } // validate_refined_big_tetra_mesh_3d

      TetraMesh* create_really_big_tetra_mesh_3d()
      {

        Index num_entities[] =
        {
          18, // vertices
          65, // edges
          84, // triangles
          36  // tetrahedron
        };

        // create mesh
        TetraMesh* mesh = new TetraMesh(num_entities);

        // set up vertex coordinates array
        static const Real vtx0[3*18] =
        {
          -2.0, -4.0,  0.0,
          -2.0,  4.0,  0.0,
           4.0,  0.0,  0.0,
           0.0,  0.0,  42.0,
           0.0,  0.0,  0.0,

           0.0,  0.0,  21.0,
          -1.0, -2.0,  0.0,
          -1.0,  2.0,  0.0,
           2.0,  0.0,  0.0,
          -1.0, -2.0,  21.0,
          -1.0,  2.0,  21.0,
           2.0,  0.0,  21.0,
          -2.0,  0.0,  0.0,
           1.0,  2.0,  0.0,
           1.0, -2.0,  0.0,

          -1.0,  0.0,  10.5,
           0.5,  1.0,  10.5,
           0.5, -1.0,  10.5
        };

        // set up vertices-at-edge array
        static const Index v_e0[65*2] =
        {
          //edges
          3, 5, //0
          5, 4,
          0, 6,
          6, 4,
          1, 7,
          7, 4, //5
          2, 8,
          8, 4,
          3, 9,
          9, 0,
          3, 10,//10
          10, 1,
          3, 11,
          11, 2,
          0, 12,
          12, 1,//15
          1, 13,
          13, 2,
          2, 14,
          14, 0,

          //triangles
          7, 10,//20
          5, 7,
          10, 5,

          13, 10,
          11, 13,
          10, 11,//25

          12, 10,
          9, 12,
          10, 9,

          9, 11,
          14, 9,//30
          11, 14,

          13, 7,
          8, 13,
          7, 8,

          14, 8,//35
          6, 14,
          8, 6,

          12, 7,
          6, 12,
          7, 6,//40

          11, 8,
          5, 11,
          8, 5,

          6, 9,
          5, 6,//45
          9, 5,

          //tetrahedron
          12, 15,
          9, 15,
          6, 15,
          10, 15,//50
          7, 15,
          5, 15,

          11, 16,
          8, 16,
          13, 16,//55
          5, 16,
          10, 16,
          7, 16,

          11, 17,
          14, 17,//60
          8, 17,
          9, 17,
          5, 17,
          6, 17 //64
        };

        // set up vertices-at-triangle array
        static const Index v_t0[84*3] =
        {
          //triangles
          1, 7, 10, //0
          7, 4, 5,
          10, 5, 3,
          5, 10, 7,

          1, 13, 10,
          13, 2, 11, //5
          10, 11, 3,
          11, 10, 13,

          1, 12, 10,
          12, 0, 9,
          10, 9, 3,//10
          9, 10, 12,

          3, 9, 11,
          9, 0, 14,
          11, 14, 2,
          14, 11, 9,//15

          1, 13, 7,
          13, 2, 8,
          7, 8, 4,
          8, 7, 13,

          2, 14, 8,//20
          14, 0, 6,
          8, 6, 4,
          6, 8, 14,

          1, 12, 7,
          12, 0, 6,//25
          7, 6, 4,
          6, 7, 12,

          2, 11, 8,
          11, 3, 5,
          8, 5, 4,//30
          5, 8, 11,

          0, 6, 9,
          6, 4, 5,
          9, 5, 3,
          5, 9, 6,//35

          //tetrahedron
          12, 9, 6,
          12, 10, 7,
          9, 10, 5,
          6, 7, 5,

          12, 9, 15,//40
          12, 6, 15,
          12, 10, 15,
          12, 7, 15,
          9, 6, 15,
          9, 10, 15,//45
          9, 5, 15,
          6, 7, 15,
          6, 5, 15,
          10, 7, 15,
          10, 5, 15,//50
          7, 5, 15,

          11, 8, 13,
          11, 5, 10,
          8, 5, 7,
          13, 10, 7,//55

          11, 8, 16,
          11, 13, 16,
          11, 5, 16,
          11, 10, 16,
          8, 13, 16,//60
          8, 5, 16,
          8, 7, 16,
          13, 10, 16,
          13, 7, 16,
          5, 10, 16,//65
          5, 7, 16,
          10, 7, 16,

          11, 14, 8,
          11, 9, 5,
          14, 9, 6,//70
          8, 5, 6,

          11, 14, 17,
          11, 8, 17,
          11, 9, 17,
          11, 5, 17,//75
          14, 8, 17,
          14, 9, 17,
          14, 6, 17,
          8, 5, 17,
          8, 6, 17,//80
          9, 5, 17,
          9, 6, 17,
          5, 6, 17
        };

        // set up vertices-at-tetrahedron array
        static const Index v_s0[36*4] =
        {
          12, 6, 9, 0,//0
          12, 9, 6, 15,
          12, 10, 7, 1,
          12, 7, 10, 15,
          9, 5, 10, 3,
          9, 10, 5, 15,//5
          6, 7, 5, 4,
          6, 5, 7, 15,

          10, 7, 5, 15,
          9, 5, 6, 15,
          12, 6, 7, 15,//10
          12, 10, 9, 15,

          11, 13, 8, 2,
          11, 8, 13, 16,
          11, 5, 10, 3,
          11, 10, 5, 16,//15
          8, 7, 5, 4,
          8, 5, 7, 16,
          13, 10, 7, 1,
          13, 7, 10, 16,

          5, 10, 7, 16,//20
          8, 7, 13, 16,
          11, 13, 10, 16,
          11, 5, 8, 16,

          11, 8, 14, 2,
          11, 14, 8, 17,//25
          11, 9, 5, 3,
          11, 5, 9, 17,
          14, 6, 9, 0,
          14, 9, 6, 17,
          8, 5, 6, 4,//30
          8, 6, 5, 17,

          9, 5, 6, 17,
          14, 6, 8, 17,
          11, 8, 5, 17,
          11, 9, 14, 17//35
        };

        copy_vtx(mesh->get_vertex_set(), vtx0);
        copy_idx(mesh->get_index_set<1,0>(), v_e0);
        copy_idx(mesh->get_index_set<2,0>(), v_t0);
        copy_idx(mesh->get_index_set<3,0>(), v_s0);

        // okay
        return mesh;
      } // create_really_big_tetra_mesh_3d

      TetraSubMesh* create_tria_submesh_3d()
      {

        Index num_entities[] =
        {
          4, // vertices
          6, // edges
          3, // tria
          0 // tetra
        };

        // create mesh
        TetraSubMesh* mesh = new TetraSubMesh(num_entities, true);
        // create a AttributeSet that holds one value for each vertex
        std::unique_ptr<TetraSubMesh::AttributeSetType> my_attrib_set(new TetraSubMesh::AttributeSetType(num_entities[0], 2));
        // Add the attribute to mesh
        mesh->add_attribute(std::move(my_attrib_set), "TriaSubAttributeSet");

        // set up vertex coordinates array
        Real attr[] =
        {
          2.0, 2.0,
          0.0, 0.0,
          2.0, 4.0,
          4.0, 0.0
        };
        copy_attr(*(mesh->find_attribute("TriaSubAttributeSet")), attr);

        // set up vertices-at-edge array
        Index v_e[] =
        {
          3, 1,
          3, 0,
          2, 3,
          0, 2,
          2, 1,
          1, 0
        };
        copy_idx(mesh->get_index_set<1,0>(), v_e);

        // set up vertices-at-tria array
        Index v_t[] =
        {
          0, 3, 1,
          2, 1, 0,
          2, 0, 3
        };
        copy_idx(mesh->get_index_set<2,0>(), v_t);

        // set up edges-at-tria array
        Index e_t[] =
        {
          0, 5, 1,
          5, 3, 4,
          1, 2, 3
        };
        copy_idx(mesh->get_index_set<2,1>(), e_t);

        // set up vertex-target indices
        Index vti[] =
        {
          4, 0, 3, 2
        };
        copy_trg(mesh->get_target_set<0>(), vti);

        // set up edge-target indices
        Index eqi[] =
        {
          9, 3, 6, 0, 4, 1
        };
        copy_trg(mesh->get_target_set<1>(), eqi);

        // set up tria-target indices
        Index tti[] =
        {
          5, 8, 7
        };
        copy_trg(mesh->get_target_set<2>(), tti);
        // okay
        return mesh;
      } // create_tria_submesh_3d()

      void validate_refined_tria_submesh_3d(const TetraSubMesh& mesh)
      {

        // validate sizes
        if(mesh.get_num_entities(0) != 10)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 21)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 12)
          throw String("Triangle count mismatch");

        // check vertex coordinates array
        Real attr[] =
        {
          2.0, 2.0,
          0.0, 0.0,
          2.0, 4.0,
          4.0, 0.0,

          2.0, 0.0,
          3.0, 1.0,
          3.0, 2.0,
          2.0, 3.0,
          1.0, 2.0,
          1.0, 1.0
        };
        if(!comp_attr(*(mesh.find_attribute("TriaSubAttributeSet")), attr))
          throw String("Attribute refinement failure");

        // check vertices-at-edge array
        Index v_e[] =
        {
          3, 4,
          4, 1,
          3, 5,
          5, 0,
          2, 6,
          6, 3,
          0, 7,
          7, 2,
          2, 8,
          8, 1,
          1, 9,
          9, 0,
          5, 9,
          4, 5,
          9, 4,
          8, 7,
          9, 8,
          7, 9,
          7, 6,
          5, 7,
          6, 5
        };
        if(!comp_idx(mesh.get_index_set<1,0>(), v_e))
          throw String("Vertex-At-Edge index set refinement failure");

        // check vertices-at-tria
        Index v_t[] =
        {
          0, 5, 9,
          5, 3, 4,
          9, 4, 1,
          4, 9, 5,
          2, 8, 7,
          8, 1, 9,
          7, 9, 0,
          9, 7, 8,
          2, 7, 6,
          7, 0, 5,
          6, 5, 3,
          5, 6, 7
        };
        if(!comp_idx(mesh.get_index_set<2,0>(), v_t))
          throw String("Vertex-At-Triangle index set refinement failure");

        // check edges-at-tria
        Index e_t[] =
        {
          12, 11, 3,
          0, 13, 2,
          1, 10, 14,
          12, 13, 14,
          15, 7, 8,
          10, 16, 9,
          11, 6, 17,
          15, 16, 17,
          18, 4, 7,
          3, 19, 6,
          2, 5, 20,
          18, 19, 20
        };
        if(!comp_idx(mesh.get_index_set<2,1>(), e_t))
          throw String("Edges-At-Triangle refinement failure");

        // check vertex-target indices
        Index vti[] =
        {
          4, 0, 3, 2, 14, 8, 11, 5, 9, 6
        };
        if(!comp_trg(mesh.get_target_set<0>(), vti))
          throw String("Vertex-Target-Indices refinement failure");

        // check edge-target indices
        Index eti[] =
        {
          18, 19, 6, 7, 12, 13, 1, 0, 8, 9, 2, 3,
          37, 35, 36, 46, 44, 45, 42, 43, 41
        };
        if(!comp_trg(mesh.get_target_set<1>(), eti))
          throw String("Edge-Target-Indices refinement failure");

        // check tria-target indices
        Index tti[] =
        {
          22, 20, 21, 23,
          34, 32, 33, 35,
          29, 30, 28, 31
        };
        if(!comp_trg(mesh.get_target_set<2>(), tti))
          throw String("Triangle-Target-Indices refinement failure");
      } //validate_refined_tria_submesh_3d

    } // namespace TestAux
  } // namespace Geometry
} // namespace FEAT

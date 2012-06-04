#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_SIMPLEX_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_SIMPLEX_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {

      typedef ConformalMesh< ConformalMeshPolicy< Shape::Tetrahedron > > TetrahedronMesh;

      TetrahedronMesh* create_tetrahedronrefinement_mesh_3d(int orientation)
      {

        Index num_entities[] =
        {
          4, // vertices
          6, // edges
          4, // triangles
          1  // tetrahedron
        };

        // create mesh
        TetrahedronMesh* mesh = new TetrahedronMesh(num_entities);

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

        //
        //                          v_3
        //                          /\\
        //                         /  \  \
        //                        /    \    \
        //                       /      \      \
        //                      /        \        \
        //                     /          \          \
        //                    /            \           v e_5
        //                   /              \              \
        //                  /                \                \
        //                 /                  \                  \
        //                /                    \                    \
        //               /                  e_3 \                      \ v_2
        //          e_4 v                        v                 /    /
        //             /                          \           /        /
        //            /                            \      /           /
        //           /                              \ /              /
        //          /                       e_0  /   \              /
        //         /              s_0        ^        \            /
        //        /                    /               \          ^ e_1
        //       /                /                     \        /
        //      /            /                           \      /
        //     /        /                                 \    /
        //    /    /                                       \  /
        //   //____________________<________________________\/
        //  v_0                    e_2                      v_1
        //
        //
        // orientation:
        //    t0: v_2 - v_1 - v_0
        //    t1: v_2 - v_3 - v_1
        //    t2: v_3 - v_0 - v_2
        //    t3: v_0 - v_3 - v_1
        // tetrahedron:
        //        v_0 - v_1 - v_2 - v_3


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

        //
        //                          v_2
        //                          /\\
        //                         /  \  \
        //                        /    \    \
        //                       /      \      \
        //                      /        \        \
        //                     /          \          \
        //                    /            \           ^ e_5
        //                   /              \              \
        //                  /                \                \
        //                 /                  \                  \
        //                /                    \                    \
        //               /                  e_3 \                      \ v_0
        //          e_0 v                        v                 /    /
        //             /                          \           /        /
        //            /                            \      /           /
        //           /                              \ /              /
        //          /                       e_4  /   \              /
        //         /              s_0        ^        \            /
        //        /                    /               \          ^ e_2
        //       /                /                     \        /
        //      /            /                           \      /
        //     /        /                                 \    /
        //    /    /                                       \  /
        //   //____________________>________________________\/
        //  v_3                    e_1                      v_1
        //
        //
        // orientation:
        //    t0: v_0 - v_3 - v_2
        //    t1: v_0 - v_1 - v_3
        //    t2: v_0 - v_1 - v_2
        //    t3: v_2 - v_1 - v_3
        // tetrahedron:
        //        v_3 - v_1 - v_0 - v_2

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
        }
        // okay
        return mesh;
      } // create_tetrahedronrefinement_mesh_3d

      void validate_refined_tetrahedronrefinement_mesh_3d(const TetrahedronMesh& mesh, int orientation)
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
          4, 5, 6, 0,
          4, 5, 6, 10,
          4, 7, 8, 1,
          4, 7, 8, 10,
          5, 7, 9, 2,
          5, 7, 9, 10,
          6, 8, 9, 3,
          6, 8, 9, 10,
          9, 8, 7, 10,
          9, 6, 5, 10,
          8, 6, 4, 10,
          7, 5, 4, 10
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
          21, 18, 0, 15, 2, 4,
          21, 18, 24, 15, 25, 26,

          22, 19, 1, 12, 6, 8,
          22, 19, 24, 12, 27, 28,

          23, 16, 3, 13, 7, 10,
          23, 16, 25, 13, 27, 29,

          20, 17, 5, 14, 9, 11,
          20, 17, 26, 14, 28, 29,

          // surfaces
          14, 13, 29, 12, 28, 27,
          17, 16, 29, 15, 26, 25,
          20, 19, 28, 18, 26, 24,
          23, 22, 27, 21, 25, 24
        };

        // triangle-at-tetrahedron
        static const Index t_s0[] =
        {
          4, 8, 12, 16,
          24, 21, 20, 16,

          0, 9, 13, 17,
          29, 23, 22, 17,

          1, 5, 14, 18,
          30, 26, 25, 18,

          2, 6, 10, 19,
          31, 28, 27, 19,

          29, 30, 31, 3,
          24, 26, 28, 7,
          21, 23, 27, 11,
          20, 22, 25, 15
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
          6, 4, 8, 0,
          6, 4, 8, 10,
          6, 5, 7, 1,
          6, 5, 7, 10,
          4, 5, 9, 2,
          4, 5, 9, 10,
          8, 7, 9, 3,
          8, 7, 9, 10,

          9, 7, 5, 10,
          9, 8, 4, 10,
          7, 8, 6, 10,
          5, 4, 6, 10
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
          14, 21, 5, 19, 0, 9,
          14, 21, 24, 19, 25, 26,
          13, 23, 4, 17, 2, 7,
          13, 23, 24, 17, 27, 28,
          12, 20, 1, 15, 3, 11,
          12, 20, 25, 15, 27, 29,
          22, 18, 8, 16, 6, 10,
          22, 18, 26, 16, 28, 29,

          // surfaces
          16, 15, 29, 17, 28, 27,
          18, 20, 29, 19, 26, 25,
          22, 23, 28, 21, 26, 24,
          12, 13, 27, 14, 25, 24
        };

        // triangle-at-tetrahedron
        static const Index t_s1[] =
        {
          9, 12, 2, 16,
          24, 21, 20, 16,

          6, 14, 1, 17,
          29, 23, 22, 17,

          4, 10, 0, 18,
          30, 26, 25, 18,

          5, 8, 13, 19,
          31, 28, 27, 19,

          29, 30, 31, 7,
          24, 26, 28, 11,
          21, 23, 27, 15,
          20, 22, 25, 3
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
          5, 8, 4, 3,
          5, 8, 4, 10,
          5, 6, 7, 1,
          5, 6, 7, 10,
          8, 6, 9, 0,
          8, 6, 9, 10,
          4, 7, 9, 2,
          4, 7, 9, 10,

          9, 7, 6, 10,
          9, 4, 8, 10,
          7, 4, 5, 10,
          6, 8, 5, 10
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
          17, 23, 2,  13, 8, 1,
          17, 23, 24, 13, 25, 26,
          16, 22, 3, 19, 4, 7,
          16, 22, 24, 19, 27, 28,
          15, 12, 9, 18, 5, 10,
          15, 12, 25, 18, 27, 29,
          21, 14, 0, 20, 6, 11,
          21, 14, 26, 20, 28, 29,

          // surfaces
          20, 18, 29, 19, 28, 27,
          14, 12, 29, 13, 26, 25,
          21, 22, 28, 23, 26, 24,
          15, 16, 27, 17, 25, 24
        };

        // triangle-at-tetrahedron
        static const Index t_s2[] =
        {
          1, 14, 6, 16,
          24, 21, 20, 16,

          9, 13, 5, 17,
          29, 23, 22, 17,

          8, 0, 4, 18,
          30, 26, 25, 18,

          10, 2, 12, 19,
          31, 28, 27, 19,

          29, 30, 31, 11,

          24, 26, 28, 3,

          21, 23, 27, 15,

          20, 22, 25, 7
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

        } //switch
      } // validate_refined_tetrahedronrefinement_mesh_3d

      TetrahedronMesh* create_bigtetrahedronrefinement_mesh_3d()
      {

        Index num_entities[] =
        {
          5, // vertices
          10, // edges
          9, // triangles
          3  // tetrahedron
        };

        // create mesh
        TetrahedronMesh* mesh = new TetrahedronMesh(num_entities);

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
      } // create_bigtetrahedronrefinement_mesh_3d

      void validate_refined_bigtetrahedronrefinement_mesh_3d(const TetrahedronMesh& mesh)
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
          12, 9, 6, 0,//0
          12, 9, 6, 15,
          12, 10, 7, 1,
          12, 10, 7, 15,
          9, 10, 5, 3,
          9, 10, 5, 15,//5
          6, 7, 5, 4,
          6, 7, 5, 15,
          5, 7, 10, 15,
          5, 6, 9, 15,
          7, 6, 12, 15,//10
          10, 9, 12, 15,

          11, 8, 13, 2,
          11, 8, 13, 16,
          11, 5, 10, 3,
          11, 5, 10, 16,//15
          8, 5, 7, 4,
          8, 5, 7, 16,
          13, 10, 7, 1,
          13, 10, 7, 16,
          7, 10, 5, 16,//20
          7, 13, 8, 16,
          10, 13, 11, 16,
          5, 8, 11, 16,

          11, 14, 8, 2,
          11, 14, 8, 17,//25
          11, 9, 5, 3,
          11, 9, 5, 17,
          14, 9, 6, 0,
          14, 9, 6, 17,
          8, 5, 6, 4,//30
          8, 5, 6, 17,
          6, 5, 9, 17,
          6, 8, 14, 17,
          5, 8, 11, 17,
          9, 14, 11, 17//35
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
          27, 39, 14, 44, 9, 2,
          27, 39, 47, 44, 48, 49,
          26, 38, 15, 20, 11, 4,
          26, 38, 47, 20, 50, 51,
          28, 46, 8, 22, 10, 0,
          28, 46, 48, 22, 50, 52,
          40, 45, 3, 21, 5, 1,
          40, 45, 49, 21, 51, 52,
          21, 22, 52, 20, 51, 50,
          45, 46, 52, 44, 49, 48,
          40, 38, 51, 39, 49, 47,
          28, 26, 50, 27, 48, 47,

          41, 24, 13, 33, 6, 17,
          41, 24, 53, 33, 54, 55,
          42, 25, 12, 22, 0, 10,
          42, 25, 53, 22, 56, 57,
          43, 34, 7, 21, 1, 5,
          43, 34, 54, 21, 56, 58,
          23, 32, 16, 20, 11, 4,
          23, 32, 55, 20, 57, 58,
          20, 21, 58, 22, 57, 56,
          32, 34, 58, 33, 55, 54,
          23, 25, 57, 24, 55, 53,
          43, 42, 56, 41, 54, 53,

          31, 41, 13, 35, 18, 6,
          31, 41, 59, 35, 60, 61,
          29, 42, 12, 46, 8, 0,
          29, 42, 59, 46, 62, 63,
          30, 36, 19, 44, 9, 2,
          30, 36, 60, 44, 62, 64,
          43, 37, 7, 45, 1, 3,
          43, 37, 61, 45, 63, 64,
          45, 44, 64, 46, 63, 62,
          37, 36, 64, 35, 61, 60,
          43, 42, 63, 41, 61, 59,
          30, 29, 62, 31, 60, 59
        };

        // triangle-at-tetrahedron
        static const Index t_s0[] =
        {
          32, 25, 9, 36,//0
          44, 41, 40, 36,
          0, 24, 8, 37,
          49, 43, 42, 37,
          2, 34, 10, 38,
          50, 46, 45, 38,//5
          1, 33, 26, 39,
          51, 48, 47, 39,
          49, 50, 51, 3,
          44, 46, 48, 35,
          41, 43, 47, 27,//10
          40, 42, 45, 11,

          17, 5, 28, 52,
          60, 57, 56, 52,
          2, 6, 29, 53,
          65, 59, 58, 53,//15
          1, 18, 30, 54,
          66, 62, 61, 54,
          0, 16, 4, 55,
          67, 64, 63, 55,
          65, 66, 67, 3,//20
          60, 62, 64, 19,
          57, 59, 63, 7,
          56, 58, 61, 31,

          20, 28, 14, 68,
          76, 73, 72, 68,//25
          34, 29, 12, 69,
          81, 75, 74, 69,
          32, 21, 13, 70,
          82, 78, 77, 70,
          33, 22, 30, 71,//30
          83, 80, 79, 71,
          81, 82, 83, 35,
          76, 78, 80, 23,
          73, 75, 79, 31,
          72, 74, 77, 15//35
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

      } // validate_refined_bigtetrahedronrefinement_mesh_3d

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_Tetrahedron_HPP
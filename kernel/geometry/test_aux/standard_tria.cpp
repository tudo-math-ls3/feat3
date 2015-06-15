#include <kernel/geometry/test_aux/standard_tria.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace TestAux
    {
      TriaMesh* create_tria_mesh_2d(int orientation)
      {

        Index num_entities[] =
        {
          3, // vertices
          3, // edges
          1, // triangles
        };

        // create mesh
        TriaMesh* mesh = new TriaMesh(num_entities);

        // first possibility (standard)

        /*
           v_2
            |\
            | \
            |  ^
            v   \
            |    \
           e_1   e_0
            |      \
            |       \
            |        ^
            v         \
            |          \
           v_0->-e_2-->-v_1

           triangle orientation: 0-1-2
        */

        // set up vertex coordinates array
        static const Real vtx0[3*2] =
        {
          0.0, 0.0,
          1.0, 0.0,
          0.0, 1.0
        };

        // set up vertices-at-edge array
        static const Index v_e0[3*2] =
        {
          1, 2,
          2, 0,
          0, 1
        };

        // set up vertices-at-triangle array
        static const Index v_t0[3] =
        {
          0, 1, 2
        };

        // set up edges-at-triangle array
        static const Index e_t0[3] =
        {
          0, 1, 2
        };

        // second possibility

        /*
           v_2
            |\
            | \
            |  v
            v   \
            |    \
           e_1   e_0
            |      \
            |       \
            |        v
            v         \
            |          \
           v_0->-e_2-->-v_1

           triangle orientation: 0-1-2
        */

        // set up vertex coordinates array
        static const Real vtx1[3*2] =
        {
          0.0, 0.0,
          1.0, 0.0,
          0.0, 1.0
        };

        // set up vertices-at-edge array
        static const Index v_e1[3*2] =
        {
          2, 1,
          2, 0,
          0, 1
        };

        // set up vertices-at-triangle array
        static const Index v_t1[3*1] =
        {
          0, 1, 2
        };

        // set up edges-at-triangle array
        static const Index e_t1[1*3] =
        {
          0, 1, 2
        };

        // third possibility

        /*
           v_2
            |\
            | \
            |  v
            v   \
            |    \
           e_1   e_2
            |      \
            |       \
            |        v
            v         \
            |          \
           v_1-<-e_0--<-v_0

           triangle orientation: 1-0-2
        */

        // set up vertex coordinates array
        static const Real vtx2[3*2] =
        {
          1.0, 0.0,
          0.0, 0.0,
          0.0, 1.0
        };

        // set up vertices-at-edge array
        static const Index v_e2[3*2] =
        {
          0, 1,
          2, 1,
          2, 0
        };

        // set up vertices-at-triangle array
        static const Index v_t2[3*1] =
        {
          1, 0, 2
        };

        // set up edges-at-triangle array
        static const Index e_t2[3*1] =
        {
          2, 1, 0
        };

        // fourth possibility

        /*
           v_0
            |\
            | \
            |  v
            v   \
            |    \
           e_0   e_2
            |      \
            |       \
            |        v
            v         \
            |          \
           v_2-<-e_1--<-v_1

           triangle orientation: 0-1-2
        */

        // set up vertex coordinates array
        static const Real vtx3[3*2] =
        {
          0.0, 1.0,
          1.0, 0.0,
          0.0, 0.0
        };

        // set up vertices-at-edge array
        static const Index v_e3[3*2] =
        {
          2, 0,
          1, 2,
          0, 1
        };

        // set up vertices-at-triangle array
        static const Index v_t3[3*1] =
        {
          2, 1, 0
        };

        // set up edges-at-triangle array
        static const Index e_t3[3] =
        {
          2, 0, 1
        };

        switch(orientation)
        {
          case 0:
            copy_vtx(mesh->get_vertex_set(), vtx0);
            copy_idx(mesh->get_index_set<1,0>(), v_e0);
            copy_idx(mesh->get_index_set<2,0>(), v_t0);
            copy_idx(mesh->get_index_set<2,1>(), e_t0);
            break;
          case 1:
            copy_vtx(mesh->get_vertex_set(), vtx1);
            copy_idx(mesh->get_index_set<1,0>(), v_e1);
            copy_idx(mesh->get_index_set<2,0>(), v_t1);
            copy_idx(mesh->get_index_set<2,1>(), e_t1);
            break;
          case 2:
            copy_vtx(mesh->get_vertex_set(), vtx2);
            copy_idx(mesh->get_index_set<1,0>(), v_e2);
            copy_idx(mesh->get_index_set<2,0>(), v_t2);
            copy_idx(mesh->get_index_set<2,1>(), e_t2);
            break;
          case 3:
            copy_vtx(mesh->get_vertex_set(), vtx3);
            copy_idx(mesh->get_index_set<1,0>(), v_e3);
            copy_idx(mesh->get_index_set<2,0>(), v_t3);
            copy_idx(mesh->get_index_set<2,1>(), e_t3);
            break;
        }
        // okay
        return mesh;
      } // create_trianglerefinement_mesh_2d

      void validate_refined_tria_mesh_2d(const TriaMesh& mesh, int orientation)
      {

        // validate sizes
        if(mesh.get_num_entities(0) != 6)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 9)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 4)
          throw String("Triangle count mismatch");


        // first possibility (standard)

        // vertex coordinates array
        static const Real vtx0[] =
        {
          0.0, 0.0,
          1.0, 0.0,
          0.0, 1.0,
          0.5, 0.5,
          0.0, 0.5,
          0.5, 0.0
        };

        // vertices-at-edge array
        static const Index v_e0[] =
        {
          1, 3,
          3, 2,
          2, 4,
          4, 0,
          0, 5,
          5, 1,
          5, 4,
          3, 5,
          4, 3
        };

        // vertices-at-triangle
        static const Index v_t0[] =
        {
          0, 5, 4,
          5, 1, 3,
          4, 3, 2,
          3, 4, 5
        };

        // edges-at-triangle
        static const Index e_t0[] =
        {
          6, 3, 4,
          0, 7, 5,
          1, 2, 8,
          6, 7, 8
        };

        // second possibility

        // vertex coordinates array
        static const Real vtx1[] =
        {
          0.0, 0.0,
          1.0, 0.0,
          0.0, 1.0,
          0.5, 0.5,
          0.0, 0.5,
          0.5, 0.0
        };

        // vertices-at-edge array
        static const Index v_e1[] =
        {
          2, 3,
          3, 1,
          2, 4,
          4, 0,
          0, 5,
          5, 1,
          5, 4,
          3, 5,
          4, 3
        };

        // vertices-at-triangle
        static const Index v_t1[] =
        {
          0, 5, 4,
          5, 1, 3,
          4, 3, 2,
          3, 4, 5
        };

        // edges-at-triangle
        static const Index e_t1[] =
        {
          6, 3, 4,
          1, 7, 5,
          0, 2, 8,
          6, 7, 8
        };

        // third possibility

        // vertex coordinates array
        static const Real vtx2[] =
        {
          1.0, 0.0,
          0.0, 0.0,
          0.0, 1.0,
          0.5, 0.0,
          0.0, 0.5,
          0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e2[] =
        {
          0, 3,
          3, 1,
          2, 4,
          4, 1,
          2, 5,
          5, 0,
          3, 4,
          5, 3,
          4, 5
        };

        // vertices-at-triangle
        static const Index v_t2[] =
        {
          1, 3, 4,
          3, 0, 5,
          4, 5, 2,
          5, 4, 3
        };

        // edges-at-triangle
        static const Index e_t2[] =
        {
          6, 3, 1,
          5, 7, 0,
          4, 2, 8,
          6, 7, 8
        };

        // fourth possibility

        // vertex coordinates array
        static const Real vtx3[] =
        {
          0.0, 1.0,
          1.0, 0.0,
          0.0, 0.0,
          0.0, 0.5,
          0.5, 0.0,
          0.5, 0.5
        };

        // vertices-at-edge array
        static const Index v_e3[] =
        {
          2, 3,
          3, 0,
          1, 4,
          4, 2,
          0, 5,
          5, 1,
          4, 3,
          5, 4,
          3, 5
        };

        // vertices-at-triangle
        static const Index v_t3[] =
        {
          2, 4, 3,
          4, 1, 5,
          3, 5, 0,
          5, 3, 4
        };

        // edges-at-triangle
        static const Index e_t3[] =
        {
          6, 0, 3,
          5, 7, 2,
          4, 1, 8,
          6, 7, 8
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

            // check edges-at-triangle
            if(!comp_idx(mesh.get_index_set<2,1>(), e_t0))
              throw String("Edge-At-Triangle index set refinement failure");
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

            // check edges-at-triangle
            if(!comp_idx(mesh.get_index_set<2,1>(), e_t1))
              throw String("Edge-At-Triangle index set refinement failure");
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

            // check edges-at-triangle
            if(!comp_idx(mesh.get_index_set<2,1>(), e_t2))
              throw String("Edge-At-Triangle index set refinement failure");
            break;

          case 3:
            // check vertex coordinates array
            if(!comp_vtx(mesh.get_vertex_set(), vtx3))
              throw String("Vertex coordinate refinement failure");

            // check vertices-at-edge array
            if(!comp_idx(mesh.get_index_set<1,0>(), v_e3))
              throw String("Vertex-At-Edge index set refinement failure");

            // check vertices-at-triangle
            if(!comp_idx(mesh.get_index_set<2,0>(), v_t3))
              throw String("Vertex-At-Triangle index set refinement failure");

            // check edges-at-triangle
            if(!comp_idx(mesh.get_index_set<2,1>(), e_t3))
              throw String("Edge-At-Triangle index set refinement failure");
            break;

        } //switch
      } // validate_refined_trianglerefinement_mesh_2d

      TriaMesh* create_patch_tria_mesh_2d()
      {

        Index num_entities[] =
        {
          5, // vertices
          8, // edges
          4,   // triangles
        };

        // create mesh
        TriaMesh* mesh = new TriaMesh(num_entities);

        //
        //   v_3________<________v_4(10,10)
        //     |\      e_2      /|
        //     | \             / |
        //     |  \   t_2     /  |
        //     |   \         /   |
        //     |    \    e_6/    |
        //     |  e_7v     v     |
        //     |      \   / t_1  |
        //     | t_3   \ /       |
        //  e_3v     v_2(6,7)    ^e_1
        //     |       / \       |
        //     |      /   \      |
        //     |     /     \e_5  |
        //     |    ^e_4    ^    |
        //     |   /         \   |
        //     |  /    t_0    \  |
        //     | /             \ |
        //     |/_______>_______\|
        //  v_0(0,0)   e_0      v_1
        //
        // triangle orientation:
        //     t_0: v_0-v_1-v_2
        //     t_1: v_2-v_4-v_1
        //     t_2: v_3-v_4-v_2
        //     t_3: v_2-v_3-v_0

        // set up vertex coordinates array
        static const Real vtx0[3*8] =
        {
          0.0, 0.0,
          10.0, 0.0,
          6.0, 7.0,
          0.0, 10.0,
          10.0, 10.0
        };

        // set up vertices-at-edge array
        static const Index v_e0[12*2] =
        {
          0, 1,
          1, 4,
          4, 3,
          3, 0,
          0, 2,
          1, 2,
          4, 2,
          3, 2
        };

        // set up vertices-at-triangle array
        static const Index v_t0[6*4] =
        {
          0, 1, 2,
          2, 4, 1,
          3, 4, 2,
          2, 3, 0
        };

        // set up edges-at-triangle array
        static const Index e_t0[6*4] =
        {
          5, 4, 0,
          1, 5, 6,
          6, 7, 2,
          3, 4, 7
        };

        copy_vtx(mesh->get_vertex_set(), vtx0);
        copy_idx(mesh->get_index_set<1,0>(), v_e0);
        copy_idx(mesh->get_index_set<2,0>(), v_t0);
        copy_idx(mesh->get_index_set<2,1>(), e_t0);

        // okay
        return mesh;
      } // create_patch_tria_mesh_2d

      void validate_refined_patch_tria_mesh_2d(const TriaMesh& mesh)
      {

        // validate sizes
        if(mesh.get_num_entities(0) != 13)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 28)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 16)
          throw String("Triangle count mismatch");

        // vertex coordinates array
        static const Real vtx0[] =
        {
          0.0, 0.0,
          10.0, 0.0,
          6.0, 7.0,
          0.0, 10.0,
          10.0, 10.0,
          5.0, 0.0,
          10.0, 5.0,
          5.0, 10.0,
          0.0, 5.0,
          3.0, 3.5,
          8.0, 3.5,
          8.0, 8.5,
          3.0, 8.5
        };

        // vertices-at-edge array
        static const Index v_e0[] =
        {
          0, 5,
          5, 1,
          1, 6,
          6, 4,
          4, 7,
          7, 3,
          3, 8,
          8, 0,
          0, 9,
          9, 2,
          1, 10,
          10, 2,
          4, 11,
          11, 2,
          3, 12,
          12, 2,
          5, 9,
          10, 5,
          9, 10,
          11, 10,
          6, 11,
          10, 6,
          7, 12,
          11, 7,
          12, 11,
          12, 9,
          8, 12,
          9, 8
        };

        // vertices-at-triangle
        static const Index v_t0[] =
        {
          0, 5, 9, //0
          5, 1, 10,
          9, 10, 2,
          10, 9, 5,
          2, 11, 10,
          11, 4, 6,//5
          10, 6, 1,
          6, 10, 11,
          3, 7, 12,
          7, 4, 11,
          12, 11, 2,//10
          11, 12, 7,
          2, 12, 9,
          12, 3, 8,
          9, 8, 0,
          8, 9, 12//15
        };

        // edges-at-triangle
        static const Index e_t0[] =
        {
          16, 8, 0,
          10, 17, 1,
          11, 9, 18,
          16, 17, 18,
          19, 11, 13,
          3, 20, 12,
          2, 10, 21,
          19, 20, 21,
          22, 14, 5,
          12, 23, 4,
          13, 15, 24,
          22, 23, 24,
          25, 9, 15,
          6, 26, 14,
          7, 8, 27,
          25, 26, 27
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

        // check edges-at-triangle
        if(!comp_idx(mesh.get_index_set<2,1>(), e_t0))
          throw String("Edge-At-Triangle index set refinement failure");

      } // validate_refined_patch_tria_mesh_2d

      TriaSubMesh* create_patch_tria_submesh_2d()
      {
        /*
           v_0(0,10)
             |\
             | \
             |  \
             |   \
             |    \
             |  e_3^
             |      \
             | t_0   \
          e_0v     v_3(6,7)
             |       / \
             |      /   \
             |     /     \e_2
             |    ve_1    ^
             |   /         \
             |  /    t_1    \
             | /             \
             |/_______<_______\
          v_1(0,0)   e_4      v_2(10,0)

         triangle orientation:
             t_0: v_0-v_3-v_1
             t_1: v_1-v_3-v_2
        */

        Index num_entities[] =
        {
          4, // vertices
          5, // edges
          2  // triangles
        };

        // create mesh
        TriaSubMesh* mesh = new TriaSubMesh(num_entities, true);
        // create a MeshAttribute that holds one value for each vertex
        TriaSubMesh::AttributeType my_vertex_set(num_entities[0],2);
        // Add the attribute to mesh
        mesh->add_attribute<0>(my_vertex_set);

        // set up vertex coordinates array
        Real vtx[] =
        {
          0.0, 10.0,
          0.0, 0.0,
          10.0, 0.0,
          6.0, 7.0
        };
        copy_vtx(mesh->get_attributes<0>()[0], vtx);

        // set up vertices-at-edge array
        Index v_e[] =
        {
          0, 1,
          3, 1,
          2, 3,
          3, 0,
          2, 1
        };
        copy_idx(mesh->get_index_set<1,0>(), v_e);

        // set up vertices-at-triangle array
        Index v_t[] =
        {
          0, 3, 1,
          1, 3, 2
        };
        copy_idx(mesh->get_index_set<2,0>(), v_t);

        // set up edges-at-triangle array
        Index e_t[] =
        {
          1, 0, 3,
          2, 4, 1
        };
        copy_idx(mesh->get_index_set<2,1>(), e_t);

        // set up vertex-target indices
        Index vti[] =
        {
          3, 0, 1, 2
        };
        copy_trg(mesh->get_target_set<0>(), vti);

        // set up edge-target indices
        Index eti[] =
        {
          3, 4, 5, 7, 0
        };
        copy_trg(mesh->get_target_set<1>(), eti);

        // set up triangle-target indices
        Index tti[] =
        {
          3, 0
        };
        copy_trg(mesh->get_target_set<2>(), tti);
        // okay
        return mesh;
      } // create_patch_tria_submesh_2d

      void validate_refined_patch_tria_submesh_2d(const TriaSubMesh& mesh)
      {

        // validate sizes
        if(mesh.get_num_entities(0) != 9)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 16)
          throw String("Edge count mismatch");
        if(mesh.get_num_entities(2) != 8)
          throw String("Triangle count mismatch");

        // check vertex coordinates array
        Real vtx[] =
        {
          0.0, 10.0,
          0.0, 0.0,
          10.0, 0.0,
          6.0, 7.0,
          0.0, 5.0,
          3.0, 3.5,
          8.0, 3.5,
          3.0, 8.5,
          5.0, 0.0
        };
        if(!comp_vtx(mesh.get_attributes<0>()[0], vtx))
          throw String("Vertex coordinate refinement failure");

        // check vertices-at-edge array
        Index v_e[] =
        {
          0, 4,
          4, 1,
          3, 5,
          5, 1,
          2, 6,
          6, 3,
          3, 7,
          7, 0,
          2, 8,
          8, 1,
          7, 4,
          5, 7,
          4, 5,
          5, 8,
          6, 5,
          8, 6
        };
        if(!comp_idx(mesh.get_index_set<1,0>(), v_e))
          throw String("Vertex-At-Edge index set refinement failure");

        // check vertices-at-triangle
        Index v_t[] =
        {
          0, 7, 4,
          7, 3, 5,
          4, 5, 1,
          5, 4, 7,
          1, 5, 8,
          5, 3, 6,
          8, 6, 2,
          6, 8, 5
        };
        if(!comp_idx(mesh.get_index_set<2,0>(), v_t))
          throw String("Vertex-At-Tria index set refinement failure");

        // check edges-at-triangle
        Index e_t[] =
        {
          10, 0, 7,
          2, 11, 6,
          3, 1, 12,
          10, 11, 12,
          13, 9, 3,
          5, 14, 2,
          4, 8, 15,
          13, 14, 15
        };
        if(!comp_idx(mesh.get_index_set<2,1>(), e_t))
          throw String("Edges-At-Tria refinement failure");

        // check vertex-target indices
        Index vti[] =
        {
          3, 0, 1, 2, 8, 9, 10, 12, 5
        };
        if(!comp_trg(mesh.get_target_set<0>(), vti))
          throw String("Vertex-Target-Indices refinement failure");

        // check edge-target indices
        Index eti[] =
        {
          6, 7, 9, 8, 10, 11, 15, 14, 1, 0, 26, 25, 27, 16, 18, 17
        };
        if(!comp_trg(mesh.get_target_set<1>(), eti))
          throw String("Edge-Target-Indices refinement failure");

        // check triangle-target indices
        Index tti[] =
        {
          13, 12, 14, 15, 0, 2, 1, 3
        };
        if(!comp_trg(mesh.get_target_set<2>(), tti))
          throw String("Tria-Target-Indices refinement failure");
      } // validate_refined_patch_tria_submesh_2d

      TriaSubMesh* create_patch_edge_submesh_2d()
      {

        Index num_entities[] =
        {
          5, // vertices
          4, // edges
          0  // triangles
        };

        // create mesh
        TriaSubMesh* mesh = new TriaSubMesh(num_entities, true);
        // create a MeshAttribute that holds one value for each vertex
        TriaSubMesh::AttributeType my_vertex_set(num_entities[0], 1, 1);
        // Add the attribute to mesh
        mesh->add_attribute<0>(my_vertex_set);

        // set up vertex coordinates array
        Real vtx[] =
        {
          0.0,
          1.0,
          2.0,
          3.0,
          4.0
        };
        copy_vtx(mesh->get_attributes<0>()[0], vtx);

        // set up vertices-at-edge array
        Index v_e[] =
        {
          0, 1,
          1, 2,
          2, 3,
          3, 4,
        };
        copy_idx(mesh->get_index_set<1,0>(), v_e);

        // set up vertex-target indices
        Index vti[] =
        {
          1, 4, 2, 3, 0
        };
        copy_trg(mesh->get_target_set<0>(), vti);

        // set up edge-target indices
        Index eti[] =
        {
          1, 6, 7, 3
        };
        copy_trg(mesh->get_target_set<1>(), eti);

        // okay
        return mesh;
      } // create_patch_edge_submesh_2d

      void validate_refined_patch_edge_submesh_2d(const TriaSubMesh& mesh)
      {

        // validate sizes
        if(mesh.get_num_entities(0) != 9)
          throw String("Vertex count mismatch");
        if(mesh.get_num_entities(1) != 8)
          throw String("Edge count mismatch");

        // check vertex coordinates array
        Real vtx[] =
        {
          0.0,
          1.0,
          2.0,
          3.0,
          4.0,
          0.5,
          1.5,
          2.5,
          3.5
        };
        if(!comp_vtx(mesh.get_attributes<0>()[0], vtx))
          throw String("Vertex coordinate refinement failure");

        // check vertices-at-edge array
        Index v_e[] =
        {
          0, 5,
          5, 1,
          1, 6,
          6, 2,
          2, 7,
          7, 3,
          3, 8,
          8, 4
        };
        if(!comp_idx(mesh.get_index_set<1,0>(), v_e))
          throw String("Vertex-At-Edge index set refinement failure");

        // check vertex-target indices
        Index vti[] =
        {
          1, 4, 2, 3, 0, 6, 11, 12, 8
        };
        if(!comp_trg(mesh.get_target_set<0>(), vti))
          throw String("Vertex-Target-Indices refinement failure");

        // check edge-target indices
        Index eti[] =
        {
          2, 3, 12, 13, 15, 14, 6, 7
        };
        if(!comp_trg(mesh.get_target_set<1>(), eti))
          throw String("Edge-Target-Indices refinement failure");

      } // validate_refined_patch_edge_submesh_2d

    } // namespace TestAux
  } // namespace Geometry
} // namespace FEAST

#include <kernel/base_mesh/cell.hpp>

namespace FEAST
{
  namespace BaseMesh
  {
    // indices of start and end vertex of the four edges in a quad
    unsigned char Numbering::quad_edge_vertices[4][2] = {{0,1}, {2,3}, {0,2}, {1,3}};

    // index of the next vertex in a quad w.r.t. ccw ordering ([1,3,0,2])
    unsigned char Numbering::quad_next_vertex_ccw[4] = {1,3,0,2};

    // index of the previous vertex in a quad w.r.t. ccw ordering ([2,0,3,1])
    unsigned char Numbering::quad_previous_vertex_ccw[4] = {2,0,3,1};

    // index of the next edge in a quad w.r.t. ccw ordering ([3,2,0,1])
    unsigned char Numbering::quad_next_edge_ccw[4] = {3,2,0,1};

    // index of the previous edge in a quad w.r.t. ccw ordering ([2,3,1,0])
    unsigned char Numbering::quad_previous_edge_ccw[4] = {2,3,1,0};

    // indices of start and end vertex of the twelve edges in a hexa
    unsigned char Numbering::hexa_edge_vertices[12][2]
      = {{0,1}, {2,3}, {4,5}, {6,7},   {0,2}, {1,3}, {4,6}, {5,7},   {0,4}, {1,5}, {2,6}, {3,7}};

    // indices of the four vertices of the six faces in a hexa
    unsigned char Numbering::hexa_face_vertices[6][4]
      = {{0,1,2,3}, {4,5,6,7}, {0,1,4,5}, {2,3,6,7}, {0,2,4,6}, {1,3,5,7}};

    // indices of the four edges of the six faces in a hexa
    unsigned char Numbering::hexa_face_edges[6][4]
      = {{0,1,4,5}, {2,3,6,7}, {0,2,8,9}, {1,3,10,11}, {4,6,8,10}, {5,7,9,11}};

    // quad-to-quad mappings for vertices
    // V0:0123    V1:1302    V2:2031    V3:3210    V4:0213    V5:1032    V6:2301    V7:3120
    unsigned char Numbering::quad_to_quad_mappings_vertices[8][4] =
      {{0,1,2,3}, {1,3,0,2}, {2,0,3,1}, {3,2,1,0}, {0,2,1,3}, {1,0,3,2}, {2,3,0,1}, {3,1,2,0}};

    // quad-to-quad mappings for edges
    // E0:0123    E1:3201    E2:2310    E3:1032    E4:2301    E5:0132    E6:1023     E7:3210
    unsigned char Numbering::quad_to_quad_mappings_edges[8][4] =
      {{0,1,2,3}, {3,2,0,1}, {2,3,1,0}, {1,0,3,2}, {2,3,0,1}, {0,1,3,2}, {1,0,2,3}, {3,2,1,0}};

  } // namespace BaseMesh
} // namespace FEAST

#pragma once
#ifndef KERNEL_FOUNDATION_MESH_CONTROL_HPP
#define KERNEL_FOUNDATION_MESH_CONTROL_HPP

#include <kernel/foundation/mesh.hpp>
#include<kernel/geometry/conformal_mesh.hpp>

using namespace FEAST::Geometry;

namespace FEAST
{
  namespace Foundation
  {
    enum Dimensions
    {
      dim_1D = 1,
      dim_2D = 2,
      dim_3D = 3
    };

    template<Dimensions dim_>
    struct MeshControl
    {
    };

    template<>
    struct MeshControl<dim_1D>
    {
      template<typename SourceMeshType_>
      static typename SourceMeshType_::index_type_ fill(SourceMeshType_& mesh, typename SourceMeshType_::index_type_* target) //TODO const param
      {
        target[0] = mesh.get_topologies().at(0).get_topology().size();
        target[1] = mesh.get_topologies().at(1).get_topology().size();
      }

      template<typename SourceMeshType_, typename TargetMeshType_>
      static void fill(SourceMeshType_& source_mesh, TargetMeshType_& target_mesh)
      {
        //edge->Vertex
        typename TargetMeshType_::template IndexSet<1, 0>::Type& target_vertex_at_edge(target_mesh.template get_index_set<1, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_edge_vertex).size() ; ++i)
        {
          //get all adjacencies edge->Vertex from source mesh
          typename SourceMeshType_::topology_type_::storage_type_ source_vertex_at_edge_i(source_mesh.get_adjacent_polytopes(pl_edge, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_edge_i.size() ; ++j)
          {
            target_vertex_at_edge[i][j] = source_vertex_at_edge_i.at(j); //edge i, adjacent vertex j
          }
        }
      }
    };

    template<>
    struct MeshControl<dim_2D>
    {
      template<typename SourceMeshType_>
      static typename SourceMeshType_::index_type_ fill(SourceMeshType_& mesh, typename SourceMeshType_::index_type_* target)
      {
        target[0] = mesh.get_topologies().at(0).get_topology().size();
        target[1] = mesh.get_topologies().at(1).get_topology().size();
        target[2] = mesh.get_topologies().at(3).get_topology().size();
      }

      template<typename SourceMeshType_, typename TargetMeshType_>
      static void fill(SourceMeshType_& source_mesh, TargetMeshType_& target_mesh)
      {
        //edge->Vertex
        typename TargetMeshType_::template IndexSet<1, 0>::Type& target_vertex_at_edge(target_mesh.template get_index_set<1, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_edge_vertex).size() ; ++i)
        {
          //get all adjacencies edge->Vertex from source mesh
          typename SourceMeshType_::topology_type_::storage_type_ source_vertex_at_edge_i(source_mesh.get_adjacent_polytopes(pl_edge, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_edge_i.size() ; ++j)
          {
            target_vertex_at_edge[i][j] = source_vertex_at_edge_i.at(j); //edge i, adjacent vertex j
          }
        }

        //face->Vertex
        typename TargetMeshType_::template IndexSet<2, 0>::Type& target_vertex_at_face(target_mesh.template get_index_set<2, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_face_vertex).size() ; ++i)
        {
          //get all adjacencies face->Vertex from source mesh
          typename SourceMeshType_::topology_type_::storage_type_ source_vertex_at_face_i(source_mesh.get_adjacent_polytopes(pl_face, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_face_i.size() ; ++j)
          {
            target_vertex_at_face[i][j] = source_vertex_at_face_i.at(j); //face i, adjacent vertex j
          }
        }
      }
    };

    template<>
    struct MeshControl<dim_3D>
    {
      template<typename SourceMeshType_>
      static typename SourceMeshType_::index_type_ fill(SourceMeshType_& mesh, typename SourceMeshType_::index_type_* target)
      {
        target[0] = mesh.get_topologies().at(0).get_topology().size();
        target[1] = mesh.get_topologies().at(1).get_topology().size();
        target[2] = mesh.get_topologies().at(3).get_topology().size();
        target[5] = mesh.get_topologies().at(5).get_topology().size();
      }

      template<typename SourceMeshType_, typename TargetMeshType_>
      static void fill(SourceMeshType_& source_mesh, TargetMeshType_& target_mesh)
      {
        //edge->Vertex
        typename TargetMeshType_::template IndexSet<1, 0>::Type& target_vertex_at_edge(target_mesh.template get_index_set<1, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_edge_vertex).size() ; ++i)
        {
          //get all adjacencies edge->Vertex from source mesh
          typename SourceMeshType_::topology_type_::storage_type_ source_vertex_at_edge_i(source_mesh.get_adjacent_polytopes(pl_edge, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_edge_i.size() ; ++j)
          {
            target_vertex_at_edge[i][j] = source_vertex_at_edge_i.at(j); //edge i, adjacent vertex j
          }
        }

        //face->Vertex
        typename TargetMeshType_::template IndexSet<2, 0>::Type& target_vertex_at_face(target_mesh.template get_index_set<2, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_face_vertex).size() ; ++i)
        {
          //get all adjacencies face->Vertex from source mesh
          typename SourceMeshType_::topology_type_::storage_type_ source_vertex_at_face_i(source_mesh.get_adjacent_polytopes(pl_face, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_face_i.size() ; ++j)
          {
            target_vertex_at_face[i][j] = source_vertex_at_face_i.at(j); //face i, adjacent vertex j
          }
        }

        //polyhedron->Vertex
        typename TargetMeshType_::template IndexSet<3, 0>::Type& target_vertex_at_polyhedron(target_mesh.template get_index_set<3, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_polyhedron_vertex).size() ; ++i)
        {
          //get all adjacencies polyhedron->Vertex from source mesh
          typename SourceMeshType_::topology_type_::storage_type_ source_vertex_at_polyhedron_i(source_mesh.get_adjacent_polytopes(pl_polyhedron, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_polyhedron_i.size() ; ++j)
          {
            target_vertex_at_polyhedron[i][j] = source_vertex_at_polyhedron_i.at(j); //polyhedron i, adjacent vertex j
          }
        }
      }
    };
  }
}

#endif

#pragma once
#ifndef KERNEL_FOUNDATION_MESH_CONTROL_HPP
#define KERNEL_FOUNDATION_MESH_CONTROL_HPP

#include <kernel/foundation/base.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/halo.hpp>
#include<kernel/geometry/conformal_mesh.hpp>
#include<kernel/geometry/index_calculator.hpp>

using namespace FEAST::Geometry;

namespace FEAST
{
  namespace Foundation
  {
    template<Dimensions dim_>
    struct MeshControl
    {
    };

    template<>
    struct MeshControl<dim_1D>
    {
      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class MeshType_>
      static void fill_sizes(const MeshType_<a_, b_, c_>& mesh,
                             typename MeshType_<a_, b_, c_>::index_type_* target)
      {
        target[0] = (typename MeshType_<a_, b_, c_>::index_type_)(mesh.get_topologies().at(0).get_topology().size());
        target[1] = (typename MeshType_<a_, b_, c_>::index_type_)(mesh.get_topologies().at(1).get_topology().size());
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_,
        typename TargetMeshType_>
      static void fill_adjacencies(const SourceMeshType_<a_, b_, c_>& source_mesh,
                                   TargetMeshType_& target_mesh)
      {
        //edge->Vertex
        typename TargetMeshType_::template IndexSet<1, 0>::Type& target_vertex_at_edge(target_mesh.template get_index_set<1, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_edge_vertex).size() ; ++i)
        {
          //get all adjacencies edge->Vertex from source mesh
          typename SourceMeshType_<a_, b_, c_>::topology_type_::storage_type_ source_vertex_at_edge_i(source_mesh.get_adjacent_polytopes(pl_edge, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_edge_i.size() ; ++j)
          {
            target_vertex_at_edge[i][j] = source_vertex_at_edge_i.at(j); //edge i, adjacent vertex j
          }
        }
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class TargetMeshType_,
        typename SourceMeshType_>
      static void fill_adjacencies(SourceMeshType_& geo_mesh,
                                   TargetMeshType_<a_, b_, c_>& found_mesh)
      {
        typename SourceMeshType_::template IndexSet<1, 0>::Type& source_vertex_at_edge(geo_mesh.template get_index_set<1, 0>());

        //create entities
        for(Index i(0) ; i < geo_mesh.get_num_entities(0) ; ++i)
        {
          found_mesh.add_polytope(pl_vertex);
        }

        for(Index i(0) ; i < geo_mesh.get_num_entities(1) ; ++i)
        {
          found_mesh.add_polytope(pl_edge);
        }

        //add adjacencies edge->vertex
        for(Index i(0) ; i < source_vertex_at_edge.get_num_entities() ; ++i)
          for(Index j(0) ; j < (Index)source_vertex_at_edge.get_num_indices() ; ++j)
          {
            found_mesh.add_adjacency(pl_edge, pl_vertex, i, source_vertex_at_edge[i][j]);
          }
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_,
        typename TargetMeshType_,
        typename AttributeType_>
      static void fill_vertex_sets(const SourceMeshType_<a_, b_, c_>& source_mesh,
                                   TargetMeshType_& target_mesh, const AttributeType_& attr)
      {
        typename TargetMeshType_::VertexSetType& vertex_coord_tuples(target_mesh.get_vertex_set());
        for(Index i(0) ; i < source_mesh.get_topologies().at(0).get_topology().size() ; ++i)
        {
          vertex_coord_tuples[i][0] = attr.get_data().at(i);
        }
      }

      ///geo -> found, overwrite version
      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class TargetMeshType_,
        typename SourceMeshType_,
        typename AttributeType_>
      static void fill_vertex_sets(SourceMeshType_& geo_mesh,
                                   const TargetMeshType_<a_, b_, c_>& found_mesh,
                                   AttributeType_& attr)
      {
        ///TODO check attribute registration status with target mesh
        typename SourceMeshType_::VertexSetType& vertex_coord_tuples(geo_mesh.get_vertex_set());

        attr.get_data().clear();
        for(Index i(0) ; i < found_mesh.get_topologies().at(0).get_topology().size() ; ++i)
        {
          attr.get_data().push_back(vertex_coord_tuples[i][0]);
        }
      }
    };

    template<>
    struct MeshControl<dim_2D>
    {
      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_>
      static void fill_sizes(const SourceMeshType_<a_, b_, c_>& mesh,
                             typename SourceMeshType_<a_, b_, c_>::index_type_* target)
      {
        typedef typename SourceMeshType_<a_, b_, c_>::index_type_ IndexType;
        target[0] = IndexType(mesh.get_topologies().at(0).get_topology().size());
        target[1] = IndexType(mesh.get_topologies().at(1).get_topology().size());
        target[2] = IndexType(mesh.get_topologies().at(3).get_topology().size());
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_,
        typename TargetMeshType_>
      static void fill_adjacencies(const SourceMeshType_<a_, b_, c_>& source_mesh,
                                   TargetMeshType_& target_mesh)
      {
        //edge->Vertex
        typename TargetMeshType_::template IndexSet<1, 0>::Type& target_vertex_at_edge(target_mesh.template get_index_set<1, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_edge_vertex).size() ; ++i)
        {
          //get all adjacencies edge->Vertex from source mesh
          typename SourceMeshType_<a_, b_, c_>::topology_type_::storage_type_ source_vertex_at_edge_i(source_mesh.get_adjacent_polytopes(pl_edge, pl_vertex, i));

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
          typename SourceMeshType_<a_, b_, c_>::topology_type_::storage_type_ source_vertex_at_face_i(source_mesh.get_adjacent_polytopes(pl_face, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_face_i.size() ; ++j)
          {
            target_vertex_at_face[i][j] = source_vertex_at_face_i.at(j); //face i, adjacent vertex j
          }
        }

        //face->edge
        typedef typename TargetMeshType_::ShapeType face_type;
        typedef typename Shape::FaceTraits<face_type, face_type::dimension - 1>::ShapeType edge_type;

        IndexTree<edge_type> it(target_mesh.get_num_entities(0));
        it.parse(target_vertex_at_edge);


        typename TargetMeshType_::template IndexSet<2, 1>::Type& target_edge_at_face(target_mesh.template get_index_set<2, 1>());
        IndexCalculator<face_type, face_type::dimension - 1>::compute(it, target_vertex_at_face, target_edge_at_face);
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class TargetMeshType_,
        typename SourceMeshType_>
      static void fill_adjacencies(SourceMeshType_& geo_mesh,
                                   TargetMeshType_<a_, b_, c_>& found_mesh)
      {
        typename SourceMeshType_::template IndexSet<1, 0>::Type& source_vertex_at_edge(geo_mesh.template get_index_set<1, 0>());
        typename SourceMeshType_::template IndexSet<2, 0>::Type& source_vertex_at_face(geo_mesh.template get_index_set<2, 0>());

        //create entities
        for(Index i(0) ; i < geo_mesh.get_num_entities(0) ; ++i)
        {
          found_mesh.add_polytope(pl_vertex);
        }

        for(Index i(0) ; i < geo_mesh.get_num_entities(1) ; ++i)
        {
          found_mesh.add_polytope(pl_edge);
        }

        for(Index i(0) ; i < geo_mesh.get_num_entities(2) ; ++i)
        {
          found_mesh.add_polytope(pl_face);
        }

        //add adjacencies edge->vertex
        for(Index i(0) ; i < source_vertex_at_edge.get_num_entities() ; ++i)
          for(Index j(0) ; j < (Index)source_vertex_at_edge.get_num_indices() ; ++j)
          {
            found_mesh.add_adjacency(pl_edge, pl_vertex, i, source_vertex_at_edge[i][j]);
          }

        //add adjacencies face->vertex
        for(Index i(0) ; i < source_vertex_at_face.get_num_entities() ; ++i)
          for(Index j(0) ; j < (Index)source_vertex_at_face.get_num_indices() ; ++j)
          {
            found_mesh.add_adjacency(pl_face, pl_vertex, i, source_vertex_at_face[i][j]);
          }
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_,
        typename TargetMeshType_,
        typename AttributeType1_,
        typename AttributeType2_>
      static void fill_vertex_sets(const SourceMeshType_<a_, b_, c_>& source_mesh,
                                   TargetMeshType_& target_mesh,
                                   const AttributeType1_& attr_0,
                                   const AttributeType2_& attr_1)
      {
        typename TargetMeshType_::VertexSetType& vertex_coord_tuples(target_mesh.get_vertex_set());
        for(Index i(0) ; i < source_mesh.get_topologies().at(0).get_topology().size() ; ++i)
        {
          vertex_coord_tuples[i][0] = attr_0.get_data().at(i);
          vertex_coord_tuples[i][1] = attr_1.get_data().at(i);
        }
      }

      ///geo -> found, overwrite version
      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class TargetMeshType_,
        typename SourceMeshType_,
        typename AttributeType1_,
        typename AttributeType2_>
      static void fill_vertex_sets(SourceMeshType_& geo_mesh,
                                   const TargetMeshType_<a_, b_, c_>& found_mesh,
                                   AttributeType1_& attr1,
                                   AttributeType2_& attr2)
      {
        ///TODO check attribute registration status with target mesh
        typename SourceMeshType_::VertexSetType& vertex_coord_tuples(geo_mesh.get_vertex_set());

        attr1.get_data().clear();
        attr2.get_data().clear();
        for(Index i(0) ; i < found_mesh.get_topologies().at(0).get_topology().size() ; ++i)
        {
          attr1.get_data().push_back(vertex_coord_tuples[i][0]);
          attr2.get_data().push_back(vertex_coord_tuples[i][1]);
        }
      }
    };

    template<>
    struct MeshControl<dim_3D>
    {
      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_>
      static void fill_sizes(const SourceMeshType_<a_, b_, c_>& mesh,
                             typename SourceMeshType_<a_, b_, c_>::index_type_* target)
      {
        target[0] = (typename SourceMeshType_<a_, b_, c_>::index_type_)(mesh.get_topologies().at(0).get_topology().size());
        target[1] = (typename SourceMeshType_<a_, b_, c_>::index_type_)(mesh.get_topologies().at(1).get_topology().size());
        target[2] = (typename SourceMeshType_<a_, b_, c_>::index_type_)(mesh.get_topologies().at(3).get_topology().size());
        target[3] = (typename SourceMeshType_<a_, b_, c_>::index_type_)(mesh.get_topologies().at(5).get_topology().size());
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_,
        typename TargetMeshType_>
      static void fill_adjacencies(const SourceMeshType_<a_, b_, c_>& source_mesh,
                                   TargetMeshType_& target_mesh)
      {
        //edge->Vertex
        typename TargetMeshType_::template IndexSet<1, 0>::Type& target_vertex_at_edge(target_mesh.template get_index_set<1, 0>());
        for(Index i(0) ; i < source_mesh.get_topologies().at(ipi_edge_vertex).size() ; ++i)
        {
          //get all adjacencies edge->Vertex from source mesh
          typename SourceMeshType_<a_, b_, c_>::topology_type_::storage_type_ source_vertex_at_edge_i(source_mesh.get_adjacent_polytopes(pl_edge, pl_vertex, i));

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
          typename SourceMeshType_<a_, b_, c_>::topology_type_::storage_type_ source_vertex_at_face_i(source_mesh.get_adjacent_polytopes(pl_face, pl_vertex, i));

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
          typename SourceMeshType_<a_, b_, c_>::topology_type_::storage_type_ source_vertex_at_polyhedron_i(source_mesh.get_adjacent_polytopes(pl_polyhedron, pl_vertex, i));

          for(Index j(0) ; j < source_vertex_at_polyhedron_i.size() ; ++j)
          {
            target_vertex_at_polyhedron[i][j] = source_vertex_at_polyhedron_i.at(j); //polyhedron i, adjacent vertex j
          }
        }

        //poly->face
        typedef typename TargetMeshType_::ShapeType poly_type;
        typedef typename Shape::FaceTraits<poly_type, poly_type::dimension - 1>::ShapeType face_type;
        typedef typename Shape::FaceTraits<poly_type, poly_type::dimension - 2>::ShapeType edge_type;

        IndexTree<face_type> it_face(target_mesh.get_num_entities(0));
        IndexTree<edge_type> it_edge(target_mesh.get_num_entities(0));

        it_face.parse(target_vertex_at_face);
        it_edge.parse(target_vertex_at_edge);


        typename TargetMeshType_::template IndexSet<2, 1>::Type& target_edge_at_face(target_mesh.template get_index_set<2, 1>());
        typename TargetMeshType_::template IndexSet<3, 1>::Type& target_edge_at_poly(target_mesh.template get_index_set<3, 1>());
        typename TargetMeshType_::template IndexSet<3, 2>::Type& target_face_at_poly(target_mesh.template get_index_set<3, 2>());

        IndexCalculator<poly_type, poly_type::dimension - 1>::compute(it_face, target_vertex_at_polyhedron, target_face_at_poly);
        IndexCalculator<face_type, face_type::dimension - 1>::compute(it_edge, target_vertex_at_face, target_edge_at_face);
        IndexCalculator<poly_type, poly_type::dimension - 2>::compute(it_edge, target_vertex_at_polyhedron, target_edge_at_poly);
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class TargetMeshType_,
        typename SourceMeshType_>
      static void fill_adjacencies(SourceMeshType_& geo_mesh,
                                   TargetMeshType_<a_, b_, c_>& found_mesh)
      {
        typename SourceMeshType_::template IndexSet<1, 0>::Type& source_vertex_at_edge(geo_mesh.template get_index_set<1, 0>());
        typename SourceMeshType_::template IndexSet<2, 0>::Type& source_vertex_at_face(geo_mesh.template get_index_set<2, 0>());
        typename SourceMeshType_::template IndexSet<3, 0>::Type& source_vertex_at_poly(geo_mesh.template get_index_set<3, 0>());

        //create entities
        for(Index i(0) ; i < geo_mesh.get_num_entities(0) ; ++i)
        {
          found_mesh.add_polytope(pl_vertex);
        }

        for(Index i(0) ; i < geo_mesh.get_num_entities(1) ; ++i)
        {
          found_mesh.add_polytope(pl_edge);
        }

        for(Index i(0) ; i < geo_mesh.get_num_entities(2) ; ++i)
        {
          found_mesh.add_polytope(pl_face);
        }

        for(Index i(0) ; i < geo_mesh.get_num_entities(3) ; ++i)
        {
          found_mesh.add_polytope(pl_polyhedron);
        }

        //add adjacencies edge->vertex
        for(Index i(0) ; i < source_vertex_at_edge.get_num_entities() ; ++i)
          for(Index j(0) ; j < (Index)source_vertex_at_edge.get_num_indices() ; ++j)
          {
            found_mesh.add_adjacency(pl_edge, pl_vertex, i, source_vertex_at_edge[i][j]);
          }

        //add adjacencies face->vertex
        for(Index i(0) ; i < source_vertex_at_face.get_num_entities() ; ++i)
          for(Index j(0) ; j < (Index)source_vertex_at_face.get_num_indices() ; ++j)
          {
            found_mesh.add_adjacency(pl_face, pl_vertex, i, source_vertex_at_face[i][j]);
          }

        //add adjacencies poly->vertex
        for(Index i(0) ; i < source_vertex_at_poly.get_num_entities() ; ++i)
          for(Index j(0) ; j < (Index)source_vertex_at_poly.get_num_indices() ; ++j)
          {
            found_mesh.add_adjacency(pl_polyhedron, pl_vertex, i, source_vertex_at_poly[i][j]);
          }
      }

      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class SourceMeshType_,
        typename TargetMeshType_,
        typename AttributeType1_,
        typename AttributeType2_,
        typename AttributeType3_>
      static void fill_vertex_sets(const SourceMeshType_<a_, b_, c_>& source_mesh,
                                   TargetMeshType_& target_mesh,
                                   const AttributeType1_& attr_0,
                                   const AttributeType2_& attr_1,
                                   const AttributeType3_& attr_2)
      {
        typename TargetMeshType_::VertexSetType& vertex_coord_tuples(target_mesh.get_vertex_set());
        for(Index i(0) ; i < source_mesh.get_topologies().at(0).get_topology().size() ; ++i)
        {
          vertex_coord_tuples[i][0] = attr_0.get_data().at(i);
          vertex_coord_tuples[i][1] = attr_1.get_data().at(i);
          vertex_coord_tuples[i][2] = attr_2.get_data().at(i);
        }
      }

      ///geo -> found, overwrite version
      template<
        typename a_,
        typename b_,
        template <typename, typename> class c_,
        template<
          typename,
          typename,
          template <typename, typename> class> class TargetMeshType_,
        typename SourceMeshType_,
        typename AttributeType1_,
        typename AttributeType2_,
        typename AttributeType3_>
      static void fill_vertex_sets(SourceMeshType_& geo_mesh,
                                   const TargetMeshType_<a_, b_, c_>& found_mesh,
                                   AttributeType1_& attr1,
                                   AttributeType2_& attr2,
                                   AttributeType3_& attr3)
      {
        ///TODO check attribute registration status with target mesh
        typename SourceMeshType_::VertexSetType& vertex_coord_tuples(geo_mesh.get_vertex_set());

        attr1.get_data().clear();
        attr2.get_data().clear();
        attr3.get_data().clear();
        for(Index i(0) ; i < found_mesh.get_topologies().at(0).get_topology().size() ; ++i)
        {
          attr1.get_data().push_back(vertex_coord_tuples[i][0]);
          attr2.get_data().push_back(vertex_coord_tuples[i][1]);
          attr3.get_data().push_back(vertex_coord_tuples[i][2]);
        }
      }
    };
  }
}

#endif

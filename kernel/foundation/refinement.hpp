#pragma once
#ifndef FOUNDATION_GUARD_REFINEMENT_HPP
#define FOUNDATION_GUARD_REFINEMENT_HPP 1

#include<kernel/foundation/attribute.hpp>
#include<kernel/foundation/mesh.hpp>
#include<kernel/foundation/halo.hpp>

using namespace FEAST;

namespace FEAST
{
  namespace Foundation
  {
    enum MeshRefinementTypes
    {
      mrt_standard = 0 //standard refinement depending on cell type TODO assumes QUAD
    };

    enum HaloRefinementTypes
    {
      hrt_refine = 0, //halos are refined -> overlaps doubled
      hrt_norefine    //halos remain unrefined -> overlaps remain unmodified
    };

    template<typename Tag_,
             typename Arch_,
             Dimensions dim_,
             PolytopeLevels pl_,
             MeshRefinementTypes mrt_,
             HaloRefinementTypes hrt_>
    struct PolytopeRefinement
    {
    };

    ///enforce interpretation of arbitrary mesh as 'residing in 1D space'
    template<typename Tag_,
             typename Arch_,
             HaloRefinementTypes hrt_>
    struct PolytopeRefinement<Tag_,
                              Arch_,
                              dim_1D,
                              pl_edge,
                              mrt_standard,
                              hrt_>
    {
      ///single edge w. vertex data, halos
      template<
        typename Dim_,
        typename t_,
        template <typename, typename> class os_,
        template <typename, typename, template<typename, typename> class > class MeshType_,
        typename DataType_>
      static void execute(
                          Index edge_index,
                          MeshType_<Dim_, t_, os_>& target,
                          os_<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> > > >* halos,
                          os_<Attribute<DataType_, os_>, std::allocator<Attribute<DataType_, os_> > >& origin_coords)
      {
        //get adjacent vertices *____*
        typename MeshType_<Dim_, t_, os_>::topology_type_::storage_type_ vertex_at_edge(target.get_adjacent_polytopes(pl_edge, pl_vertex, edge_index));

        //create new vertex *--*--*
        //                  a  j  b
        Index new_vertex_id(target.get_topologies().at(ipi_vertex_edge).size());
        target.add_polytope(pl_vertex);

        //add new coordinate
        const typename Attribute<DataType_, os_>::data_type_ x_0(origin_coords.at(0).at(vertex_at_edge.at(0)));
        const typename Attribute<DataType_, os_>::data_type_ x_1(origin_coords.at(0).at(vertex_at_edge.at(1)));
        if(x_1 > x_0)
          origin_coords.at(0).push_back(x_0 + DataType_(1./2.)*(x_1 - x_0));
        else
          origin_coords.at(0).push_back(x_1 + DataType_(1./2.)*(x_0 - x_1));

        //make old edges' second adjacency be the new vertex  *__*--*
        //                                                    a  j  b
        target.get_topologies().at(ipi_edge_vertex).at(edge_index).at(1) = new_vertex_id;
        target.get_topologies().at(ipi_vertex_edge).at(new_vertex_id).push_back(edge_index);

        //remove edge number from second vertex' adjacency list *__*  *
        //                                                      a  j  b
        target.get_topologies().at(ipi_vertex_edge).at(vertex_at_edge.at(1)).erase(std::find(
              target.get_topologies().at(ipi_vertex_edge).at(vertex_at_edge.at(1)).begin(),
              target.get_topologies().at(ipi_vertex_edge).at(vertex_at_edge.at(1)).end(),
              edge_index));

        //create additional edge
        Index new_edge_id(target.get_topologies().at(ipi_edge_vertex).size());
        target.add_polytope(pl_edge);

        //connect new edge  *__*__*
        //                  a  j  b
        target.add_adjacency(pl_edge, pl_vertex, new_edge_id, new_vertex_id);
        target.add_adjacency(pl_edge, pl_vertex, new_edge_id, vertex_at_edge.at(1));
        //target.add_adjacency(pl_vertex, pl_edge, new_vertex_id, edge_index);

        if(hrt_ == hrt_refine)
        {
          //refine halo edges
          if(halos != nullptr)
          {
            for(Index i(0) ; i < halos->size() ; ++i)
            {
              //only do this for edge halos with delta > 0
              if(halos->at(i)->get_level() == pl_edge && halos->at(i)->get_overlap() > 0)
              {
                typename HaloBase<MeshType_<Dim_, t_, os_>, os_>::compound_storage_type_::iterator iter(std::find(halos->at(i)->get_elements().begin(), halos->at(i)->get_elements().end(), edge_index));
                if(iter != halos->at(i)->get_elements().end())
                {
                  if(iter + 1 == halos->at(i)->get_elements().end())
                  {
                    halos->at(i)->get_elements().push_back(new_edge_id);
                  }
                  else
                  {
                    halos->at(i)->get_elements().insert(iter + 1, new_edge_id); //insertion after
                  }
                }
              }
            }
          }
        }
      }
    };

    ///enforce interpretation of arbitrary mesh as 'residing in 2D space'
    template<typename Tag_,
             typename Arch_,
             HaloRefinementTypes hrt_>
    struct PolytopeRefinement<Tag_,
                              Arch_,
                              dim_2D,
                              pl_edge,
                              mrt_standard,
                              hrt_>
    {
      ///single edge w. vertex data, halos
      template<
        typename Dim_,
        typename t_,
        template <typename, typename> class os_,
        template <typename, typename, template<typename, typename> class > class MeshType_,
        typename DataType_>
      static void execute(
                          Index edge_index,
                          MeshType_<Dim_, t_, os_>& target,
                          os_<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> > > >* halos,
                          os_<Attribute<DataType_, os_>, std::allocator<Attribute<DataType_, os_> > >& origin_coords)
      {
        //get adjacent vertices *____*
        typename MeshType_<Dim_, t_, os_>::topology_type_::storage_type_ vertex_at_edge(target.get_adjacent_polytopes(pl_edge, pl_vertex, edge_index));

        //store needed old adjacency situation temporarily
        typename MeshType_<Dim_, t_, os_>::topology_type_::storage_type_ vertex_at_edge_old(target.get_adjacent_polytopes(pl_edge, pl_vertex, edge_index));

        //create new vertex *--*--*
        //                  a  j  b
        Index new_vertex_id(target.get_topologies().at(ipi_vertex_edge).size());
        target.add_polytope(pl_vertex);

        //make old edges' second adjacency be the new vertex  *__*--*
        //                                                    a  j  b
        target.get_topologies().at(ipi_edge_vertex).at(edge_index).at(1) = new_vertex_id;
        target.get_topologies().at(ipi_vertex_edge).at(new_vertex_id).push_back(edge_index);

        //remove edge number from second vertex' adjacency list *__*  *
        //                                                      a  j  b
        target.get_topologies().at(ipi_vertex_edge).at(vertex_at_edge.at(1)).erase(std::find(
              target.get_topologies().at(ipi_vertex_edge).at(vertex_at_edge.at(1)).begin(),
              target.get_topologies().at(ipi_vertex_edge).at(vertex_at_edge.at(1)).end(),
              edge_index));

        //create additional edge
        Index new_edge_id(target.get_topologies().at(ipi_edge_vertex).size());
        target.add_polytope(pl_edge);

        //connect new edge  *__*__*
        //                  a  j  b
        target.add_adjacency(pl_edge, pl_vertex, new_edge_id, new_vertex_id);
        target.add_adjacency(pl_edge, pl_vertex, new_edge_id, vertex_at_edge.at(1));
        //target.add_adjacency(pl_vertex, pl_edge, new_vertex_id, edge_index);

        if(hrt_ == hrt_refine)
        {
          //refine halo edges
          if(halos != nullptr)
          {
            for(Index i(0) ; i < halos->size() ; ++i)
            {
              //only do this for edge halos with delta > 0
              if(halos->at(i)->get_level() == pl_edge) //&&halos->at(i)->get_overlap() > 0)
              {
                typename HaloBase<MeshType_<Dim_, t_, os_>, os_>::compound_storage_type_::iterator iter(std::find(halos->at(i)->get_elements().begin(), halos->at(i)->get_elements().end(), edge_index));
                if(iter != halos->at(i)->get_elements().end())
                {
                  if(iter + 1 == halos->at(i)->get_elements().end())
                  {
                    halos->at(i)->get_elements().push_back(new_edge_id);
                  }
                  else
                  {
                    halos->at(i)->get_elements().insert(iter + 1, new_edge_id); //insertion after
                  }
                }
              }
            }
          }
        }

        if(Dim_::required_num_topologies == rnt_2D)
        {
          //set face adjacencies for new vertex
          typename MeshType_<Dim_, t_, os_>::topology_type_::storage_type_ faces_at_start_vertex_old(target.get_adjacent_polytopes(pl_vertex, pl_face, vertex_at_edge_old.at(0)));
          typename MeshType_<Dim_, t_, os_>::topology_type_::storage_type_ faces_at_end_vertex_old(target.get_adjacent_polytopes(pl_vertex, pl_face, vertex_at_edge_old.at(1)));
          for(Index i(0) ; i < faces_at_end_vertex_old.size() ; ++i)
            if(std::find(faces_at_start_vertex_old.begin(), faces_at_start_vertex_old.end(), faces_at_end_vertex_old.at(i)) != faces_at_start_vertex_old.end())
              target.add_adjacency(pl_vertex,
                                   pl_face,
                                   target.get_topologies().at(ipi_vertex_face).size() - 1,
                                   faces_at_end_vertex_old.at(i));
        }

        //add new coordinates
        const typename Attribute<DataType_, os_>::data_type_ x_0(origin_coords.at(0).at(vertex_at_edge.at(0)));
        const typename Attribute<DataType_, os_>::data_type_ x_1(origin_coords.at(0).at(vertex_at_edge.at(1)));

        if(std::abs(x_1 - x_0) <= std::numeric_limits<DataType_>::epsilon())
          origin_coords.at(0).push_back(x_0);
        else if(x_1 > x_0)
          origin_coords.at(0).push_back(x_0 + DataType_(1./2.)*(x_1 - x_0));
        else
          origin_coords.at(0).push_back(x_1 + DataType_(1./2.)*(x_0 - x_1));

        const typename Attribute<DataType_, os_>::data_type_ y_0(origin_coords.at(1).at(vertex_at_edge.at(0)));
        const typename Attribute<DataType_, os_>::data_type_ y_1(origin_coords.at(1).at(vertex_at_edge.at(1)));

        if(std::abs(y_1 - y_0) <= std::numeric_limits<DataType_>::epsilon())
          origin_coords.at(1).push_back(y_0);
        else if(x_1 > x_0) //must use same direction as x !!
          origin_coords.at(1).push_back(y_0 + DataType_(1./2.)*(y_1 - y_0));
        else
          origin_coords.at(1).push_back(y_1 + DataType_(1./2.)*(y_0 - y_1));

      }
    };

    template<typename Tag_,
             typename Arch_,
             MeshRefinementTypes mrt_ = mrt_standard,
             HaloRefinementTypes hrt_ = hrt_refine>
    struct Refinement
    {
      ///1D
      ///complete w. vertex data, halos
      template<
        typename t_,
        template <typename, typename> class os_,
        template <typename, typename, template<typename, typename> class > class MeshType_,
        typename DataType_>
      static void execute(
                          MeshType_<Dim1D, t_, os_>& origin,
                          os_<std::shared_ptr<HaloBase<MeshType_<Dim1D, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_<Dim1D, t_, os_>, os_> > > >* halos,
                          os_<Attribute<DataType_, os_>, std::allocator<Attribute<DataType_, os_> > >& origin_coords)
      {
        const Index num_edges(origin.get_topologies().at(ipi_edge_vertex).size());
        for(Index i(0) ; i < num_edges; ++i)
        {
          PolytopeRefinement<Tag_,
                             Arch_,
                             dim_1D,
                             pl_edge,
                             mrt_,
                             hrt_>::execute(i, origin, halos, origin_coords);

        ///TODO convert Halo<k> to Halo<2k> if hrt_refine
        }
      }

      ///2D
      ///complete w. vertex data
      template<
        typename t_,
        template <typename, typename> class os_,
        template <typename, typename, template<typename, typename> class > class MeshType_,
        typename DataType_>
      static void execute(
                          MeshType_<Dim2D, t_, os_>& origin,
                          os_<std::shared_ptr<HaloBase<MeshType_<Dim2D, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_<Dim2D, t_, os_>, os_> > > >* halos,
                          os_<Attribute<DataType_, os_>, std::allocator<Attribute<DataType_, os_> > >& origin_coords)
      {
        ///TODO assumes QUAD, needs coarse copy
        MeshType_<Dim2D, t_, os_> coarse(origin);

        //refine all edges first (sets all face/vertex adjacencies correctly)
        const Index num_edges(origin.get_topologies().at(ipi_edge_vertex).size());
        for(Index i(0) ; i < num_edges; ++i)
        {
          PolytopeRefinement<Tag_,
                             Arch_,
                             dim_2D,
                             pl_edge,
                             mrt_,
                             hrt_>::execute(
                                           i,
                                           origin,
                                           halos,
                                           origin_coords
                                           );
        }

        //for all faces, finalize the refinement
        typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ vertices_added;
        const Index num_faces(origin.get_topologies().at(ipi_face_vertex).size());
        for(Index i(0) ; i < num_faces; ++i)
        {
          //get all adjacent vertices
          typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ vertex_at_face(origin.get_adjacent_polytopes(pl_face, pl_vertex, i));
          typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ vertexold_at_face(coarse.get_adjacent_polytopes(pl_face, pl_vertex, i));
          //sort
          std::sort(vertexold_at_face.begin(), vertexold_at_face.end());

          //polygonal sorting
          typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ vertex_at_polygon;
          vertex_at_polygon.push_back(vertexold_at_face.at(0));
          typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ edgeold_at_face(coarse.get_adjacent_polytopes(pl_face, pl_edge, i));

          while(vertex_at_polygon.size() < 4) //QUAD
          {
            Index current(vertex_at_polygon.at(vertex_at_polygon.size() - 1));
            for(Index j(0) ; j < edgeold_at_face.size() ; ++j)
            {
              typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ vertex_at_edge_current(coarse.get_adjacent_polytopes(pl_edge, pl_vertex, edgeold_at_face.at(j)));
              //outgoing edge?
              if(std::find(vertex_at_edge_current.begin(), vertex_at_edge_current.end(), current) != vertex_at_edge_current.end())
              {
                //target selection
                Index next(Index(vertex_at_edge_current.at(Index(vertex_at_edge_current.begin() + 1 - std::find(vertex_at_edge_current.begin(), vertex_at_edge_current.end(), current)))));
                //not yet added?
                if(std::find(vertex_at_polygon.begin(), vertex_at_polygon.end(), next) == vertex_at_polygon.end())
                {
                  vertex_at_polygon.push_back(next);
                  break;
                }
              }
            }
          }

          //for all old vertices: use them as starting point for creation of new faces
          for(Index j(0) ; j < 4 ; ++j) //QUAD!
          {
            //get all vertex neighbours
            typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ vertex_neighbours(origin.get_adjacent_polytopes(pl_vertex, pl_vertex, vertexold_at_face.at(j)));
            std::sort(vertex_neighbours.begin(), vertex_neighbours.end());

            typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_ result(100);
            typename MeshType_<Dim2D, t_, os_>::topology_type_::storage_type_::iterator it;
            it = std::set_intersection(vertex_at_face.begin(), vertex_at_face.end(), vertex_neighbours.begin(), vertex_neighbours.end(), result.begin());
            result.resize(Index(it - result.begin()));

            if(std::find(result.begin(), result.end(), vertexold_at_face.at(j)) != result.end())
              result.erase(std::find(result.begin(), result.end(), vertexold_at_face.at(j))); //remove self

            //use old face number for first new face
            if(j == 0)
            {
              //remove all occurances of i in the vertex-face lists of all adjacent vertices
              for(Index k(0) ; k < vertex_at_face.size() ; ++k)
                origin.get_topologies().at(ipi_vertex_face).at(k).erase(std::find(origin.get_topologies().at(ipi_vertex_face).at(k).begin(), origin.get_topologies().at(ipi_vertex_face).at(k).end(), i));

              //clear the face's vertex list
              origin.get_topologies().at(ipi_face_vertex).at(i).clear();

              //create new vertex and add all vertices
              origin.add_polytope(pl_vertex); //create
              origin.add_adjacency(pl_face, pl_vertex, i, origin.get_topologies().at(ipi_vertex_face).size() - 1); //new one
              origin.add_adjacency(pl_face, pl_vertex, i, vertexold_at_face.at(j)); //starting vertex
              for(Index k(0) ; k < result.size() ; ++k) //result vertices
                origin.add_adjacency(pl_face, pl_vertex, i, result.at(k));

              //create and connect new edges
              for(Index k(0) ; k < result.size() ; ++k)
              {
                origin.add_polytope(pl_edge);
                origin.add_adjacency(pl_edge, pl_vertex, origin.get_topologies().at(ipi_edge_vertex).size() - 1, result.at(k));
                vertices_added.push_back(result.at(k));
                origin.add_adjacency(pl_edge, pl_vertex, origin.get_topologies().at(ipi_edge_vertex).size() - 1, origin.get_topologies().at(ipi_vertex_face).size() - 1);
              }
            }
            else
            {
              //create new face
              origin.add_polytope(pl_face);

              //add adjacencies
              origin.add_adjacency(pl_face, pl_vertex, origin.get_topologies().at(ipi_face_vertex).size() - 1, origin.get_topologies().at(ipi_vertex_face).size() - 1); //new one
              origin.add_adjacency(pl_face, pl_vertex, origin.get_topologies().at(ipi_face_vertex).size() - 1, vertexold_at_face.at(j)); //starting vertex
              for(Index k(0) ; k < result.size() ; ++k) //result vertices
                origin.add_adjacency(pl_face, pl_vertex, origin.get_topologies().at(ipi_face_vertex).size() - 1, result.at(k));

              //create and connect new edges
              for(Index k(0) ; k < result.size() ; ++k)
              {
                if(std::find(vertices_added.begin(), vertices_added.end(), result.at(k)) == vertices_added.end())
                {
                  origin.add_polytope(pl_edge);
                  origin.add_adjacency(pl_edge, pl_vertex, origin.get_topologies().at(ipi_edge_vertex).size() - 1, result.at(k));
                  vertices_added.push_back(result.at(k));
                  origin.add_adjacency(pl_edge, pl_vertex, origin.get_topologies().at(ipi_edge_vertex).size() - 1, origin.get_topologies().at(ipi_vertex_face).size() - 1);
                }
              }
            }
          }

          //finally: for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
          if(hrt_ == hrt_refine)
          {
            //refine halo faces
            if(halos != nullptr)
            {
              for(Index hi(0) ; hi < halos->size() ; ++hi)
              {
                //only do this for face halos with delta > 0
                if(halos->at(hi)->get_level() == pl_face && halos->at(hi)->get_overlap() > 0)
                {
                  typename HaloBase<MeshType_<Dim2D, t_, os_>, os_>::compound_storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                  if(iter != halos->at(hi)->get_elements().end())
                  {
                    if(iter + 1 == halos->at(hi)->get_elements().end())
                    {
                      halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 3); ///TODO QUAD
                      halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 2); ///TODO QUAD
                      halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1); ///TODO QUAD
                    }
                    else
                    {
                      halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 3); //insertion after ///TODO QUAD
                      halos->at(hi)->get_elements().insert(iter + 2, origin.get_topologies().at(ipi_face_vertex).size() - 2); //insertion after ///TODO QUAD
                      halos->at(hi)->get_elements().insert(iter + 3, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                    }
                  }
                }
              }
            }
          }

          //determine new vertex coordinates for mid-point
          //area
          DataType_ area(0), c_x(0), c_y(0);
          for(Index l(0); l < vertex_at_polygon.size() ; ++l)
          {
            Index lp1(l == vertex_at_polygon.size() - 1 ? 0 : l + 1);
            area += origin_coords.at(0).at(vertex_at_polygon.at(l)) * origin_coords.at(1).at(vertex_at_polygon.at(lp1)) - origin_coords.at(0).at(vertex_at_polygon.at(lp1)) * origin_coords.at(1).at(vertex_at_polygon.at(l));
          }
          area =  DataType_(1./2.) * std::abs(area);

          //area multiplier
          area = DataType_(1) / (DataType_(6) * area);

          //coords
          for(Index l(0); l < vertex_at_polygon.size() - 1 ; ++l)
          {
            Index lp1(l == vertex_at_polygon.size() - 1 ? 0 : l + 1);
            c_x += (origin_coords.at(0).at(vertex_at_polygon.at(l)) + origin_coords.at(0).at(vertex_at_polygon.at(lp1)))
              * (origin_coords.at(0).at(vertex_at_polygon.at(l)) * origin_coords.at(1).at(vertex_at_polygon.at(lp1)) - origin_coords.at(0).at(vertex_at_polygon.at(lp1)) * origin_coords.at(1).at(vertex_at_polygon.at(l)));
          }
          c_x *= area;

          for(Index l(0); l < vertex_at_polygon.size() - 1 ; ++l)
          {
            Index lp1(l == vertex_at_polygon.size() - 1 ? 0 : l + 1);
            c_y += (origin_coords.at(1).at(vertex_at_polygon.at(l)) + origin_coords.at(1).at(vertex_at_polygon.at(lp1)))
              * (origin_coords.at(0).at(vertex_at_polygon.at(l)) * origin_coords.at(1).at(vertex_at_polygon.at(lp1)) - origin_coords.at(0).at(vertex_at_polygon.at(lp1))*origin_coords.at(1).at(vertex_at_polygon.at(l)));
          }
          c_y *= area;

          origin_coords.at(0).push_back(c_x);
          origin_coords.at(1).push_back(c_y);
        }
      }
    };
  }
}

#endif

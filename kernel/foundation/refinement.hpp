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
            {
              target.add_adjacency(pl_vertex,
                                   pl_face,
                                   target.get_topologies().at(ipi_vertex_face).size() - 1,
                                   faces_at_end_vertex_old.at(i));

            }
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


    ///enforce interpretation of arbitrary mesh as 'residing in 3D space'
    template<typename Tag_,
             typename Arch_,
             HaloRefinementTypes hrt_>
    struct PolytopeRefinement<Tag_,
                              Arch_,
                              dim_3D,
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

        // TODO halo for 3D
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

        if(Dim_::required_num_topologies == rnt_3D)
        {
          //set face adjacencies for new vertex
          typename MeshType_<Dim_, t_, os_>::topology_type_::storage_type_ faces_at_start_vertex_old(target.get_adjacent_polytopes(pl_vertex, pl_face, vertex_at_edge_old.at(0)));
          typename MeshType_<Dim_, t_, os_>::topology_type_::storage_type_ faces_at_end_vertex_old(target.get_adjacent_polytopes(pl_vertex, pl_face, vertex_at_edge_old.at(1)));
          for(Index i(0) ; i < faces_at_end_vertex_old.size() ; ++i)
            if(std::find(faces_at_start_vertex_old.begin(), faces_at_start_vertex_old.end(), faces_at_end_vertex_old.at(i)) != faces_at_start_vertex_old.end())
            {
              target.add_adjacency(pl_vertex,
                                   pl_face,
                                   target.get_topologies().at(ipi_vertex_face).size() - 1,
                                   faces_at_end_vertex_old.at(i));

            }
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
        else if(y_1 > y_0) //must use same direction as x !!
          origin_coords.at(1).push_back(y_0 + DataType_(1./2.)*(y_1 - y_0));
        else
          origin_coords.at(1).push_back(y_1 + DataType_(1./2.)*(y_0 - y_1));

        const typename Attribute<DataType_, os_>::data_type_ z_0(origin_coords.at(2).at(vertex_at_edge.at(0)));
        const typename Attribute<DataType_, os_>::data_type_ z_1(origin_coords.at(2).at(vertex_at_edge.at(1)));

        if(std::abs(z_1 - z_0) <= std::numeric_limits<DataType_>::epsilon())
          origin_coords.at(2).push_back(z_0);
        else if(x_1 > x_0) //must use same direction as x !!
          origin_coords.at(2).push_back(z_0 + DataType_(1./2.)*(z_1 - z_0));
        else
          origin_coords.at(2).push_back(z_1 + DataType_(1./2.)*(z_0 - z_1));
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

        typename t_::storage_type_ vertex_at_polygon;
        typename t_::storage_type_ vf_0; //new vertices for new face with old index i ; need to find them in adj-list updates later

        const Index num_faces(origin.num_polytopes(pl_face));
        for(Index i(0) ; i < num_faces; ++i)
        {
          ///vertex set T*f_i <- EMPTYSET
          t_ t_fi;

          /// operator V_f->v
          typename t_::storage_type_ v_fc(coarse.get_adjacent_polytopes(pl_face, pl_vertex, i));
          vertex_at_polygon = v_fc;
          std::sort(vertex_at_polygon.begin(), vertex_at_polygon.end());

          /// operator E_v->e
          for(Index j(0) ; j < v_fc.size(); ++j)
          {
            typename t_::storage_type_ e_fo_tmp(origin.get_adjacent_polytopes(pl_vertex, pl_edge, v_fc.at(j)));

            /// operator V_e->v
            typename t_::storage_type_ v_e_v;
            for(Index k(0) ; k < e_fo_tmp.size(); ++k)
            {
              typename t_::storage_type_ v_fo_tmp(origin.get_adjacent_polytopes(pl_edge, pl_vertex, e_fo_tmp.at(k)));
              std::sort(v_fo_tmp.begin(), v_fo_tmp.end());
              typename t_::storage_type_ res(v_e_v.size() + v_fo_tmp.size());
              typename t_::storage_type_::iterator it(std::set_union(v_e_v.begin(), v_e_v.end(), v_fo_tmp.begin(), v_fo_tmp.end(), res.begin()));
              res.resize(Index(it - res.begin()));
              v_e_v = res;
            }
            t_fi.push_back(v_e_v);
          }

          /// add new vertex v_i
          origin.add_polytope(pl_vertex);

          /// operator {[T*f_i] /intersect [T*f_j]} (aka only new vertices for a new face)
          typename t_::storage_type_ t_union;
          for(Index l(0) ; l < t_fi.size() ; ++l)
          {
            ///[T*f_l] /intersect [T*f_j]
            t_ t_fi_l_intersect_t_fi_j;

            for(Index j(0) ; j < t_fi.size() ; ++j)
            {
              if(l != j)
              {
                ///T_fi_l
                typename t_::storage_type_ t_l(t_fi.at(l));
                std::sort(t_l.begin(), t_l.end());
                ///T_fi_j
                typename t_::storage_type_ t_j(t_fi.at(j));
                std::sort(t_j.begin(), t_j.end());
                ///T_fi_lj
                typename t_::storage_type_ t_lj(std::max(t_l.size(), t_j.size()));
                typename t_::storage_type_::iterator it(std::set_intersection(t_l.begin(), t_l.end(), t_j.begin(), t_j.end(), t_lj.begin()));
                t_lj.resize(Index(it - t_lj.begin()));

                if(t_lj.size() > 0)
                  t_fi_l_intersect_t_fi_j.push_back(t_lj);

                typename t_::storage_type_ union_tmp(t_union.size() + t_lj.size());
                typename t_::storage_type_::iterator itu(std::set_union(t_union.begin(), t_union.end(), t_lj.begin(), t_lj.end(), union_tmp.begin()));
                union_tmp.resize(Index(itu - union_tmp.begin()));
                t_union = union_tmp;
              }
            }

            ///here, we already have all vertices for the new face originating from old corner vertex v_fc_l aka [v_fc_l, v_i, all in t_fi_l_intersect_t_fi_j]
            ///add faces {fo_i}
            if(l == 0) //use face index i
            {
              vf_0.clear();
              vf_0.push_back(v_fc.at(l));
              vf_0.push_back(t_fi_l_intersect_t_fi_j.at(0).at(0));
              vf_0.push_back(t_fi_l_intersect_t_fi_j.at(1).at(0));
              vf_0.push_back(origin.num_polytopes(pl_vertex)- 1);

              //inject
              origin.get_topologies().at(ipi_face_vertex).at(i) = vf_0; //vertices are already adjacent to face i and  will be removed from all other later => done here
              //add new vertex adjacence to f_i
              origin.get_topologies().at(ipi_vertex_face).at(origin.num_polytopes(pl_vertex) - 1).push_back(i);
            }
            else if( l == 1)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
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
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if( l == 2)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
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
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if( l == 3)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
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
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if(l % 2 == 0)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
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
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if(l % 2 != 0)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
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
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          ///all vertices at old face
          typename t_::storage_type_ vertex_at_polygon_tmp(vertex_at_polygon.size() + t_union.size());
          typename t_::storage_type_::iterator itu(std::set_union(vertex_at_polygon.begin(), vertex_at_polygon.end(), t_union.begin(), t_union.end(), vertex_at_polygon_tmp.begin()));
          vertex_at_polygon_tmp.resize(Index(itu - vertex_at_polygon_tmp.begin()));
          vertex_at_polygon = vertex_at_polygon_tmp;

          ///add new edges {e_i_j}
          for(Index j(0) ; j < t_union.size() ; ++j)
          {
            //(., v_i)
            origin.add_polytope(pl_edge);
            origin.add_adjacency(pl_edge, pl_vertex, origin.num_polytopes(pl_edge) - 1, t_union.at(j));
            origin.add_adjacency(pl_edge, pl_vertex, origin.num_polytopes(pl_edge) - 1, origin.num_polytopes(pl_vertex) - 1);
          }

          ///remove v->f adjacencies from all vertices in vertex_at_polygon - vf_0
          std::sort(vertex_at_polygon.begin(), vertex_at_polygon.end());
          std::sort(vf_0.begin(), vf_0.end());
          typename t_::storage_type_ diff(vertex_at_polygon.size());
          typename t_::storage_type_::iterator itd(std::set_difference(vertex_at_polygon.begin(), vertex_at_polygon.end(), vf_0.begin(), vf_0.end(), diff.begin()));
          diff.resize(Index(itd - diff.begin()));

          for(Index j(0) ; j < diff.size() ; ++j)
          {
            typename t_::storage_type_::iterator ite(std::find(origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).begin(),
                                                               origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).end(),
                                                               i
                                                              ));

            if(ite != origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).end())
            {
              origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).erase(ite);
            }
          }

          //determine new vertex coordinates for mid-point
          //sorting
          vertex_at_polygon = coarse.get_adjacent_polytopes(pl_face, pl_vertex, i);
          Index v_current(*vertex_at_polygon.begin());
          Index v_start(*vertex_at_polygon.begin());
          typename t_::storage_type_ e_way;
          typename t_::storage_type_ v_way;
          do
          {
            //get all edges
            typename t_::storage_type_ e_current(coarse.get_adjacent_polytopes(pl_vertex, pl_edge, v_current));

            //restrict to those not used yet
            std::sort(e_way.begin(), e_way.end());
            std::sort(e_current.begin(), e_current.end());
            typename t_::storage_type_ e_valid(e_current.size());
            typename t_::storage_type_::iterator itv(std::set_difference(e_current.begin(), e_current.end(), e_way.begin(), e_way.end(), e_valid.begin()));
            e_valid.resize(Index(itv - e_valid.begin()));

            //search next candidate
            for(Index j(0) ; j < e_valid.size() ; ++j)
            {
              typename t_::storage_type_ v_valid_j(coarse.get_adjacent_polytopes(pl_edge, pl_vertex, e_valid.at(j)));
              Index candidate_pos(v_valid_j.at(0) == v_current ? Index(1) : Index(0));
              Index candidate(v_valid_j.at(candidate_pos));

              if(std::find(vertex_at_polygon.begin(), vertex_at_polygon.end(), candidate) != vertex_at_polygon.end())
              {
                e_way.push_back(e_valid.at(j));
                v_way.push_back(candidate);
                v_current = candidate;
                break;
              }
            }
          }
          while(v_current != v_start);
          vertex_at_polygon = v_way;

          //area
          DataType_ area(0), c_x(0), c_y(0);
          for(Index l(0); l < vertex_at_polygon.size() ; ++l)
          {
            Index lp1(l == vertex_at_polygon.size() - 1 ? 0 : l + 1);
            area += origin_coords.at(0).at(vertex_at_polygon.at(l)) * origin_coords.at(1).at(vertex_at_polygon.at(lp1)) - origin_coords.at(0).at(vertex_at_polygon.at(lp1)) * origin_coords.at(1).at(vertex_at_polygon.at(l));
          }
          area =  DataType_(1./2.) * fabs(area);

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
              * (origin_coords.at(0).at(vertex_at_polygon.at(l)) * origin_coords.at(1).at(vertex_at_polygon.at(lp1)) - origin_coords.at(0).at(vertex_at_polygon.at(lp1)) * origin_coords.at(1).at(vertex_at_polygon.at(l)));
          }
          c_y *= area;

          origin_coords.at(0).push_back(c_x);
          origin_coords.at(1).push_back(c_y);
        }
      }

      ///3D
      ///complete w. vertex data
      template<
        typename t_,
        template <typename, typename> class os_,
        template <typename, typename, template<typename, typename> class > class MeshType_,
        typename DataType_>
      static void execute(
                          MeshType_<Dim3D, t_, os_>& origin,
                          os_<std::shared_ptr<HaloBase<MeshType_<Dim3D, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<MeshType_<Dim3D, t_, os_>, os_> > > >* halos,
                          os_<Attribute<DataType_, os_>, std::allocator<Attribute<DataType_, os_> > >& origin_coords)
      {
        ///TODO assumes QUAD, needs coarse copy
        MeshType_<Dim3D, t_, os_> coarse(origin);

        //refine all edges first (sets all face/vertex adjacencies correctly)
        const Index num_edges(origin.get_topologies().at(ipi_edge_vertex).size());
        for(Index i(0) ; i < num_edges; ++i)
        {
          PolytopeRefinement<Tag_,
                             Arch_,
                             dim_3D,
                             pl_edge,
                             mrt_,
                             hrt_>::execute(
                                           i,
                                           origin,
                                           halos,
                                           origin_coords
                                           );
        }

        typename t_::storage_type_ vertex_at_polygon;
        typename t_::storage_type_ vf_0; //new vertices for new face with old index i ; need to find them in adj-list updates later

        // refine all faces
        const Index num_faces(origin.num_polytopes(pl_face));
        for(Index i(0) ; i < num_faces; i++)
        {
          ///vertex set T*f_i <- EMPTYSET
          t_ t_fi;

          /// operator V_f->v
          typename t_::storage_type_ v_fc(coarse.get_adjacent_polytopes(pl_face, pl_vertex, i));
          vertex_at_polygon = v_fc;
          std::sort(vertex_at_polygon.begin(), vertex_at_polygon.end());

          /// operator E_v->e
          for(Index j(0) ; j < v_fc.size(); ++j)
          {
            typename t_::storage_type_ e_fo_tmp(origin.get_adjacent_polytopes(pl_vertex, pl_edge, v_fc.at(j)));

            /// operator V_e->v
            typename t_::storage_type_ v_e_v;
            for(Index k(0) ; k < e_fo_tmp.size(); ++k)
            {
              typename t_::storage_type_ v_fo_tmp(origin.get_adjacent_polytopes(pl_edge, pl_vertex, e_fo_tmp.at(k)));
              std::sort(v_fo_tmp.begin(), v_fo_tmp.end());
              typename t_::storage_type_ res(v_e_v.size() + v_fo_tmp.size());
              typename t_::storage_type_::iterator it(std::set_union(v_e_v.begin(), v_e_v.end(), v_fo_tmp.begin(), v_fo_tmp.end(), res.begin()));
              res.resize(Index(it - res.begin()));
              v_e_v = res;
            }
            t_fi.push_back(v_e_v);
          }

          /// add new vertex v_i
          origin.add_polytope(pl_vertex);

          /// operator {[T*f_i] /intersect [T*f_j]} (aka only new vertices for a new face)
          typename t_::storage_type_ t_union;
          for(Index l(0) ; l < t_fi.size() ; ++l)
          {
            ///[T*f_l] /intersect [T*f_j]
            t_ t_fi_l_intersect_t_fi_j;

            for(Index j(0) ; j < t_fi.size() ; ++j)
            {
              if(l != j)
              {
                ///T_fi_l
                typename t_::storage_type_ t_l(t_fi.at(l));
                std::sort(t_l.begin(), t_l.end());
                ///T_fi_j
                typename t_::storage_type_ t_j(t_fi.at(j));
                std::sort(t_j.begin(), t_j.end());
                ///T_fi_lj
                typename t_::storage_type_ t_lj(std::max(t_l.size(), t_j.size()));
                typename t_::storage_type_::iterator it(std::set_intersection(t_l.begin(), t_l.end(), t_j.begin(), t_j.end(), t_lj.begin()));
                t_lj.resize(Index(it - t_lj.begin()));

                if(t_lj.size() > 0)
                  t_fi_l_intersect_t_fi_j.push_back(t_lj);

                typename t_::storage_type_ union_tmp(t_union.size() + t_lj.size());
                typename t_::storage_type_::iterator itu(std::set_union(t_union.begin(), t_union.end(), t_lj.begin(), t_lj.end(), union_tmp.begin()));
                union_tmp.resize(Index(itu - union_tmp.begin()));
                t_union = union_tmp;
              }
            }

            ///here, we already have all vertices for the new face originating from old corner vertex v_fc_l aka [v_fc_l, v_i, all in t_fi_l_intersect_t_fi_j]
            ///add faces {fo_i}
            if(l == 0) //use face index i
            {
              vf_0.clear();
              vf_0.push_back(v_fc.at(l));
              vf_0.push_back(t_fi_l_intersect_t_fi_j.at(0).at(0));
              vf_0.push_back(t_fi_l_intersect_t_fi_j.at(1).at(0));
              vf_0.push_back(origin.num_polytopes(pl_vertex)- 1);

              //inject
              origin.get_topologies().at(ipi_face_vertex).at(i) = vf_0; //vertices are already adjacent to face i and  will be removed from all other later => done here
              //add new vertex adjacence to f_i
              origin.get_topologies().at(ipi_vertex_face).at(origin.num_polytopes(pl_vertex) - 1).push_back(i);
            }
            else if( l == 1)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
              if(hrt_ == hrt_refine)
              {
                //refine halo faces
                if(halos != nullptr)
                {
                  for(Index hi(0) ; hi < halos->size() ; ++hi)
                  {
                    //only do this for face halos with delta > 0
                    if(halos->at(hi)->get_level() == pl_face)// && halos->at(hi)->get_overlap() > 0)
                    {
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if( l == 2)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
              if(hrt_ == hrt_refine)
              {
                //refine halo faces
                if(halos != nullptr)
                {
                  for(Index hi(0) ; hi < halos->size() ; ++hi)
                  {
                    //only do this for face halos with delta > 0
                    if(halos->at(hi)->get_level() == pl_face)// && halos->at(hi)->get_overlap() > 0)
                    {
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if( l == 3)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
              if(hrt_ == hrt_refine)
              {
                //refine halo faces
                if(halos != nullptr)
                {
                  for(Index hi(0) ; hi < halos->size() ; ++hi)
                  {
                    //only do this for face halos with delta > 0
                    if(halos->at(hi)->get_level() == pl_face)// && halos->at(hi)->get_overlap() > 0)
                    {
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if(l % 2 == 0)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
              if(hrt_ == hrt_refine)
              {
                //refine halo faces
                if(halos != nullptr)
                {
                  for(Index hi(0) ; hi < halos->size() ; ++hi)
                  {
                    //only do this for face halos with delta > 0
                    if(halos->at(hi)->get_level() == pl_face)// && halos->at(hi)->get_overlap() > 0)
                    {
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
            else if(l % 2 != 0)
            {
              origin.add_polytope(pl_face);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(0).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_fi_l_intersect_t_fi_j.at(1).at(0));
              origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, v_fc.at(l));

              //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
              if(hrt_ == hrt_refine)
              {
                //refine halo faces
                if(halos != nullptr)
                {
                  for(Index hi(0) ; hi < halos->size() ; ++hi)
                  {
                    //only do this for face halos with delta > 0
                    if(halos->at(hi)->get_level() == pl_face)// && halos->at(hi)->get_overlap() > 0)
                    {
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_face_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_face_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          ///all vertices at old face
          typename t_::storage_type_ vertex_at_polygon_tmp(vertex_at_polygon.size() + t_union.size());
          std::sort(vertex_at_polygon.begin(), vertex_at_polygon.end());
          typename t_::storage_type_::iterator itu(std::set_union(vertex_at_polygon.begin(), vertex_at_polygon.end(), t_union.begin(), t_union.end(), vertex_at_polygon_tmp.begin()));
          vertex_at_polygon_tmp.resize(Index(itu - vertex_at_polygon_tmp.begin()));
          vertex_at_polygon = vertex_at_polygon_tmp;

          ///add new edges {e_i_j}
          for(Index j(0) ; j < t_union.size() ; ++j)
          {
            //(., v_i)
            origin.add_polytope(pl_edge);
            origin.add_adjacency(pl_edge, pl_vertex, origin.num_polytopes(pl_edge) - 1, t_union.at(j));
            origin.add_adjacency(pl_edge, pl_vertex, origin.num_polytopes(pl_edge) - 1, origin.num_polytopes(pl_vertex) - 1);
          }

          ///remove v->f adjacencies from all vertices in vertex_at_polygon - vf_0
          std::sort(vertex_at_polygon.begin(), vertex_at_polygon.end());
          std::sort(vf_0.begin(), vf_0.end());
          typename t_::storage_type_ diff(vertex_at_polygon.size());
          typename t_::storage_type_::iterator itd(std::set_difference(vertex_at_polygon.begin(), vertex_at_polygon.end(), vf_0.begin(), vf_0.end(), diff.begin()));
          diff.resize(Index(itd - diff.begin()));

          for(Index j(0) ; j < diff.size() ; ++j)
          {
            typename t_::storage_type_::iterator ite(std::find(origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).begin(),
                                                               origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).end(),
                                                               i
                                                              ));

            if(ite != origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).end())
            {
              origin.get_topologies().at(ipi_vertex_face).at(diff.at(j)).erase(ite);
            }
          }

          //determine new vertex coordinates for mid-point
          //sorting
          vertex_at_polygon = coarse.get_adjacent_polytopes(pl_face, pl_vertex, i);
          Index v_current(*vertex_at_polygon.begin());
          Index v_start(*vertex_at_polygon.begin());
          typename t_::storage_type_ e_way;
          typename t_::storage_type_ v_way;
          do
          {
            //get all edges
            typename t_::storage_type_ e_current(coarse.get_adjacent_polytopes(pl_vertex, pl_edge, v_current));

            //restrict to those not used yet
            std::sort(e_way.begin(), e_way.end());
            std::sort(e_current.begin(), e_current.end());
            typename t_::storage_type_ e_valid(e_current.size());
            typename t_::storage_type_::iterator itv(std::set_difference(e_current.begin(), e_current.end(), e_way.begin(), e_way.end(), e_valid.begin()));
            e_valid.resize(Index(itv - e_valid.begin()));

            //search next candidate
            for(Index j(0) ; j < e_valid.size() ; ++j)
            {
              typename t_::storage_type_ v_valid_j(coarse.get_adjacent_polytopes(pl_edge, pl_vertex, e_valid.at(j)));
              Index candidate_pos(v_valid_j.at(0) == v_current ? Index(1) : Index(0));
              Index candidate(v_valid_j.at(candidate_pos));

              if(std::find(vertex_at_polygon.begin(), vertex_at_polygon.end(), candidate) != vertex_at_polygon.end())
              {
                e_way.push_back(e_valid.at(j));
                v_way.push_back(candidate);
                v_current = candidate;
                break;
              }
            }
          }
          while(v_current != v_start);
          vertex_at_polygon = v_way;

          //area, center, midpoint
          DataType_ area(0), c_x(0), c_y(0), c_z(0), m_x(0), m_y(0), m_z(0);

          //calculate a midpoint
          for(Index l(0) ; l < vertex_at_polygon.size() ; ++l)
          {
            m_x += origin_coords.at(0).at(vertex_at_polygon.at(l));
            m_y += origin_coords.at(1).at(vertex_at_polygon.at(l));
            m_z += origin_coords.at(2).at(vertex_at_polygon.at(l));
          }
          m_x /= double(vertex_at_polygon.size());
          m_y /= double(vertex_at_polygon.size());
          m_z /= double(vertex_at_polygon.size());

          //Subdivide polygon into triangles
          for(Index l(0) ; l < vertex_at_polygon.size() ; ++l)
          {
            double a, b, c, s, areal; // edge lengths a,b,c of triangle, perimeter/2 s
            Index lp1(l == vertex_at_polygon.size() - 1 ? 0 : l + 1);
            a = sqrt( (origin_coords.at(0).at(vertex_at_polygon.at(l)) - origin_coords.at(0).at(vertex_at_polygon.at(lp1)))
              * (origin_coords.at(0).at(vertex_at_polygon.at(l)) - origin_coords.at(0).at(vertex_at_polygon.at(lp1))) + (origin_coords.at(1).at(vertex_at_polygon.at(l)) - origin_coords.at(1).at(vertex_at_polygon.at(lp1)))
              * (origin_coords.at(1).at(vertex_at_polygon.at(l)) - origin_coords.at(1).at(vertex_at_polygon.at(lp1))) + (origin_coords.at(2).at(vertex_at_polygon.at(l)) - origin_coords.at(2).at(vertex_at_polygon.at(lp1)))
              * (origin_coords.at(2).at(vertex_at_polygon.at(l)) - origin_coords.at(2).at(vertex_at_polygon.at(lp1))));
            b = sqrt( (origin_coords.at(0).at(vertex_at_polygon.at(l)) - m_x)
              * (origin_coords.at(0).at(vertex_at_polygon.at(l)) - m_x) + (origin_coords.at(1).at(vertex_at_polygon.at(l)) - m_y)
              * (origin_coords.at(1).at(vertex_at_polygon.at(l)) - m_y) + (origin_coords.at(2).at(vertex_at_polygon.at(l)) - m_z)
              * (origin_coords.at(2).at(vertex_at_polygon.at(l)) - m_z) );
            c = sqrt( (origin_coords.at(0).at(vertex_at_polygon.at(lp1)) - m_x)
              * (origin_coords.at(0).at(vertex_at_polygon.at(lp1)) - m_x) + (origin_coords.at(1).at(vertex_at_polygon.at(lp1)) - m_y)
              * (origin_coords.at(1).at(vertex_at_polygon.at(lp1)) - m_y) + (origin_coords.at(2).at(vertex_at_polygon.at(lp1)) - m_z)
              * (origin_coords.at(2).at(vertex_at_polygon.at(lp1)) - m_z) );
            s = 0.5 * (a + b + c);
            areal = sqrt( s * (s - a) * (s - b) * (s - c) );

            c_x += (origin_coords.at(0).at(vertex_at_polygon.at(l)) + origin_coords.at(0).at(vertex_at_polygon.at(lp1)) + m_x) * areal;
            c_y += (origin_coords.at(1).at(vertex_at_polygon.at(l)) + origin_coords.at(1).at(vertex_at_polygon.at(lp1)) + m_y) * areal;
            c_z += (origin_coords.at(2).at(vertex_at_polygon.at(l)) + origin_coords.at(2).at(vertex_at_polygon.at(lp1)) + m_z) * areal;
            area += areal;
          }
          c_x /= 3.0 * area;
          c_y /= 3.0 * area;
          c_z /= 3.0 * area;

          origin_coords.at(0).push_back(c_x);
          origin_coords.at(1).push_back(c_y);
          origin_coords.at(2).push_back(c_z);

        }


        // here, we already refine all edges and faces. now, we need the center of each polyhedron and the inner edges, faces to add the new polyhedron
        typename t_::storage_type_ vertex_at_polyhedron;
        typename t_::storage_type_ vp_0; //new vertices for new face with old index i ; need to find them in adj-list updates later

        const Index num_polyhedron(origin.num_polytopes(pl_polyhedron));
        for(Index i(0) ; i < num_polyhedron; i++)
        {
          // add new vertex v_i
          origin.add_polytope(pl_vertex);

          ///vertex set T*p_i <- EMPTYSET
          t_ t_pi;

          /// operator V_p->v
          typename t_::storage_type_ v_pc(coarse.get_adjacent_polytopes(pl_polyhedron, pl_vertex, i));
          vertex_at_polyhedron = v_pc;
          std::sort(vertex_at_polyhedron.begin(), vertex_at_polyhedron.end());

          /// operator E_v->e
          for(Index j(0) ; j < v_pc.size(); ++j)
          {
            typename t_::storage_type_ e_po_tmp(origin.get_adjacent_polytopes(pl_vertex, pl_edge, v_pc.at(j)));

            /// operator V_e->v
            typename t_::storage_type_ v_e_v;
            for(Index k(0) ; k < e_po_tmp.size(); ++k)
            {
              typename t_::storage_type_ v_po_tmp(origin.get_adjacent_polytopes(pl_edge, pl_vertex, e_po_tmp.at(k)));
              std::sort(v_po_tmp.begin(), v_po_tmp.end());
              typename t_::storage_type_ res(v_e_v.size() + v_po_tmp.size());
              typename t_::storage_type_::iterator it(std::set_union(v_e_v.begin(), v_e_v.end(), v_po_tmp.begin(), v_po_tmp.end(), res.begin()));
              res.resize(Index(it - res.begin()));
              v_e_v = res;
            }
            t_pi.push_back(v_e_v);
          }

          /// operator {[T*p_i] /intersect [T*p_j]} (aka only new vertices without center)
          typename t_::storage_type_ t_union;
          for(Index l(0) ; l < t_pi.size() ; ++l)
          {
            ///[T*p_l] /intersect [T*p_j]
            t_ t_pi_l_intersect_t_pi_j;

            for(Index j(0) ; j < t_pi.size() ; ++j)
            {
              if(l != j)
              {
                ///T_pi_l
                typename t_::storage_type_ t_l(t_pi.at(l));
                std::sort(t_l.begin(), t_l.end());
                ///T_pi_j
                typename t_::storage_type_ t_j(t_pi.at(j));
                std::sort(t_j.begin(), t_j.end());
                ///T_pi_lj
                typename t_::storage_type_ t_lj(std::max(t_l.size(), t_j.size()));
                typename t_::storage_type_::iterator it(std::set_intersection(t_l.begin(), t_l.end(), t_j.begin(), t_j.end(), t_lj.begin()));
                t_lj.resize(Index(it - t_lj.begin()));

                if(t_lj.size() > 0)
                  t_pi_l_intersect_t_pi_j.push_back(t_lj);

                typename t_::storage_type_ union_tmp(t_union.size() + t_lj.size());
                typename t_::storage_type_::iterator itu(std::set_union(t_union.begin(), t_union.end(), t_lj.begin(), t_lj.end(), union_tmp.begin()));
                union_tmp.resize(Index(itu - union_tmp.begin()));
                t_union = union_tmp;
              }
            }
          }

          ///t_center: center of coarse faces
          t_ t_unioni;

          /// operator E_v->e
          for(Index j(0) ; j < t_union.size(); ++j)
          {
            typename t_::storage_type_ e_union_tmp(origin.get_adjacent_polytopes(pl_vertex, pl_edge, t_union.at(j)));

            /// operator V_e->v
            typename t_::storage_type_ v_e_union;
            for(Index k(0) ; k < e_union_tmp.size(); ++k)
            {
              typename t_::storage_type_ v_union_tmp(origin.get_adjacent_polytopes(pl_edge, pl_vertex, e_union_tmp.at(k)));
              std::sort(v_union_tmp.begin(), v_union_tmp.end());
              typename t_::storage_type_ res(v_e_union.size() + v_union_tmp.size());
              typename t_::storage_type_::iterator it(std::set_union(v_e_union.begin(), v_e_union.end(), v_union_tmp.begin(), v_union_tmp.end(), res.begin()));
              res.resize(Index(it - res.begin()));
              v_e_union = res;
            }
            t_unioni.push_back(v_e_union);
          }

          /// operator {[T*p_i] /intersect [T*p_j]} (aka only new vertices for a new face)
          typename t_::storage_type_ t_center;
          for(Index l(0) ; l < t_unioni.size() ; ++l)
          {
            ///[T*p_l] /intersect [T*p_j]
            t_ t_unioni_l_intersect_t_unioni_j;

            for(Index j(0) ; j < t_unioni.size() ; ++j)
            {
              if(l != j)
              {
                ///T_pi_l
                typename t_::storage_type_ t_l(t_unioni.at(l));
                std::sort(t_l.begin(), t_l.end());
                ///T_pi_j
                typename t_::storage_type_ t_j(t_unioni.at(j));
                std::sort(t_j.begin(), t_j.end());
                ///T_pi_lj
                typename t_::storage_type_ t_lj(std::max(t_l.size(), t_j.size()));
                typename t_::storage_type_::iterator it(std::set_intersection(t_l.begin(), t_l.end(), t_j.begin(), t_j.end(), t_lj.begin()));
                t_lj.resize(Index(it - t_lj.begin()));

                if(t_lj.size() > 0)
                  t_unioni_l_intersect_t_unioni_j.push_back(t_lj);

                typename t_::storage_type_ union_tmp(t_center.size() + t_lj.size());
                typename t_::storage_type_::iterator itu(std::set_union(t_center.begin(), t_center.end(), t_lj.begin(), t_lj.end(), union_tmp.begin()));
                union_tmp.resize(Index(itu - union_tmp.begin()));
                t_center = union_tmp;
              }
            }
          }

          // first, add new inner edges
          typename t_::storage_type_ union_tmp(t_center.size());
          std::sort(vertex_at_polyhedron.begin(), vertex_at_polyhedron.end());
          typename t_::storage_type_::iterator itu(std::set_difference(t_center.begin(), t_center.end(), vertex_at_polyhedron.begin(), vertex_at_polyhedron.end(), union_tmp.begin()));
          union_tmp.resize(Index(itu - union_tmp.begin()));
          t_center = union_tmp;


          for(Index j(0) ; j < t_center.size() ; ++j)
          {
            origin.add_polytope(pl_edge);
            origin.add_adjacency(pl_vertex, pl_edge, t_center.at(j), origin.num_polytopes(pl_edge) - 1);
            origin.add_adjacency(pl_vertex, pl_edge, origin.num_polytopes(pl_vertex) - 1, origin.num_polytopes(pl_edge) - 1);
          }

          // here, we already have all points and all edges
          // now, add new inner faces
          for(Index j(0) ; j < t_center.size() ; ++j)
          {
            // V_E_vj
            typename t_::storage_type_ e_center_tmpj(origin.get_adjacent_polytopes(pl_vertex, pl_edge, t_center.at(j)));

            typename t_::storage_type_ v_e_vj;
            for(Index k(0) ; k < e_center_tmpj.size(); ++k)
            {
              typename t_::storage_type_ v_center_tmp(origin.get_adjacent_polytopes(pl_edge, pl_vertex, e_center_tmpj.at(k)));
              std::sort(v_center_tmp.begin(), v_center_tmp.end());
              typename t_::storage_type_ res(v_e_vj.size() + v_center_tmp.size());
              typename t_::storage_type_::iterator it(std::set_union(v_e_vj.begin(), v_e_vj.end(), v_center_tmp.begin(), v_center_tmp.end(), res.begin()));
              res.resize(Index(it - res.begin()));
              v_e_vj = res;
            }

            for(Index l(j + 1) ; l < t_center.size() ; ++l)// add new inner faces
            {
              // V_E_vl
              typename t_::storage_type_ e_center_tmpl(origin.get_adjacent_polytopes(pl_vertex, pl_edge, t_center.at(l)));

              typename t_::storage_type_ v_e_vl;
              for(Index k(0) ; k < e_center_tmpl.size(); ++k)
              {
                typename t_::storage_type_ v_center_tmp(origin.get_adjacent_polytopes(pl_edge, pl_vertex, e_center_tmpl.at(k)));
                std::sort(v_center_tmp.begin(), v_center_tmp.end());
                typename t_::storage_type_ res(v_e_vl.size() + v_center_tmp.size());
                typename t_::storage_type_::iterator it(std::set_union(v_e_vl.begin(), v_e_vl.end(), v_center_tmp.begin(), v_center_tmp.end(), res.begin()));
                res.resize(Index(it - res.begin()));
                v_e_vl = res;
              }

              typename t_::storage_type_ vjl(v_e_vl.size());
              typename t_::storage_type_ vjl_tmp(v_e_vl.size());
              typename t_::storage_type_::iterator it(std::set_intersection(v_e_vl.begin(), v_e_vl.end(), v_e_vj.begin(), v_e_vj.end(), vjl_tmp.begin()));
              vjl_tmp.resize(Index(it - vjl_tmp.begin()));
              typename t_::storage_type_ v_center_tmp;
              v_center_tmp.push_back(origin.num_polytopes(pl_vertex) - 1);
              typename t_::storage_type_::iterator itu1(std::set_difference(vjl_tmp.begin(), vjl_tmp.end(), v_center_tmp.begin(), v_center_tmp.end(), vjl.begin()));
              vjl.resize(Index(itu1 - vjl.begin()));

              if(vjl.size() > 0)
              {
                origin.add_polytope(pl_face);
                origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, origin.num_polytopes(pl_vertex) - 1);
                origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_center.at(j));
                origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, t_center.at(l));
                origin.add_adjacency(pl_face, pl_vertex, origin.num_polytopes(pl_face) - 1, vjl.at(0) );
              }

            }
          }

          // here, we already have all faces. now, add new polyhedron
          /// operator F_v->f
          typename t_::storage_type_ f_vcenter(origin.get_adjacent_polytopes(pl_vertex,pl_face,origin.num_polytopes(pl_vertex) - 1));

          for(Index j(0) ; j < vertex_at_polyhedron.size() ; ++j)
          {
            typename t_::storage_type_ f_pnew;

            /// operator F_v->f
            typename t_::storage_type_ f_vj(origin.get_adjacent_polytopes(pl_vertex,pl_face,vertex_at_polyhedron.at(j)));

            for(Index k(0) ; k < f_vj.size() ; ++k)
            {
              for(Index l(0); l < f_vcenter.size() ; ++l)
              {
                /// operator E_f->e
                typename t_::storage_type_ e_fvc(origin.get_adjacent_polytopes(pl_face, pl_edge, f_vcenter.at(l)));
                typename t_::storage_type_ e_fvj(origin.get_adjacent_polytopes(pl_face, pl_edge, f_vj.at(k)));

                typename t_::storage_type_ intersect_cj(e_fvc.size());
                std::sort(e_fvc.begin(), e_fvc.end());
                std::sort(e_fvj.begin(), e_fvj.end());
                typename t_::storage_type_::iterator it(std::set_intersection(e_fvc.begin(), e_fvc.end(), e_fvj.begin(), e_fvj.end(), intersect_cj.begin()));
                intersect_cj.resize(Index(it - intersect_cj.begin()));

                if( intersect_cj.size() > 0)
                {
                  typename t_::storage_type_ f_pnew_add;
                  f_pnew_add.push_back(f_vj.at(k));
                  f_pnew_add.push_back(f_vcenter.at(l));
                  std::sort(f_pnew_add.begin(),f_pnew_add.end());

                  typename t_::storage_type_ f_pnew_tmp(f_pnew.size()+f_pnew_add.size());
                  typename t_::storage_type_::iterator itu2(std::set_union(f_pnew_add.begin(), f_pnew_add.end(), f_pnew.begin(), f_pnew.end(), f_pnew_tmp.begin()));
                  f_pnew_tmp.resize(Index(itu2 - f_pnew_tmp.begin()));

                  f_pnew = f_pnew_tmp;
                }
              }
            }

            //
            typename t_::storage_type_ v_pnew;
            for(Index k(0) ; k < f_pnew.size() ; ++k)
            {
              typename t_::storage_type_ v_fk(origin.get_adjacent_polytopes(pl_face,pl_vertex,f_pnew.at(k)));
              typename t_::storage_type_ v_pnew_tmp(v_pnew.size()+v_fk.size());
              std::sort(v_fk.begin(), v_fk.end());
              typename t_::storage_type_::iterator it(std::set_union(v_pnew.begin(), v_pnew.end(), v_fk.begin(), v_fk.end(), v_pnew_tmp.begin()));
              v_pnew_tmp.resize(Index(it - v_pnew_tmp.begin()));
              v_pnew = v_pnew_tmp;
            }

            // add polyhedron p_oi
            if(j == 0) // use polyhedron index i
            {
              vp_0.clear();
              for(Index k(0) ; k < v_pnew.size() ; ++k)
              {
                vp_0.push_back(v_pnew.at(k));
              }
              //inject
              origin.get_topologies().at(ipi_polyhedron_vertex).at(i) = vp_0; //vertices are already adjacent to polyhedron i and  will be removed from all other later => done here
            }
            else
            {
              origin.add_polytope(pl_polyhedron);
              for(Index k(0) ; k < v_pnew.size() ; ++k)
              {
                origin.add_adjacency(pl_vertex, pl_polyhedron, v_pnew.at(k), origin.num_polytopes(pl_polyhedron) - 1);
              }

             //for each refined face, insert the new face numbers into halos with delta > 0 AFTER the old face (covers Halo<k>, k > 0)
              if(hrt_ == hrt_refine)
              {
                //refine halo faces
                if(halos != nullptr)
                {
                  for(Index hi(0) ; hi < halos->size() ; ++hi)
                  {
                    //only do this for face halos with delta > 0
                    if(halos->at(hi)->get_level() == pl_polyhedron && halos->at(hi)->get_overlap() > 0)
                    {
                      typename t_::storage_type_::iterator iter(std::find(halos->at(hi)->get_elements().begin(), halos->at(hi)->get_elements().end(), i));
                      if(iter != halos->at(hi)->get_elements().end())
                      {
                        if(iter + 1 == halos->at(hi)->get_elements().end())
                        {
                          halos->at(hi)->get_elements().push_back(origin.get_topologies().at(ipi_polyhedron_vertex).size() - 1);
                        }
                        else
                        {
                          halos->at(hi)->get_elements().insert(iter + 1, origin.get_topologies().at(ipi_polyhedron_vertex).size() - 1); //insertion after ///TODO QUAD
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          ///remove v->p adjacencies from all vertices in vertex_at_polygon - vp_0
          std::sort(vertex_at_polyhedron.begin(), vertex_at_polyhedron.end());
          for(Index j(0) ; j < vertex_at_polyhedron.size() ; ++j)
          {
            typename t_::storage_type_::iterator ite(std::find(origin.get_topologies().at(ipi_vertex_polyhedron).at(vertex_at_polyhedron.at(j)).begin(),
                                                               origin.get_topologies().at(ipi_vertex_polyhedron).at(vertex_at_polyhedron.at(j)).end(),
                                                               i
                                                              ));

            if(ite != origin.get_topologies().at(ipi_vertex_polyhedron).at(vertex_at_polyhedron.at(j)).end())
            {
              origin.get_topologies().at(ipi_vertex_polyhedron).at(vertex_at_polyhedron.at(j)).erase(ite);
            }
          }

          for(Index j(0); j < origin.get_adjacent_polytopes(pl_polyhedron,pl_vertex,i).size(); ++j)
          {
            //add new vertex adjacence to f_i
            origin.get_topologies().at(ipi_vertex_polyhedron).at(origin.get_adjacent_polytopes(pl_polyhedron,pl_vertex,i).at(j)).push_back(i);
          }

          DataType_ vol(0), c_x(0), c_y(0), c_z(0);
          Index num_faces1((Index)coarse.get_adjacent_polytopes(pl_polyhedron, pl_face, i).size());
          for(Index j(0) ; j < num_faces1; j++)
          {
            //sorting
            vertex_at_polygon = coarse.get_adjacent_polytopes(pl_face, pl_vertex, j);
            Index v_current(*vertex_at_polygon.begin());
            Index v_start(*vertex_at_polygon.begin());
            typename t_::storage_type_ e_way;
            typename t_::storage_type_ v_way;
            do
            {
              //get all edges
              typename t_::storage_type_ e_current(coarse.get_adjacent_polytopes(pl_vertex, pl_edge, v_current));

              //restrict to those not used yet
              std::sort(e_way.begin(), e_way.end());
              std::sort(e_current.begin(), e_current.end());
              typename t_::storage_type_ e_valid(e_current.size());
              typename t_::storage_type_::iterator itv(std::set_difference(e_current.begin(), e_current.end(), e_way.begin(), e_way.end(), e_valid.begin()));
              e_valid.resize(Index(itv - e_valid.begin()));

              //search next candidate
              for(Index k(0) ; k < e_valid.size() ; ++k)
              {
                typename t_::storage_type_ v_valid_k(coarse.get_adjacent_polytopes(pl_edge, pl_vertex, e_valid.at(k)));
                Index candidate_pos(v_valid_k.at(0) == v_current ? Index(1) : Index(0));
                Index candidate(v_valid_k.at(candidate_pos));

                if(std::find(vertex_at_polygon.begin(), vertex_at_polygon.end(), candidate) != vertex_at_polygon.end())
                {
                  e_way.push_back(e_valid.at(k));
                  v_way.push_back(candidate);
                  v_current = candidate;
                  break;
                }
              }
            }
            while(v_current != v_start);
            vertex_at_polygon = v_way;

            // Calculating centroid of polyhedron
            Index num_vertices((Index)coarse.get_adjacent_polytopes(pl_face, pl_vertex, j).size());
            for(Index k(0) ; k < num_vertices ; k++)
            {
              double n[3];
              DataType_ areak(0);
              n[0] = (origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices)) - origin_coords.at(1).at(vertex_at_polygon.at(k)))
                   * (origin_coords.at(2).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) - origin_coords.at(2).at(vertex_at_polygon.at((k + 1) % num_vertices))) - (origin_coords.at(1).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) - origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                   * (origin_coords.at(2).at(vertex_at_polygon.at((k + 1) % num_vertices)) - origin_coords.at(2).at(vertex_at_polygon.at(k)));
              n[1] = (origin_coords.at(2).at(vertex_at_polygon.at((k + 1) % num_vertices)) - origin_coords.at(2).at(vertex_at_polygon.at(k)))
                   * (origin_coords.at(0).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) - origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices))) - (origin_coords.at(2).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) - origin_coords.at(2).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                   * (origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices)) - origin_coords.at(0).at(vertex_at_polygon.at(k)));
              n[2] = (origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices)) - origin_coords.at(0).at(vertex_at_polygon.at(k)))
                   * (origin_coords.at(1).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) - origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices))) - (origin_coords.at(0).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) - origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                   * (origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices)) - origin_coords.at(1).at(vertex_at_polygon.at(k)));

              areak= origin_coords.at(0).at(vertex_at_polygon.at(k)) * n[0] + origin_coords.at(1).at(vertex_at_polygon.at(k)) * n[1] + origin_coords.at(2).at(vertex_at_polygon.at(k)) * n[2];
              if (areak < 0) //outer normal
                for(Index m(0); m < 3 ;m++) n[m] *= DataType_(-1);
              vol += fabs(areak);

              c_x += n[0] * (((origin_coords.at(0).at(vertex_at_polygon.at(k)) + origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                            * (origin_coords.at(0).at(vertex_at_polygon.at(k)) + origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices))))
                           + ((origin_coords.at(0).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) + origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                            * (origin_coords.at(0).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) + origin_coords.at(0).at(vertex_at_polygon.at((k + 1) % num_vertices))))
                           + ((origin_coords.at(0).at(vertex_at_polygon.at(k)) + origin_coords.at(0).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)))
                            * (origin_coords.at(0).at(vertex_at_polygon.at(k)) + origin_coords.at(0).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)))) );

              c_y += n[1] * (((origin_coords.at(1).at(vertex_at_polygon.at(k)) + origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                            * (origin_coords.at(1).at(vertex_at_polygon.at(k)) + origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices))))
                           + ((origin_coords.at(1).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) + origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                            * (origin_coords.at(1).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) + origin_coords.at(1).at(vertex_at_polygon.at((k + 1) % num_vertices))))
                           + ((origin_coords.at(1).at(vertex_at_polygon.at(k)) + origin_coords.at(1).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)))
                            * (origin_coords.at(1).at(vertex_at_polygon.at(k)) + origin_coords.at(1).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)))) );

              c_z += n[2] * (((origin_coords.at(2).at(vertex_at_polygon.at(k)) + origin_coords.at(2).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                            * (origin_coords.at(2).at(vertex_at_polygon.at(k)) + origin_coords.at(2).at(vertex_at_polygon.at((k + 1) % num_vertices))))
                           + ((origin_coords.at(2).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) + origin_coords.at(2).at(vertex_at_polygon.at((k + 1) % num_vertices)))
                            * (origin_coords.at(2).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)) + origin_coords.at(2).at(vertex_at_polygon.at(( k + 1) % num_vertices))))
                           + ((origin_coords.at(2).at(vertex_at_polygon.at(k)) + origin_coords.at(2).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)))
                            * (origin_coords.at(2).at(vertex_at_polygon.at(k)) + origin_coords.at(2).at(origin.num_polytopes(pl_vertex) - (num_faces1 - j + 1)))) );
             }

           }
           origin_coords.at(0).push_back(c_x / (DataType_(8) * vol));
           origin_coords.at(1).push_back(c_y / (DataType_(8) * vol));
           origin_coords.at(2).push_back(c_z / (DataType_(8) * vol));
        }
      }
    };
  }
}

#endif

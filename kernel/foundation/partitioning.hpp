#pragma once
#ifndef FOUNDATION_GUARD_PARTITIONING_HPP
#define FOUNDATION_GUARD_PARTITIONING_HPP 1

#include<kernel/foundation/data.hpp>
#include<kernel/foundation/refinement.hpp>
#include<kernel/foundation/sub_mesh.hpp>
#include<algorithm>

using namespace FEAST;

namespace FEAST
{
  namespace Foundation
  {
    template<typename Dim_,
             MeshRefinementTypes mrt_>
    struct _OnPartitioning
    {
    };

    template<typename Dim_>
    struct _OnPartitioning<Dim_, mrt_standard>
    {
      static inline Index num_refinements_needed(Index num_cells, Index num_procs)
      {
        return Index(std::ceil(std::log(double(num_procs)/double(num_cells))/std::log(std::pow(double(2), double(Dim_::tag_value)))));
      }

      template<typename t_,
               template <typename, typename> class os_,
               template <typename, typename, template<typename, typename> class > class MeshType_
               >
      static inline void increase_partition(const MeshType_<Dim_, t_, os_>& m,
                                            t_& target,
                                            Index rank,
                                            Index num_des_elems_local,
                                            Index& num_dist_elems_global,
                                            Index& num_dist_elems_local,
                                            Index& last)
      {
        //recursion end
        if(num_dist_elems_global == m.num_polytopes(Dim_::ElementPolytopeType_::tag_value) || num_dist_elems_local == num_des_elems_local)
          return;

        typename t_::storage_type_ nb(m.get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::tag_value, last));
        for(Index j(0) ; j < nb.size() ; ++j) //for every element adjacent to last
        {
          //add neighbours found
          bool found(false);
          for(Index k(0) ; k < target.size() ; ++k)
            if(std::find(target.at(k).begin() , target.at(k).end(), nb.at(j)) < target.at(k).end())
            {
              found = true;
              break;
            }

          if(!found)
          {
            if(!(num_dist_elems_global == m.num_polytopes(Dim_::ElementPolytopeType_::tag_value) || num_dist_elems_local == num_des_elems_local))
            {
              target.at(rank).push_back(nb.at(j));
              ++num_dist_elems_global;
              ++num_dist_elems_local;
            }

            if(nb.at(j) != last)
            {
              last = nb.at(j);
              increase_partition(m,
                  target,
                  rank,
                  num_des_elems_local,
                  num_dist_elems_global,
                  num_dist_elems_local,
                  last);
            }
          }
        }
      }
    };

    template<typename Tag_,
             typename Arch_,
             typename Dim_,
             unsigned delta_,
             PolytopeLevels level_,
             MeshRefinementTypes mrt_ = mrt_standard,
             HaloRefinementTypes hrt_ = hrt_refine>
    struct Partitioning
    {
      template<
        typename t_,
        template <typename, typename> class os_,
        template <typename, typename, template<typename, typename> class > class MeshType_,
        typename DT_>
      static PData<Dim_, t_, os_, MeshType_, DT_> execute(MeshType_<Dim_, t_, os_>& mesh,
                                                          os_<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, Mesh<Dim_, t_, os_>, os_>,
                                                              std::allocator<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, Mesh<Dim_, t_, os_>, os_> > >& boundaries,
                                                          const Index num_procs,
                                                          const Index proc_rank,
                                                          os_<Attribute<DT_, os_>, std::allocator<Attribute<DT_, os_> > >& origin_coords)
      {
        PData<Dim_, t_, os_, MeshType_, DT_> result;

        ///-->Each process for all
        ///Phase 1: refine to match process count
        //only case: more procs than polytopes: refine num_refinements_needed times

        if(num_procs > mesh.get_topologies().at(mesh.get_topologies().size() - 1).size())
        {

          Index num_refs(_OnPartitioning<Dim_, mrt_>::num_refinements_needed(mesh.get_topologies().at(mesh.get_topologies().size() - 1).size(), num_procs));

          ///TODO opt:
          os_<std::shared_ptr<HaloBase<Mesh<Dim_, t_, os_>, os_> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim_, t_, os_>, os_> > > > boundaries_copy;
          for(Index i(0) ; i < boundaries.size() ; ++i)
            boundaries_copy.push_back(std::shared_ptr<HaloBase<Mesh<Dim_, t_, os_>, os_> >(new Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, Mesh<Dim_, t_, os_>, os_>(boundaries.at(i))));

          for(Index i(0) ; i < num_refs ; ++i)
          {
            Refinement<Tag_, Arch_, mrt_, hrt_refine>::execute(mesh,
                                                               &boundaries_copy,
                                                               origin_coords);
          }

          boundaries.clear();
          for(Index i(0) ; i < boundaries_copy.size() ; ++i)
            boundaries.push_back(*((Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, Mesh<Dim_, t_, os_>, os_>*)(boundaries_copy.at(i).get())));
        }

        ///Phase 2: initial geometric partitioning
        t_ elements_for_process(num_procs);
        if(num_procs < mesh.num_polytopes(Dim_::ElementPolytopeType_::tag_value))
        {
          Index num_elems_distributed_local(0);
          Index num_elems_distributed_global(0);
          Index num_elems_global(mesh.num_polytopes(Dim_::ElementPolytopeType_::tag_value));
          os_<Index, std::allocator<Index> > num_desired_elems_local(num_procs, Index(0));

          Index elems_dist(0), iter(0);
          while(elems_dist < num_elems_global)
          {
            ++elems_dist;
            ++num_desired_elems_local.at(iter);

            iter = iter < num_procs - 1 ? iter + 1 : 0;
          }

          Index last(0);

          for(Index i(0) ; i < num_procs ; ++i) //for every process
          {
            num_elems_distributed_local = 0;
            _OnPartitioning<Dim_, mrt_>::increase_partition(mesh,
                elements_for_process,
                i,
                num_desired_elems_local.at(i),
                num_elems_distributed_global,
                num_elems_distributed_local,
                last);
          }
        }
        else if(num_procs == mesh.num_polytopes(Dim_:: ElementPolytopeType_::tag_value))
        {
          for(Index i(0) ; i < elements_for_process.size() ; ++i)
            elements_for_process.at(i).push_back(i);
        }

        //create commhalos (global view)
        std::vector<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_> > comm_halos_globalview;
        std::vector<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_> > comm_halos_globalview_sub;
        std::vector<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_> > comm_halos_globalview_sub_sub;

        typename t_::storage_type_ comm_halo_originating_from_element;
        typename t_::storage_type_ comm_halo_originating_from_element_sub;
        typename t_::storage_type_ comm_halo_originating_from_element_sub_sub;
        typename t_::storage_type_ comm_halo_originating_from_process;
        typename t_::storage_type_ comm_halo_originating_from_process_sub;
        typename t_::storage_type_ comm_halo_originating_from_process_sub_sub;

        for(Index i(0) ; i < elements_for_process.size() ; ++i)
          for(Index j(0) ; j < elements_for_process.at(i).size() ; ++j)
            for(Index k(0) ; k < elements_for_process.size() ; ++k)
              for(Index l(0) ; l < elements_for_process.at(k).size() ; ++l)
              {
                if(i != k)
                {
                  typename t_::storage_type_ tmp(mesh.get_comm_intersection(Dim_::ElementPolytopeType_::tag_value,
                                                                            Dim_::ElementPolytopeType_::SubElementPolytopeType_::tag_value,
                                                                            elements_for_process.at(i).at(j),
                                                                            elements_for_process.at(k).at(l)));
                  typename t_::storage_type_ tmp_sub(mesh.get_comm_intersection(Dim_::ElementPolytopeType_::tag_value,
                                                                                Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::tag_value,
                                                                                elements_for_process.at(i).at(j),
                                                                                elements_for_process.at(k).at(l)));
                  typename t_::storage_type_ tmp_sub_sub(mesh.get_comm_intersection(Dim_::ElementPolytopeType_::tag_value,
                                                                                    Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::tag_value,
                                                                                    elements_for_process.at(i).at(j),
                                                                                    elements_for_process.at(k).at(l)));
                  if(tmp.size() > 0)
                  {
                    comm_halos_globalview.push_back(Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_>(mesh, k));
                    comm_halo_originating_from_element.push_back(elements_for_process.at(i).at(j));
                    comm_halo_originating_from_process.push_back(i);

                    for(Index m(0) ; m < tmp.size() ; ++m)
                      comm_halos_globalview.at(comm_halos_globalview.size() - 1).push_back(tmp.at(m));
                  }

                  //sublevel 1
                  if(tmp.size() == 0)
                  {
                    if(tmp_sub.size() > 0)
                    {
                      comm_halos_globalview_sub.push_back(Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_>(mesh, k));
                      comm_halo_originating_from_element_sub.push_back(elements_for_process.at(i).at(j));
                      comm_halo_originating_from_process_sub.push_back(i);

                      for(Index m(0) ; m < tmp_sub.size() ; ++m)
                        comm_halos_globalview_sub.at(comm_halos_globalview_sub.size() - 1).push_back(tmp_sub.at(m));
                    }
                  }

                  //sublevel 2
                  if(tmp_sub.size() == 0 && tmp.size() == 0)
                  {
                    if(tmp_sub_sub.size() > 0)
                    {
                      comm_halos_globalview_sub_sub.push_back(Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_>(mesh, k));
                      comm_halo_originating_from_element_sub_sub.push_back(elements_for_process.at(i).at(j));
                      comm_halo_originating_from_process_sub_sub.push_back(i);

                      for(Index m(0) ; m < tmp_sub_sub.size() ; ++m)
                        comm_halos_globalview_sub_sub.at(comm_halos_globalview_sub_sub.size() - 1).push_back(tmp_sub_sub.at(m));
                    }
                  }
                }
              }

        ///-->Each process for its patch
        ///submesh and halo generation; here not delta_ may be used but rather delta=0, delta_ applies later for the micro-mesh
        //create submesh
        Halo<0, typename Dim_::ElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_ > submeshproxy(mesh);
        for(unsigned long i(0) ; i < elements_for_process.at(proc_rank).size() ; ++i)
          submeshproxy.push_back(elements_for_process.at(proc_rank).at(i));

        result.submesh = std::shared_ptr<SubMesh<Dim_, t_, os_> >(new SubMesh<Dim_, t_, os_>(&submeshproxy));


        //restrict boundary components belonging to this patch to submesh
        std::vector<Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_ > > local_boundaries;

        const Index i_end(Index(boundaries.size()));
        for(unsigned long i(0) ; i < i_end; ++i)
        {
          local_boundaries.push_back(Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_ >(mesh));

          const Index j_end(boundaries.at(i).size());
          for(unsigned long j(0) ; j <  j_end ; ++j)
          {
            Index boundary_element(boundaries.at(i).get_element(j));

            const Index k_end(result.submesh->num_polytopes(Dim_::ElementPolytopeType_::tag_value));
            for(unsigned long k(0) ; k < k_end ; ++k)
            {
              typename t_::storage_type_ adj_elements_origin(mesh.get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::tag_value, result.submesh->get_map().at(k)));
              typename t_::storage_type_ adj_elements_submesh(result.submesh->get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::tag_value, k));

              typename t_::storage_type_::iterator found(std::find(adj_elements_origin.begin(), adj_elements_origin.end(), boundary_element));

              if(found != adj_elements_origin.end())
              {
                Index found_diff(Index(found - adj_elements_origin.begin()));
                local_boundaries.at(i).push_back(adj_elements_submesh.at(found_diff));
              }
            }
          }
        }

        //kick out boundary-components not belonging to process
        for(Index i(0) ; i < local_boundaries.size() ; ++i)
          if(local_boundaries.at(i).size() > 0)
          {
            result.boundaries.push_back(Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_>(local_boundaries.at(i)));
            result.boundaries.at(result.boundaries.size() - 1).reset_mesh((MeshType_<Dim_, t_, os_>*)&result.submesh); ///TODO conceptual: only result.comm_halos of Mesh (and not SubMesh) needed?
          }

        //create commhalos (local view)
        for(Index i(0) ; i < comm_halos_globalview.size() ; ++i)
        {
          if(comm_halo_originating_from_process.at(i) == proc_rank)
          {
            typename t_::storage_type_ tmp_globalview(mesh.get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::tag_value, comm_halo_originating_from_element.at(i)));

            const typename t_::storage_type_& map(result.submesh->get_map());
            typename t_::storage_type_::const_iterator f(std::find(map.begin(),
                                                         map.end(),
                                                         comm_halo_originating_from_element.at(i)));
            Index element_in_submesh(Index(f - map.begin()));

            typename t_::storage_type_ tmp_localview(result.submesh->get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::tag_value, element_in_submesh));

            Index halo_position_in_adj_list(Index(std::find(tmp_globalview.begin(), tmp_globalview.end(), comm_halos_globalview.at(i).get_element(0)) - tmp_globalview.begin()));
            //create new halo
            result.comm_halos.push_back(std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> >(new Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_ >()));
            result.comm_halos.at(result.comm_halos.size() -1)->reset_mesh((MeshType_<Dim_, t_, os_>*)&result.submesh);
            result.comm_halos.at(result.comm_halos.size() -1)->reset_other(comm_halos_globalview.at(i).get_other());
            result.comm_halos.at(result.comm_halos.size() -1)->push_back(tmp_localview.at(halo_position_in_adj_list));
          }
        }

        for(Index i(0) ; i < comm_halos_globalview_sub.size() ; ++i)
        {
          if(comm_halo_originating_from_process_sub.at(i) == proc_rank)
          {
            typename t_::storage_type_ tmp_globalview_sub(mesh.get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::tag_value, comm_halo_originating_from_element_sub.at(i)));

            const typename t_::storage_type_& map(result.submesh->get_map());
            typename t_::storage_type_::const_iterator f(std::find(map.begin(),
                                                         map.end(),
                                                         comm_halo_originating_from_element_sub.at(i)));
            Index element_in_submesh(Index(f - map.begin()));

            typename t_::storage_type_ tmp_localview(result.submesh->get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::tag_value, element_in_submesh));

            Index halo_position_in_adj_list(Index(std::find(tmp_globalview_sub.begin(), tmp_globalview_sub.end(), comm_halos_globalview_sub.at(i).get_element(0)) - tmp_globalview_sub.begin()));
            //create new halo
            result.comm_halos.push_back(std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> >(new Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_ >()));
            result.comm_halos.at(result.comm_halos.size() -1)->reset_mesh((MeshType_<Dim_, t_, os_>*)&result.submesh);
            result.comm_halos.at(result.comm_halos.size() -1)->reset_other(comm_halos_globalview_sub.at(i).get_other());
            result.comm_halos.at(result.comm_halos.size() -1)->push_back(tmp_localview.at(halo_position_in_adj_list));
          }
        }

        for(Index i(0) ; i < comm_halos_globalview_sub_sub.size() ; ++i)
        {
          if(comm_halo_originating_from_process_sub_sub.at(i) == proc_rank)
          {
            typename t_::storage_type_ tmp_globalview_sub_sub(mesh.get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::tag_value, comm_halo_originating_from_element_sub_sub.at(i)));

            const typename t_::storage_type_& map(result.submesh->get_map());
            typename t_::storage_type_::const_iterator f(std::find(map.begin(),
                                                         map.end(),
                                                         comm_halo_originating_from_element_sub_sub.at(i)));
            Index element_in_submesh(Index(f - map.begin()));

            typename t_::storage_type_ tmp_localview(result.submesh->get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::tag_value, element_in_submesh));

            Index halo_position_in_adj_list(Index(std::find(tmp_globalview_sub_sub.begin(), tmp_globalview_sub_sub.end(), comm_halos_globalview_sub_sub.at(i).get_element(0)) - tmp_globalview_sub_sub.begin()));
            //create new halo
            result.comm_halos.push_back(std::shared_ptr<HaloBase<MeshType_<Dim_, t_, os_>, os_> >(new Halo<0, typename Dim_::ElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_::SubElementPolytopeType_, MeshType_<Dim_, t_, os_>, os_ >()));
            result.comm_halos.at(result.comm_halos.size() -1)->reset_mesh((MeshType_<Dim_, t_, os_>*)&result.submesh);
            result.comm_halos.at(result.comm_halos.size() -1)->reset_other(comm_halos_globalview_sub_sub.at(i).get_other());
            result.comm_halos.at(result.comm_halos.size() -1)->push_back(tmp_localview.at(halo_position_in_adj_list));
          }
        }
        std::sort(result.comm_halos.begin(), result.comm_halos.end(), compare_other<MeshType_<Dim_, t_, os_>, os_>);

        //restrict x_coords attr to submesh
        result.attrs.push_back(std::shared_ptr<AttributeBase<os_> >(new Attribute<DT_, os_>())); //local parts, empty at first
        if(Dim_::tag_value > dim_1D)
          result.attrs.push_back(std::shared_ptr<AttributeBase<os_> >(new Attribute<DT_, os_>())); //local parts, empty at first
        if(Dim_::tag_value > dim_2D)
          result.attrs.push_back(std::shared_ptr<AttributeBase<os_> >(new Attribute<DT_, os_>())); //local parts, empty at first

        for(Index i(0) ; i < result.submesh->num_polytopes(pl_vertex) ; ++i)
        {
            //local elements adjacent to vertex i
            typename t_::storage_type_ elements_i_local(result.submesh->get_adjacent_polytopes(pl_vertex, Dim_::ElementPolytopeType_::tag_value, i));

            //which position has vertex i in adj list of specific element elements_i_local.at(0)?
            typename t_::storage_type_ tmp_local(result.submesh->get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, pl_vertex, elements_i_local.at(0)));
            Index pos(Index(std::find(tmp_local.begin(), tmp_local.end(), i) - tmp_local.begin()));

            //now get element in corresponding supermesh
            Index element_global(result.submesh->get_map().at(elements_i_local.at(0)));

            //get adjacent vertices for that element in supermesh
            typename t_::storage_type_ tmp_global(mesh.get_adjacent_polytopes(Dim_::ElementPolytopeType_::tag_value, pl_vertex, element_global));

            //access adj list right at position pos
            Index target_vertex_index(tmp_global.at(pos));

            //push back corresponding global value into local mesh
            ((Attribute<DT_, os_>*)(result.attrs.at(0).get()))->push_back(origin_coords.at(0).at(target_vertex_index));
            if(Dim_::tag_value > dim_1D)
            {
              ((Attribute<DT_, os_>*)(result.attrs.at(1).get()))->push_back(origin_coords.at(1).at(target_vertex_index));
            }

            if(Dim_::tag_value > dim_2D)
            {
              ((Attribute<DT_, os_>*)(result.attrs.at(2).get()))->push_back(origin_coords.at(2).at(target_vertex_index));
            }
        }

        return result;
      }
    };
  }
}

#endif

#pragma once
#ifndef KERNEL_FOUNDATION_HALO_CONTROL_HPP
#define KERNEL_FOUNDATION_HALO_CONTROL_HPP

#include <kernel/foundation/halo.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template<Dimensions dim_, typename DT_ = double>
    struct HaloControl
    {
    };

    template<typename DT_>
    struct HaloControl<dim_1D, DT_>
    {
        ///delta = 0 case: in 1D, zero-overlap halos can only be given in terms of vertices
        ///0, pl_vertex case
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, PLVertex, b_, DT_, c_>& halo, typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_* target)
        {
          if (halo.size() != 1)
            throw InternalError("Error: Halo with 0-overlap may not contain more than one vertex in 1D!");
          typedef typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_ IndexType;
          target[0] = IndexType(1);
          target[1] = IndexType(0);
        }

        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLVertex, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<1> >& target)
        {
          ASSERT(halo.size() == 1, "Error: Halo with 0-overlap may not contain more than one vertex in 1D!");
          target.template get_target_set<0>()[0] = halo.get_element(0);
        }

        //reverse
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<1> >& source, HaloType_<0, PLVertex, b_, DT_, c_>& halo)
        {
          halo.get_elements().clear();
          halo.push_back(source.template get_target_set<0>()[0]);
        }

        ///TODO remove these and move them to their actual dimension!
        ///overload for Hypercube<2> (in order to get point-diagonal halos)
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLVertex, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<2> >& target)
        {
          ASSERT(halo.size() == 1, "Error: Halo with 0-overlap may not contain more than one vertex in 1D!");
          target.template get_target_set<0>()[0] = halo.get_element(0);
        }

        //reverse
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<2> >& source, HaloType_<0, PLVertex, b_, DT_, c_>& halo)
        {
          halo.get_elements().clear();
          halo.push_back(source.template get_target_set<0>()[0]);
        }

        ///overload for Hypercube<3>
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLEdge, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<3> >& target)
        {
          ASSERT(halo.size() == 1, "Error: Halo with 0-overlap may not contain more than one vertex in 1D!");
          target.template get_target_set<1>()[0] = halo.get_element(0);
        }

        //reverse
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<3> >& source, HaloType_<0, PLVertex, b_, DT_, c_>& halo)
        {
          halo.get_elements().clear();
          halo.push_back(source.template get_target_set<0>()[0]);
        }

        ///delta = i case: in 1D, overlapping meshes have halos given in terms of edges
        ///i, pl_edge case
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<a_, PLEdge, b_, DT_, c_>& halo, typename HaloType_<a_, PLEdge, b_, DT_, c_>::index_type_* target)
        {
          ASSERT(a_ != 0, "Error: Halos with 0-overlap may not contain edges in 1D!");

          typedef typename HaloType_<a_, PLEdge, b_, DT_, c_>::index_type_ IndexType;
          target[0] = IndexType(halo.size() + 1);
          target[1] = IndexType(halo.size());
        }

        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<a_, PLEdge, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<1> >& target)
        {
          ASSERT(a_ != 0, "Error: Halos with 0-overlap may not contain edges in 1D!");

          typedef typename HaloType_<a_, PLEdge, b_, DT_, c_>::index_type_ IndexType;

          typename HaloType_<a_, PLEdge, b_, DT_, c_>::mesh_type_::storage_type_ vertices_performed;

          IndexType vertex_index(0);

          ///for any edge add vertices
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            target.template get_target_set<1>()[i] = halo.get_element(i);

            typename HaloType_<a_, PLEdge, b_, DT_, c_>::mesh_type_::storage_type_ adjacent_vertices(halo.get_mesh()->get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, halo.get_element(i)));
            for(IndexType j(0) ; j < adjacent_vertices.size() ; ++j)
            {
              if(std::find(vertices_performed.begin(), vertices_performed.end(), adjacent_vertices.at(j)) == vertices_performed.end()) //not yet performed
              {
                target.template get_target_set<0>()[vertex_index] = adjacent_vertices.at(j);
                vertices_performed.push_back(adjacent_vertices.at(j));
                ++vertex_index;
              }
            }

          }
        }

        //reverse
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<1> >& source, HaloType_<a_, PLEdge, b_, DT_, c_>& halo)
        {
          typedef typename HaloType_<a_, PLEdge, b_, DT_, c_>::index_type_ IndexType;

          Index size(source.template get_target_set<1>().get_num_entities());

          for(IndexType i(0) ; i < size ; ++i)
            halo.push_back(source.template get_target_set<1>()[i]);
        }

    };

    template<typename DT_>
    struct HaloControl<dim_2D, DT_>
    {
      public:

        ///delta = 0 case: in 2D, zero-overlap halos can be given in terms of edges
        ///0, pl_edge case
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, PLEdge, b_, DT_, c_>& halo, typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_* target)
        {
          typedef typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_ IndexType;
          target[0] = IndexType(halo.size()  + 1);
          target[1] = IndexType(halo.size());
          target[2] = IndexType(0);
        }

        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLEdge, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<2> >& target)
        {
          typedef typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_ IndexType;
          typename HaloType_<0, PLEdge, b_, DT_, c_>::mesh_type_::storage_type_ vertices_performed;
          IndexType vertex_index(0);

          ///for any edge add vertices
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            typename HaloType_<0, PLEdge, b_, DT_, c_>::mesh_type_::storage_type_ adjacent_vertices(halo.get_mesh()->get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, halo.get_element(i)));

            target.template get_target_set<1>()[i] = halo.get_element(i);

            for(IndexType j(0) ; j < adjacent_vertices.size() ; ++j)
            {
              if(std::find(vertices_performed.begin(), vertices_performed.end(), adjacent_vertices.at(j)) == vertices_performed.end())
              {
                target.template get_target_set<0>()[vertex_index] = adjacent_vertices.at(j);
                vertices_performed.push_back(adjacent_vertices.at(j));
                ++vertex_index;
              }
            }
          }
        }

        //reverse
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<2> >& source, HaloType_<0, PLEdge, b_, DT_, c_>& halo)
        {
          typedef typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_ IndexType;

          Index size(source.template get_target_set<1>().get_num_entities());

          for(IndexType i(0) ; i < size ; ++i)
            halo.push_back(source.template get_target_set<1>()[i]);
        }

        ///overload for Hypercube<3>
        /*template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLEdge, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<3> >& target)
        {
          typedef typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_ IndexType;

          ///for any edge add vertices
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            typename HaloType_<0, PLEdge, b_, DT_, c_>::mesh_type_::storage_type_ adjacent_vertices(halo.get_mesh()->get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, halo.get_element(i)));

            target.template get_target_set<1>()[i] = halo.get_element(i);
            for(IndexType j(0) ; j < adjacent_vertices.size() ; ++j)
            {
              target.template get_target_set<0>()[j] = adjacent_vertices.at(j);
            }
          }
        }

        //reverse
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<3> >& source, HaloType_<0, PLEdge, b_, DT_, c_>& halo)
        {
          typedef typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_ IndexType;

          IndexType size(source.template get_target_set<1>().get_num_entities());

          for(IndexType i(0) ; i < size ; ++i)
            halo.push_back(source.template get_target_set<1>()[i]);
        }*/

        ///delta = 0 case: in 2D, zero-overlap halos can be given in terms of vertices (diagonal case)
        ///0, pl_vertex case
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, PLVertex, b_, DT_, c_>& halo,
                               typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_* target)
        {
          typedef typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_ IndexType;
          target[0] = IndexType(halo.size());
          target[1] = IndexType(0);
          target[2] = IndexType(0);
        }

        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLVertex, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<2> >& target)
        {
          typedef typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_ IndexType;
          typename HaloType_<0, PLVertex, b_, DT_, c_>::mesh_type_::storage_type_ vertices_performed;
          IndexType vertex_index(0);

          ///add all vertices
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            if(std::find(vertices_performed.begin(), vertices_performed.end(), halo.get_element(i)) == vertices_performed.end())
            {
              target.template get_target_set<0>()[vertex_index] = halo.get_element(i);
              vertices_performed.push_back(halo.get_element(i));
              ++vertex_index;
            }
          }
        }


        ///delta = i case: in 2D, delta > 0 halos must be given in terms of faces
        ///i, pl_face case
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<a_, PLFace, b_, DT_, c_>& halo,
                               typename HaloType_<a_, PLFace, b_, DT_, c_>::index_type_* target)
        {
          ASSERT(a_ != 0, "Error: Halos with 0-overlap may not contain faces in 2D!");

          typedef typename HaloType_<a_, PLFace, b_, DT_, c_>::index_type_ IndexType;

          IndexType num_edges(0);
          IndexType num_vertices(0);
          typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh()->get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)));
            IndexType count(0);
            for(IndexType j(0) ; j < edges.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_edges.size() ; ++k)
              {
                if(all_edges.at(k) == edges.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_edges.push_back(edges.at(j));
                ++count;
              }
            }
            num_edges += count;

          }

          ///for any edge computed before, check vertex and add if not yet added
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh()->get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
            IndexType count(0);
            for(IndexType j(0) ; j < vertices.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_vertices.size() ; ++k)
              {
                if(all_vertices.at(k) == vertices.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_vertices.push_back(vertices.at(j));
                ++count;
              }
            }
            num_vertices += count;
          }

          target[0] = IndexType(num_vertices);
          target[1] = IndexType(num_edges);
          target[2] = IndexType(halo.size());
        }

        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<a_, PLFace, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<2> >& target)
        {
          ASSERT(a_ != 0, "Error: Halos with 0-overlap may not contain faces in 2D!");

          typedef typename HaloType_<a_, PLFace, b_, DT_, c_>::index_type_ IndexType;

          typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh()->get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)));
            for(IndexType j(0) ; j < edges.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_edges.size() ; ++k)
              {
                if(all_edges.at(k) == edges.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_edges.push_back(edges.at(j));
              }
            }
          }

          ///for any edge computed before, check vertex and add if not yet added
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<a_, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh()->get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
            for(IndexType j(0) ; j < vertices.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_vertices.size() ; ++k)
              {
                if(all_vertices.at(k) == vertices.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_vertices.push_back(vertices.at(j));
              }
            }
          }

          ///transfer precomputed sets
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            target.template get_target_set<2>()[i] = halo.get_element(i);
          }
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            target.template get_target_set<1>()[i] = all_edges.at(i);
          }
          for(IndexType i(0) ; i < all_vertices.size() ; ++i)
          {
            target.template get_target_set<0>()[i] = all_vertices.at(i);
          }
        }

        //reverse
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<2> >& source, HaloType_<a_, PLFace, b_, DT_, c_>& halo)
        {
          typedef typename HaloType_<a_, PLFace, b_, DT_, c_>::index_type_ IndexType;

          IndexType size(source.template get_target_set<2>().get_num_entities());

          for(IndexType i(0) ; i < size ; ++i)
            halo.push_back(source.template get_target_set<2>()[i]);
        }
    };

    template<typename DT_>
    struct HaloControl<dim_3D, DT_>
    {
        ///delta = 0 case: in 3D, zero-overlap halos can only be given disambigously by faces
        ///0, pl_face case
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, PLFace, b_, DT_, c_>& halo,
                              typename HaloType_<0, PLFace, b_, DT_, c_>::index_type_* target)
        {
          typedef typename HaloType_<0, PLFace, b_, DT_, c_>::index_type_ IndexType;

          IndexType num_edges(0);
          IndexType num_vertices(0);
          typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh()->get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)));
            IndexType count(0);
            for(IndexType j(0) ; j < edges.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_edges.size() ; ++k)
              {
                if(all_edges.at(k) == edges.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_edges.push_back(edges.at(j));
                ++count;
              }
            }
            num_edges += count;

          }

          ///for any edge computed before, check vertex and add if not yet added
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh()->get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
            IndexType count(0);
            for(IndexType j(0) ; j < vertices.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_vertices.size() ; ++k)
              {
                if(all_vertices.at(k) == vertices.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_vertices.push_back(vertices.at(j));
                ++count;
              }
            }
            num_vertices += count;
          }

          target[0] = IndexType(num_vertices);
          target[1] = IndexType(num_edges);
          target[2] = IndexType(halo.size());
          target[3] = IndexType(0);
        }

        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLFace, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<3> >& target)
        {
          typedef typename HaloType_<0, PLFace, b_, DT_, c_>::index_type_ IndexType;

          typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh()->get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)));
            for(IndexType j(0) ; j < edges.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_edges.size() ; ++k)
              {
                if(all_edges.at(k) == edges.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_edges.push_back(edges.at(j));
              }
            }
          }

          ///for any edge computed before, check vertex and add if not yet added
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<0, PLFace, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh()->get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
            for(IndexType j(0) ; j < vertices.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_vertices.size() ; ++k)
              {
                if(all_vertices.at(k) == vertices.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_vertices.push_back(vertices.at(j));
              }
            }
          }

          ///transfer precomputed sets
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            target.template get_target_set<2>()[i] = halo.get_element(i);
          }
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            target.template get_target_set<1>()[i] = all_edges.at(i);
          }
          for(IndexType i(0) ; i < all_vertices.size() ; ++i)
          {
            target.template get_target_set<0>()[i] = all_vertices.at(i);
          }
        }

        //reverse
        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<3> >& source, HaloType_<0, PLFace, b_, DT_, c_>& halo)
        {
          typedef typename HaloType_<0, PLFace, b_, DT_, c_>::index_type_ IndexType;

          IndexType size(source.template get_target_set<2>().get_num_entities());

          for(Index i(0) ; i < size ; ++i)
            halo.push_back(source.template get_target_set<2>()[i]);
        }

        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, PLVertex, b_, DT_, c_>& halo,
                              typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_* target)
        {
          typedef typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_ IndexType;

          target[0] = IndexType(halo.size());
          target[1] = IndexType(0);
          target[2] = IndexType(0);
          target[3] = IndexType(0);
        }

        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLVertex, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<3> >& target)
        {
          typedef typename HaloType_<0, PLVertex, b_, DT_, c_>::index_type_ IndexType;
          ///add all vertices
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
              target.template get_target_set<0>()[i] = halo.get_element(i);
          }

        }


        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, PLEdge, b_, DT_, c_>& halo,
                              typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_* target)
        {
          typedef typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_ IndexType;
          target[0] = IndexType(halo.size()  + 1);
          target[1] = IndexType(halo.size());
          target[2] = IndexType(0);
          target[3] = IndexType(0);
        }

        template<
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<0, PLEdge, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<3> >& target)
        {
          typedef typename HaloType_<0, PLEdge, b_, DT_, c_>::index_type_ IndexType;

          ///for any edge add vertices
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            typename HaloType_<0, PLEdge, b_, DT_, c_>::mesh_type_::storage_type_ adjacent_vertices(halo.get_mesh()->get_adjacent_polytopes(Foundation::pl_edge, Foundation::pl_vertex, halo.get_element(i)));

            target.template get_target_set<1>()[i] = halo.get_element(i);
            for(IndexType j(0) ; j < adjacent_vertices.size() ; ++j)
            {
              target.template get_target_set<0>()[j] = adjacent_vertices.at(j);
            }
          }

        }

        ///delta = i case: in 3D, halos with overlap = i can only be given by polyhedrons
        ///i, pl_polyhedron case
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_sizes(const HaloType_<a_, PLPolyhedron, b_, DT_, c_>& halo,
                               typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::index_type_* target)
        {
          typedef typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::index_type_ IndexType;

          IndexType num_faces(0);
          IndexType num_edges(0);
          IndexType num_vertices(0);
          typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_faces;
          typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ faces(halo.get_mesh()->get_adjacent_polytopes(pl_polyhedron, pl_face, halo.get_element(i)));
            IndexType count(0);
            for(IndexType j(0) ; j < faces.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_faces.size() ; ++k)
              {
                if(all_faces.at(k) == faces.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_faces.push_back(faces.at(j));
                ++count;
              }
            }
            num_faces += count;
          }

          ///for any face computed before, check edge and add if not yet added
          for(IndexType i(0) ; i < all_faces.size() ; ++i)
          {
            typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh()->get_adjacent_polytopes(pl_face, pl_edge, all_faces.at(i)));
            IndexType count(0);
            for(IndexType j(0) ; j < edges.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_edges.size() ; ++k)
              {
                if(all_edges.at(k) == edges.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_edges.push_back(edges.at(j));
                ++count;
              }
            }
            num_edges += count;
          }

          ///for any edge computed before, check vertex and add if not yet added
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh()->get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
            IndexType count(0);
            for(IndexType j(0) ; j < vertices.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_vertices.size() ; ++k)
              {
                if(all_vertices.at(k) == vertices.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_vertices.push_back(vertices.at(j));
                ++count;
              }
            }
            num_vertices += count;
          }

          target[0] = IndexType(num_vertices);
          target[1] = IndexType(num_edges);
          target[2] = IndexType(num_faces);
          target[3] = IndexType(halo.size());
        }

        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const HaloType_<a_, PLPolyhedron, b_, DT_, c_>& halo, Geometry::CellSubSet<Shape::Hypercube<3> >& target)
        {
          typedef typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::index_type_ IndexType;

          typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_faces;
          typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ faces(halo.get_mesh()->get_adjacent_polytopes(pl_polyhedron, pl_face, halo.get_element(i)));
            for(IndexType j(0) ; j < faces.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_faces.size() ; ++k)
              {
                if(all_faces.at(k) == faces.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_faces.push_back(faces.at(j));
              }
            }
          }

          ///for any face computed before, check edge and add if not yet added
          for(IndexType i(0) ; i < all_faces.size() ; ++i)
          {
            typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh()->get_adjacent_polytopes(pl_face, pl_edge, all_faces.at(i)));
            for(IndexType j(0) ; j < edges.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_edges.size() ; ++k)
              {
                if(all_edges.at(k) == edges.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_edges.push_back(edges.at(j));
              }
            }
          }

          ///for any edge computed before, check vertex and add if not yet added
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh()->get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
            for(IndexType j(0) ; j < vertices.size() ; ++j)
            {
              bool already_in(false);
              for(IndexType k(0) ; k < all_vertices.size() ; ++k)
              {
                if(all_vertices.at(k) == vertices.at(j))
                  already_in = true;
              }
              if(!already_in)
              {
                all_vertices.push_back(vertices.at(j));
              }
            }
          }

          ///transfer precomputed sets
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            target.template get_target_set<3>()[i] = halo.get_element(i);
          }
          for(IndexType i(0) ; i < all_faces.size() ; ++i)
          {
            target.template get_target_set<2>()[i] = all_faces.at(i);
          }
          for(IndexType i(0) ; i < all_edges.size() ; ++i)
          {
            target.template get_target_set<1>()[i] = all_edges.at(i);
          }
          for(IndexType i(0) ; i < all_vertices.size() ; ++i)
          {
            target.template get_target_set<0>()[i] = all_vertices.at(i);
          }
        }

        //reverse
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          template<unsigned,
            typename,
            typename,
            typename,
            template<typename, typename> class>
          class HaloType_>
        static void fill_target_set(const Geometry::CellSubSet<Shape::Hypercube<3> >& source, HaloType_<a_, PLPolyhedron, b_, DT_, c_>& halo)
        {
          typedef typename HaloType_<a_, PLPolyhedron, b_, DT_, c_>::index_type_ IndexType;

          IndexType size(source.template get_target_set<3>().get_num_entities());

          for(IndexType i(0) ; i < size ; ++i)
            halo.push_back(source.template get_target_set<3>()[i]);
        }
    };
  }
}

#endif

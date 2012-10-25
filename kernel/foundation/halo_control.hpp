#pragma once
#ifndef KERNEL_FOUNDATION_HALO_CONTROL_HPP
#define KERNEL_FOUNDATION_HALO_CONTROL_HPP

#include <kernel/foundation/halo.hpp>
#include<kernel/geometry/conformal_mesh.hpp>

using namespace FEAST::Geometry;

namespace FEAST
{
  namespace Foundation
  {
    template<Dimensions dim_>
    struct HaloControl
    {
    };

    template<>
    struct HaloControl<dim_1D>
    {
        ///delta = 0 case: in 1D, zero-overlap halos can only be given in terms of vertices
        ///0, pl_vertex case
        template<
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, pl_vertex, b_, c_, d_>& halo, typename HaloType_<0, pl_vertex, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<0, pl_vertex, b_, c_, d_>::index_type_ IndexType;
          target[0] = IndexType(halo.size());
          target[1] = IndexType(halo.size() - 1);
        }

        ///delta = i case: in 1D, overlapping meshes have halos given in terms of edges
        ///i, pl_edge case
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<a_, pl_edge, b_, c_, d_>& halo, typename HaloType_<a_, pl_edge, b_, c_, d_>::index_type_* target)
        {
          ASSERT(a_ != 0, "Error: Halos with 0-overlap may not contain edges in 1D!");

          typedef typename HaloType_<0, pl_edge, b_, c_, d_>::index_type_ IndexType;
          target[0] = IndexType(halo.size() + 1);
          target[1] = IndexType(halo.size());
        }
    };

    template<>
    struct HaloControl<dim_2D>
    {
      public:

        ///delta = 0 case: in 2D, zero-overlap halos can be given in terms of edges OR in terms of vertices
        ///0, pl_vertex case
        template<
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, pl_vertex, b_, c_, d_>& halo, typename HaloType_<0, pl_vertex, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<0, pl_vertex, b_, c_, d_>::index_type_ IndexType;
          target[0] = IndexType(halo.size());
          target[1] = IndexType(halo.size() - 1);
          target[2] = IndexType(0);
        }

        ///0, pl_edge case
        template<
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, pl_edge, b_, c_, d_>& halo, typename HaloType_<0, pl_edge, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<0, pl_edge, b_, c_, d_>::index_type_ IndexType;
          target[0] = IndexType(halo.size()  + 1);
          target[1] = IndexType(halo.size());
          target[2] = IndexType(0);
        }

        ///delta = i case: in 2D, delta > 0 halos must be given in terms of faces
        ///i, pl_face case
        template<
          unsigned a_,
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<a_, pl_face, b_, c_, d_>& halo, typename HaloType_<a_, pl_face, b_, c_, d_>::index_type_* target)
        {
          ASSERT(a_ != 0, "Error: Halos with 0-overlap may not contain faces in 2D!");

          typedef typename HaloType_<a_, pl_face, b_, c_, d_>::index_type_ IndexType;

          IndexType num_edges(0);
          IndexType num_vertices(0);
          typename HaloType_<0, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<0, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<a_, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh().get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)));
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
            typename HaloType_<0, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh().get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
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
    };

    template<>
    struct HaloControl<dim_3D>
    {
        ///delta = 0 case: in 3D, zero-overlap halos can only be given disambigously by faces
        ///0, pl_face case
        template<
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, pl_face, b_, c_, d_>& halo, typename HaloType_<0, pl_face, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<0, pl_face, b_, c_, d_>::index_type_ IndexType;

          IndexType num_edges(0);
          IndexType num_vertices(0);
          typename HaloType_<0, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<0, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            typename HaloType_<0, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh().get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)));
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
            typename HaloType_<0, pl_face, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh().get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
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

        ///delta = i case: in 3D, halos with overlap = i can only be given by polyhedrons
        ///0, pl_polyhedron case
        template<
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<0, pl_polyhedron, b_, c_, d_>& halo, typename HaloType_<0, pl_polyhedron, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<0, pl_polyhedron, b_, c_, d_>::index_type_ IndexType;

          IndexType num_faces(0);
          IndexType num_edges(0);
          IndexType num_vertices(0);
          typename HaloType_<0, pl_polyhedron, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ all_faces;
          typename HaloType_<0, pl_polyhedron, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ all_edges;
          typename HaloType_<0, pl_polyhedron, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ all_vertices;
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            typename HaloType_<0, pl_polyhedron, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ faces(halo.get_mesh().get_adjacent_polytopes(pl_polyhedron, pl_face, halo.get_element(i)));
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
            typename HaloType_<0, pl_polyhedron, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ edges(halo.get_mesh().get_adjacent_polytopes(pl_face, pl_edge, all_faces.at(i)));
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
            typename HaloType_<0, pl_polyhedron, b_, c_, d_>::mesh_type_::topology_type_::storage_type_ vertices(halo.get_mesh().get_adjacent_polytopes(pl_edge, pl_vertex, all_edges.at(i)));
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
    };

  }
}

#endif

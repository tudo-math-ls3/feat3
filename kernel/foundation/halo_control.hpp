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
          typedef typename HaloType_<a_, pl_edge, b_, c_, d_>::index_type_ IndexType;
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
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            num_edges += halo.get_mesh().get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)).size();
            num_vertices += halo.get_mesh().get_adjacent_polytopes(pl_face, pl_vertex, halo.get_element(i)).size();
          }

          target[0] = IndexType(num_vertices);
          target[1] = IndexType(num_edges);
          target[2] = IndexType(halo.size());
        }
    };

    template<>
    struct HaloControl<dim_3D>
    {
      ///TODO
    };

  }
}

#endif

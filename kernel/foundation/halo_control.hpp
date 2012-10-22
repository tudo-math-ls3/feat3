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
    };

    template<>
    struct HaloControl<dim_2D>
    {
      public:

        /*template<
          unsigned delta_,
          PolytopeLevels a_,
          typename b_,
          template<typename, typename> class c_,
          typename d_,
          template<unsigned,
            PolytopeLevels,
            typename,
            template<typename, typename> class,
            typename>
          class HaloType_>
        static void fill_sizes(const HaloType_<delta_, a_, b_, c_, d_>& halo, typename HaloType_<delta_, a_, b_, c_, d_>::index_type_* target)
        {
        }*/

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
        static void fill_sizes(const HaloType_<0, pl_vertex, b_, c_, d_>& halo, typename HaloType_<delta_, a_, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<delta_, a_, b_, c_, d_>::index_type_ IndexType;
          target[0] = IndexType(halo.size());
          target[1] = IndexType(halo.size() / 2);
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
        static void fill_sizes(const HaloType_<0, pl_edge, b_, c_, d_>& halo, typename HaloType_<delta_, a_, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<delta_, a_, b_, c_, d_>::index_type_ IndexType;
          target[0] = IndexType(halo.size() * 2);
          target[1] = IndexType(halo.size());
          target[2] = IndexType(0);
        }

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
        static void fill_sizes(const HaloType_<0, pl_face, b_, c_, d_>& halo, typename HaloType_<delta_, a_, b_, c_, d_>::index_type_* target)
        {
          typedef typename HaloType_<delta_, a_, b_, c_, d_>::index_type_ IndexType;

          IndexType num_edges(0);
          for(IndexType i(0) ; i < halo.size() ; ++i)
          {
            ///for any face count edges
            num_edges += halo.get_mesh().get_adjacent_polytopes(pl_face, pl_edge, halo.get_element(i)).size();
          }

          target[0] = IndexType(num_edges * 2);
          target[1] = IndexType(num_edges);
          target[2] = IndexType(halo.size());
        }
    };

    template<>
    struct HaloControl<dim_3D>
    {
    };

  }
}

#endif

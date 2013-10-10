#pragma once
#ifndef KERNEL_FOUNDATION_LOAD_BALANCING_HPP
#define KERNEL_FOUNDATION_LOAD_BALANCING_HPP 1

#include<kernel/foundation/base.hpp>
#include<kernel/foundation/data.hpp>

namespace FEAST
{
  namespace Foundation
  {

    ///Example of a LBPolicy class, simply reads out comm_structures of PData object
    template<typename DT_>
    class LBPUniformCompScaledComm
    {
      public:
        ///decides how the current local comm situation has to be weighted
        template<
          typename Dim_,
          typename t_,
          template <typename, typename> class os_,
          template <typename, typename, template<typename, typename> class > class MeshType_>
        static const DT_ patch_comm_cost(const PData<Dim_, t_, os_, MeshType_, DT_>& data)
        {
          DT_ result(0);

          for(Index i(0) ; i < data.comm_halos.size() ; ++i)
            result += _latency + data.comm_halos.at(i)->size() * (DT_(Dim_::ElementPolytopeType_::tag_value)) * _bandwidth; ///estimated n edges => n + 1 verts, ...

          return result;
        }

        ///decides how the current local computational cost situation has to be weighted
        template<
          typename Dim_,
          typename t_,
          template <typename, typename> class os_,
          template <typename, typename, template<typename, typename> class > class MeshType_>
        static const DT_ patch_comp_cost(const PData<Dim_, t_, os_, MeshType_, DT_>& data)
        {
          return DT_(data.submesh->num_polytopes(Dim_::ElementPolytopeType_::tag_value));
        }

      private:
        ///T_comm = latency + N / bandwidth
        ///later these should be retrieved from the system
        constexpr static const DT_ _latency = DT_(10);  ///maybe ms
        constexpr static const DT_ _bandwidth = DT_(1); ///maybe byte/ms
    };

    template<typename LBPT_>
    struct LoadBalancing
    {
      template<
        typename Dim_,
        typename t_,
        template <typename, typename> class os_,
        template <typename, typename, template<typename, typename> class > class MeshType_,
        typename DT_>
      static void execute(PData<Dim_, t_, os_, MeshType_, DT_>& data)
      {
        DT_ comm_cost(LBPT_::patch_comm_cost(data));
        DT_ comp_cost(LBPT_::patch_comp_cost(data));

        data.comm_cost = comm_cost;
        data.comp_cost = comp_cost;

        ///TODO transform PData
      }
    };
  }
}
#endif

#pragma once
#ifndef KERNEL_FOUNDATION_LOAD_BALANCING_HH
#define KERNEL_FOUNDATION_LOAD_BALANCING_HH 1

#include<kernel/base_header.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template<typename ReturnType_>
    struct LBConfig
    {
      LBConfig(const ReturnType_& n, const ReturnType_& p) :
        network(n),
        patch_mesh(p),
        patch_process_map(ReturnType_())
      {
      }

      const ReturnType_& network;
      const ReturnType_& patch_mesh;
      ReturnType_ patch_process_map;
    };

    struct SimpleLoadBalancingPolicy
    {
      template<typename ReturnType_>
      static ReturnType_ execute(const LBConfig<ReturnType_>& lbconf)
      {

#ifndef FEAST_SERIAL_MODE

        ///TODO only lb procs, scatter, apply real LB
        int numprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

        ReturnType_ patch_process_map;

        for(Index i(0) ; i < lbconf.patch_mesh.size() ; ++i)
        {
          patch_process_map.push_back();
        }

        Index j(0);
        for(Index i(0) ; i < numprocs ; ++i)
        {
          patch_process_map.at(j).push_back(i);
          ++j;
          j = j < lbconf.patch_mesh.size() ? j : Index(0);
        }

#else

        patch_process_map.at(0).push_back(IndexType_(0));

#endif

        return patch_process_map;
      }
    };
  }
}
#endif

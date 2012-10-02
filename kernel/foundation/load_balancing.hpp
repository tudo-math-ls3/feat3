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
      ReturnType_ process_patch_map;
    };

    struct SimpleLoadBalancingPolicy
    {
      template<typename ReturnType_>
      static void execute(LBConfig<ReturnType_>& lbconf)
      {

#ifndef FEAST_SERIAL_MODE

        ///TODO only lb procs, scatter, apply real LB
        int numprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

        for(Index i(0) ; i < lbconf.patch_mesh.size() ; ++i)
        {
          lbconf.patch_process_map.push_back();
        }

        for(Index i(0) ; i < numprocs ; ++i)
        {
          lbconf.process_patch_map.push_back();
        }

        Index j(0);
        for(Index i(0) ; i < numprocs ; ++i)
        {
          lbconf.patch_process_map.at(j).push_back(i);
          lbconf.process_patch_map.at(i).push_back(j);
          ++j;
          j = j < lbconf.patch_mesh.size() ? j : Index(0);
        }

#else
        ///every patch is processed by our only process 0
        lbconf.process_patch_map.push_back();
        for(Index i(0) ; i < lbconf.patch_mesh.size() ; ++i)
        {
          lbconf.patch_process_map.push_back();
          lbconf.patch_process_map.at(i).push_back(Index(0));
          lbconf.process_patch_map.at(0).push_back(i);
        }

#endif
      }
    };
  }
}
#endif

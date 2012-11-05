#pragma once
#ifndef KERNEL_FOUNDATION_CONTROL_HH
#define KERNEL_FOUNDATION_CONTROL_HH 1

#include<kernel/archs.hpp>
#include<kernel/foundation/topology.hpp>
#include<kernel/foundation/load_balancing.hpp>
#include<iostream>

#ifndef SERIAL
#include<mpi.h>
#endif


using namespace FEAST;
using namespace Archs;

namespace FEAST
{
  namespace Foundation
  {
    template<typename Tag_, typename LBPolicy_, typename DataFillPolicy_>
    struct Control
    {
      template<typename ConfigType_, typename PatchDataType_>
      static void init(ConfigType_ & lbconf,
                       PatchDataType_ & pd,
                       int rank)
      {
        LBPolicy_::execute(lbconf);

        DataFillPolicy_::execute(pd, rank);
      }
    };
  }
}
#endif

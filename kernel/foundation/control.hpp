#pragma once
#ifndef KERNEL_FOUNDATION_CONTROL_HH
#define KERNEL_FOUNDATION_CONTROL_HH 1

#include<kernel/archs.hpp>
#include<kernel/foundation/topology.hpp>
#include<kernel/foundation/load_balancing.hpp>
#include<iostream>

#ifndef FEAST_SERIAL_MODE
#include<mpi.h>
#endif


using namespace FEAST;
using namespace Archs;

namespace FEAST
{
  namespace Foundation
  {
    template<typename Tag_, typename LBPolicy_>
    struct Control
    {
      template<typename IndexType_, template<typename, typename> class OuterStorageType_, typename StorageType_>
      static void init(Config<Topology<IndexType_, OuterStorageType_, StorageType_> > & lbconf)
      {
        LBPolicy_::execute(lbconf);
      }
    };
  }
}
#endif

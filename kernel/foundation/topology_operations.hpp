#pragma once
#ifndef KERNEL_FOUNDATION_TOPOLOGY_OPERATIONS_HPP
#define KERNEL_FOUNDATION_TOPOLOGY_OPERATIONS_HPP 1

#include <vector>
#include<kernel/base_header.hpp>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief Erase element in Topology TODO: use different policies
     *
     * \author Markus Geveler
     */
    struct TopologyElementErasure
    {
      template<typename TopologyType_>
      static inline void execute(TopologyType_& target, typename TopologyType_::index_type_ position)
      {
        target.erase(target.begin() + position);
      }

      template<typename TopologyType_>
      static inline void execute(TopologyType_& target)
      {
        target.erase(target.end() - 1);
      }
    };
  }
}
#endif

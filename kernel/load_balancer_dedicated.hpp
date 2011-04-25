#pragma once
#ifndef KERNEL_LOAD_BAL_DEDICATED_HPP
#define KERNEL_LOAD_BAL_DEDICATED_HPP 1

// includes, system

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/process_group.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief represents a dedicated load balancer
  *
  * Currently just a stub. Not clear yet, how a dedicated load balancer process is used.
  *
  * \author Hilmar Wobker
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class LoadBalancerDedicated
  {

  private:

    /* *****************
    * member variables *
    *******************/
    /// pointer to the main process group
    ProcessGroup* _process_group;


  public:

    /* ****************************
    * constructors and destructor *
    ******************************/
    /**
    * \brief CTOR for the case the load balancer lives on the coordinator process
    *
    * In this case, the base mesh can be simply passed as pointer
    */
    LoadBalancerDedicated(ProcessGroup* process_group)
      : _process_group(process_group)
    {
      std::cout << Process::rank << " is dedicated load balancer!" << std::endl;
    }
  };
} // namespace FEAST

#endif // guard KERNEL_LOAD_BAL_DEDICATED_HPP

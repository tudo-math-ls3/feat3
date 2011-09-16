/* GENERAL_REMARK_BY_HILMAR:
 * As we discussed in our last meeting, it is not clear yet how this dedicated load balancing process is to be
 * integrated into the program. My idea was that only the (normal) load balancer "knows" this dedicated load balancer
 * and provides it with special tasks... but maybe it turns out that this is not so clever... I really don't know
 * since its not really clear to me, how these tasks exactly look like... (predicting and calculating some possible
 * scenarios based on previous iterations during the actual computation and before the current data arrives...)
 * So, Luis... this is your playground! Do whatever you like! :-)
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */

#pragma once
#ifndef KERNEL_LOAD_BALANCER_DEDICATED_HPP
#define KERNEL_LOAD_BALANCER_DEDICATED_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/process_group.hpp>

// includes, system

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

#endif // KERNEL_LOAD_BALANCER_DEDICATED_HPP

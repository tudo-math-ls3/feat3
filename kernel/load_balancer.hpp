#pragma once
#ifndef KERNEL_MANAGER_HPP
#define KERNEL_MANAGER_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>
#include <kernel/process_subgroup.hpp>
#include <kernel/base_mesh/file_parser.hpp>
#include <kernel/base_mesh/bm.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief class defining the load balancer
  *
  * One process of process group is responsible for performing load balancing (collect and process matrix patch
  * statistics, compute partitioning, define work groups, ...). This is either the coordinator process of the process
  * group (being a compute process at the same time), or a dedicated load balancer process (which does nothing else).
  * If the user decides to use such a decicated load balancer process, it is the process with the highest rank in the
  * process group.
  *
  * Currently, the plan is to create the load balancer instance only on one process (namely, the coordinator of the
  * process group, or the dedicated load balancer process if there is one). This may change in future, when it turns
  * out that it is simpler/more clever to let the load balancer live on each process of its process group.
  *
  * Tasks of the load balancer are:
  * -# Define which base mesh cells (BMCs) build a matrix patch (MP) and which MPs build a process patch (PP).
  * -# Define with how many processes (usually equals the number of virtual coarse grid matrix patches, CGMPs)
  *    to solve the global coarse grid problem and how to distribute the MPs to these processes/CGMPs. The standard
  *    case will be that the coarse grid problem is solved on one processor, i.e. one CGMP containing all MPs.
  * -# Define the necessary WorkGroup objects (e.g., one WorkGroup for the fine mesh problems and one for the coarse
  *    mesh problem) and creates corresponding MPI communicators. The workers with rank 0 within the work group (being
  *    itself a process group) are the coordinators of these work groups, which communicate with
  *    the master or the dedicated load balancer.
  * -# Eventually optimise the process topologies of the work groups by building corresponding graph structures.
  *
  * There are two cases:
  * -# There is a dedicated load balancer: The dedicated load balancer gets all necessary information from the
  *    coordinator of the process group (e.g., graph describing the neighourhood of the base mesh cells) and defines
  *    MP and PP partitioning. Then it returns the necessary information to the coordinator of the process group
  *    which organises the actual setup of work groups and process distribution.
  * -# There is no dedicated load balancer: Same as case a), but the tasks of the dedicated load balancer are
  *    performed by the coordinator of the process group.
  *
  * Example:
  *
  * Distribution of submeshes to processes A-G on different levels (note that processes A-G are not necessarily
  * disjunct, i.e., several of them can refer to the same physical process, see cases a and b):
  *
  \verbatim
  ---------------      ---------------      ---------------
  |             |      |      |      |      |      |      |
  |             |      |      |      |      |  D   |  G   |
  |             |      |      |      |      |      |      |
  |      A      |      |  B   |  C   |      ---------------
  |             |      |      |      |      |      |      |
  |             |      |      |      |      |  E   |  F   |
  |             |      |      |      |      |      |      |
  ---------------      ---------------      ---------------
    level 0               level 1              levels 2-L
  \endverbatim
  *
  * -# case a, four physical processes:
  \verbatim
  process group rank:  0  1  2  3
         WorkGroup 2:  D  E  F  G         (four WorkGroup processes for the problems on level 2-L)
         WorkGroup 1:  B     C            (two WorkGroup processes for the problem on level 1)
         WorkGroup 0:  A                  (one WorkGroup process for the coarse mesh problem on level 0)
  \endverbatim
  *
  *   Communication:
  *   A <--> B (internal, rank 0) A <--> C (external, ranks 0+2)
  *   B <--> D (internal, rank 0) B <--> E (external, ranks 0+1)
  *   C <--> F (internal, rank 2) C <--> G (external, ranks 2+3)
  *
  * -# case b, five physical processes::
  \verbatim
  process group rank:  0  1  2  3  4
         WorkGroup 2:     D  E  F  G
         WorkGroup 1:     B     C
         WorkGroup 0:  A
  \endverbatim
  *
  *   Communication:
  *   A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+3)
  *   B <--> D (internal, rank 1)    B <--> E (external, ranks 1+2)
  *   C <--> F (internal, rank 3)    C <--> G (external, ranks 3+4)
  *
  * -# case c, seven physical processes:
  \verbatim
  process group rank:  0  1  2  3  4  5  6
         WorkGroup 2:           D  E  F  G
         WorkGroup 1:     B  C
         WorkGroup 0:  A
  \endverbatim
  *
  *   Communication:
  *   A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+2)
  *   B <--> D (external, ranks 1+3) B <--> E (external, ranks 1+4)
  *   C <--> F (external, ranks 2+5) C <--> G (external, ranks 2+6)
  *
  * \tparam space_dim_
  * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
  *
  * \tparam world_dim_
  * world dimension (determines the number of coordinates)
  *
  * \todo This class is only a dummy with some comments! No real load balancing functionalities implemented yet!
  *
  * \author Hilmar Wobker
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class LoadBalancer
  {

  private:

    /* *****************
    * member variables *
    *******************/

    // ...

    /**
    * \brief dummy function in preparation of a function for defining work groups
    *
    * As long as we cannot do this in an automatic way, the test driver calls a corresponding function, which sets
    * work groups manually.
    */
    void define_work_groups()
    {
      CONTEXT("LoadBalancer::define_work_groups()");
    }
  };
} // namespace FEAST

#endif // guard KERNEL_MANAGER_HPP

#pragma once
#ifndef KERNEL_INTER_LEVEL_GROUP_HPP
#define KERNEL_INTER_LEVEL_GROUP_HPP 1

// includes, system
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <mpi.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/logger.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief represents the group of processes from work groups of different levels that interact with each other
  *
  * Illustration of inter-level work group communication:
  *
  *   ---------------        ---------------
  *   |      |      |        |      |      |
  *   |      |      |        |  E   |  F   |
  *   |      |      |        |      |      |
  *   |  A   |  B   |        ---------------
  *   |      |      |        |      |      |
  *   |      |      |        |  C   |  D   |
  *   |      |      |        |      |      |
  *   ---------------        ---------------
  *   coarse grid (CG)       fine grid (FG)
  *
  *   A - F = compute tasks, each performed on one process
  *
  *   coarse grid work group: A,B
  *   fine grid work group: C, D, E, F
  *
  *   necessary data exchange for restriction and prolongation:   A <-> (C|E),     B <-> (D|F)
  *
  *   exemplary distribution of processes and corresponding program flow:
  *
  *       processes:  0   1   2   3   4
  *   CG work group:  A   B
  *   FG work group:  C       D   E   F
  *
  * step  0(A+C)                 1(B)                  2(D),                3(E)                 4(F)
  * 0     solver.start()         solver.start()        solver.start()       solver.start()       solver.start()
  * 1     smooth()               recv_defect(2,4)      smooth()             smooth()             smooth()
  * 2     restrict(FG->CG)         *wait*              restrict(FG->CG)     restrict(FG->CG)     restrict(FG->CG)
  * 3     send_defect(0)           *wait*              send_defect(1)       send_defect(0)       send_defect(1)
  * 4     recv_defect(0,3)         *wait*              recv_corr(1)         recv_corr(0)         recv_corr(1)
  * 5     solve_coarse()         solve_coarse()          *wait*               *wait*               *wait*
  * 6     send_corr(0,3)         send_corr(2,4)          *wait*               *wait*               *wait*
  * 7     recv_corr(0)           recv_defect(2,4)      prolong(FG->CG)      prolong(FG->CG)      prolong(FG->CG)
  * 8     prolong(CG -> FG)        *wait*                *wait*               *wait*               *wait*
  * 9     smooth()                 *wait*              smooth()             smooth()             smooth()
  *
  * send_defect(i) means 'send defect vector to process i'
  * send_corr(i,j) means 'send parts of the correction vector corresponding to the region scheduled to processes
  *                       i and j, resp., to these processes'
  *
  * Example (see sketch above): C and E compute their local portions of the global defect and sent it to A. A assembles
  * these two portions and computes a correction vector. The two parts of this correction vector corresponding to the
  * lower left and upper left subregion, resp., are then sent back to C and E.
  *
  * We require that a region, which is scheduled to one process on grid level L is also scheduled to exactly one process
  * on all coarser levels l < L. On a finer level l > L, however, the region may be scheduled to more than one process.
  *
  * Consequence: A process sends its portion of the defect vector to exactly one process (possibly the same). It never
  * occurs that the defect vector to be sent has to be distributed over more than one process. (The same is analogously
  * true for receiving correction vectors.)
  * On the other hand, defects can be received from several processes (where one of these processes may be the same as
  * the receiving one) and corrections can be sent to several processes.
  *
  * All MPI communication is done via extra communicators containing only the processes of the different work groups
  * that have to exchange data with each other. In the example above, these are:
  *   communicator 1:  0 (CG-WG + FG-WG), 3 (FG-WG)
  *   communicator 2:  1 (CG-WG), 2, 4 (FG-WG)
  * This way, all necessary data exchange can be performed efficiently via corresponding collective communication
  * routines (e.g., receive portions of the defect vector via MPI_Gather and send portions of the correction vector via
  * MPI_Scatter, where the CG-WG process is the root, resp.).
  * This is possible since we always have a 1:n communication pattern (one processor on the coarser level communicates
  * with n processes on the finer level).
  * We don't have to care about general m:n communication patterns since each fine grid work group process is related to
  * exactly one coarse grid work group process. Hence, each m:n relation decomposes into m 1:n_k, k=1,..,m, relations
  * where n_1 + ... + n_k = n.
  * (In the example above, we can treat two 1:2 relations in parallel instead of one 2:4 relation.)
  *
  * The operations (send|recv)_(defect|corr) will probably be part of the grid transfer routines restrict/prolong
  * (and will be called send_restricted_vector() etc. or similar). When two compute tasks are performed on the same
  * process (e.g., process 0(A+C) in the example above), then the send/recv routines have to perform corresponding
  * copy operations instead of calling MPI. This sort of data transfer is not performed via the class InterLevelGroup.
  * COMMENT_HILMAR: Or maybe it will? Collective routines can also send data to itself... let's see what is more
  *   appropriate...
  *
  * \note When realising the n-layer-ScaRC concept, we will probably need another sort of inter work group
  * communication. Imagine a unitsquare consisting of 4x4 BMCs and a 3-layer-ScaRC scheme.
  *   Layer 1: global MG, acting on all 16 BMCs, one work group consisting of 16 processes
  *   Layer 2: four MGs acting on 2x2 BMCs each, four work groups consisting of 4 processes each
  *   Layer 3: 16 local MGs acting on 1 BMC each
  * The four work groups on layer 2 have to communicate with each other over their common boundary. For this, other
  * communication patterns are needed than for the inter level communication considered here.
  */
  class InterLevelGroup
    : public ProcessGroup
  {

  private:
    /**
    * \brief rank of the root process for all MPI routines
    *
    * rank of the process from the coarser grid level work group which is the
    */
    int _root;

//    ...

  public:
    InterLevelGroup(...)
    {
//      ...
    }

    ~InterLevelGroup()
    {
//      ...
    }

//    ...
  };
} // namespace FEAST

#endif // guard KERNEL_INTER_LEVEL_GROUP_HPP

/* GENERAL_REMARK_BY_HILMAR:
 * Here, one has to find out if it *really* is no problem, when one process is member of two interlevel groups
 * (see COMMENT_HILMAR below). The main challenge regarding interlevel groups is to dynamically determine the fine-
 * and coarse level groups and the corresponding interlevel groups.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_INTERLEVEL_GROUP_HPP
#define KERNEL_INTERLEVEL_GROUP_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/logger.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>

// includes, system
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <mpi.h>

namespace FEAST
{
  /**
  * \brief represents the group of processes from work groups of different levels that interact with each other
  *
  * Illustration of interlevel work group communication:
  * A - F = compute tasks, each performed on one process
  *
  * <bf>Simplemost example:</bf>
  *
  \verbatim
  ---------------        ---------------
  |             |        |      |      |
  |             |        |  C   |  D   |
  |             |        |      |      |
  |      A      |        ---------------
  |             |        |      |      |
  |             |        |  B   |  E   |
  |             |        |      |      |
  ---------------        ---------------
  coarse grid (CG)       fine grid (FG)
  \endverbatim
  *
  *   coarse grid work group: A
  *   fine grid work group: B,C,D,E
  *
  * necessary data exchange for restriction and prolongation: A <-> (B,C,D,E)
  *
  *   exemplary distribution of processes and corresponding program flow:
  *
  *       processes:  0   1   2   3   4
  *   CG work group:  A
  *   FG work group:      B   C   D   E
  *
  * step  0(A)                    1(B)                  2(C),                3(D)                 4(E)
  * 0     recv_defect(1,2,3,4)    solver.start()        solver.start()       solver.start()       solver.start()
  * 1     *wait*                  smooth()              smooth()             smooth()             smooth()
  * 2     *wait*                  restrict(FG->CG)      restrict(FG->CG)     restrict(FG->CG)     restrict(FG->CG)
  * 3     *wait*                  send_defect(0)        send_defect(0)       send_defect(0)       send_defect(0)
  * 4     solve_coarse()          recv_corr(0)          recv_corr(0)         recv_corr(0)         recv_corr(0)
  * 5     send_corr(0,1,2,3)      *wait*                *wait*               *wait*               *wait*
  * 6     recv_defect(1,2,3,4)(a) prolong(CG->FG)       prolong(CG->FG)      prolong(CG->FG)      prolong(CG->FG)
  * 7     *wait*                  smooth()              smooth()             smooth()             smooth()
  *
  * (a) already for the next coarse grid solution
  *
  * send_defect(i) means 'send defect vector to process i'
  * recv_defect(i,...) means 'receive defect vector from processes i,...
  * send_corr(i,j) means 'send parts of the correction vector corresponding to the region scheduled to processes
  *                       i and j, resp., to these processes'
  *
  * Explanation: one dedicated process to solve the coarse grid problem serially, while all fine grid processes
  * are idle and vice versa. More technical details below in the more complex example.
  *
  *
  * <bf>More complicated example:</bf>
  *
  \verbatim
  ---------------        ---------------
  |      |      |        |      |      |
  |      |      |        |  E   |  F   |
  |      |      |        |      |      |
  |  A   |  B   |        ---------------
  |      |      |        |      |      |
  |      |      |        |  C   |  D   |
  |      |      |        |      |      |
  ---------------        ---------------
  coarse grid (CG)       fine grid (FG)
  \endverbatim
  *
  *   coarse grid work group: A, B
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
  * 7     recv_corr(0)           recv_defect(2,4)(a)   prolong(CG->FG)      prolong(CG->FG)      prolong(CG->FG)
  * 8     prolong(CG -> FG)        *wait*                *wait*(b)             *wait*(b)            *wait*(b)
  * 9     smooth()                 *wait*              smooth()             smooth()             smooth()
  *
  * (a) already for the next coarse grid solution
  * (b) waiting because process 0 is not yet ready with prolongation, and neighbour exchange is necessary
  *     after each prolongation step
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
  * copy operations instead of calling MPI. This sort of data transfer is not performed via the class InterlevelGroup.
  * COMMENT_HILMAR: Or maybe it will? Collective routines can also send data to itself... let's see what is more
  *   appropriate...
  *
  * COMMENT_HILMAR: Imagine the following modification of the example above:
  *       processes:  0   1   2   3   4
  *   CG work group:  A   B
  *   FG work group:  D       C   E   F
  * i.e., A and D are on the same process, but they don't have talk to each other in interlevel communication. Then
  * the two interlevel groups are given by:
  *   communicator 1:  0 (CG-WG), 2, 3 (FG-WG)
  *   communicator 2:  1 (CG-WG), 0, 4 (FG-WG)
  * i.e., they are not disjunct (process 0 appears in both). As far as I understand MPI, it should not be a problem to
  * start communication in both groups simultaneously: As long as the first 'blocks' process 0, the second has to
  * wait. Of course, this is not efficient, since the communications have to wait on each other. Hence, it is advisable
  * to create disjunct interlevel groups when possible.
  *
  * \note When realising the n-layer-ScaRC concept, we will probably need another sort of inter work group
  * communication. Imagine a unitsquare consisting of 4x4 BMCs and a 3-layer-ScaRC scheme.
  *   Layer 1: global MG, acting on all 16 BMCs, one work group consisting of 16 processes
  *   Layer 2: four MGs acting on 2x2 BMCs each, four work groups consisting of 4 processes each
  *   Layer 3: 16 local MGs acting on 1 BMC each
  * The four work groups on layer 2 have to communicate with each other over their common boundary. For this, other
  * communication patterns are needed than for the interlevel communication considered here.
  *
  * \note (dom, Jul 20, 2011) This is closely related to the matrix patch layer, and not to the coarse/fine distribution.
  * Depending on how we implement this, the ScaRC layers define the "outer" communication, and depending on the sizes
  * and on the available resources, each such layer may or may not choose to perform their coarse grid solution using
  * separate process subsets.
  */
  class InterlevelGroup
    : public ProcessGroup
  {

  private:
    /**
    * \brief rank of the root process (w.r.t. to the InterlevelGroup's communicator) for collective MPI routines
    *
    * The root process is always the process from the coarser grid level work group.
    */
    int _rank_root;

//    ...

  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /**
    * \brief CTOR
    *
    * \param[in] comm
    * communicator shared by the group processes
    *
    * \param[in] group_id
    * ID of this group
    *
    * \param[in] rank_root
    * rank of the root process for collective MPI routines
    *
    */
    InterlevelGroup(
      MPI_Comm comm,
      unsigned int const group_id,
      int rank_root)
      : ProcessGroup(comm, group_id),
        _rank_root(rank_root)
    {
//      ...

    }

    ~InterlevelGroup()
    {
    }

    /* *****************
    * member functions *
    *******************/
    /**
    * \brief checks whether this process is the root of the interlevel group
    *
    * \return true if this process is the root of the interlevel group, otherwise false
    */
    inline bool is_root() const
    {
      CONTEXT("InterlevelGroup::is_root()");
      return _rank == _rank_root;
    }


    /**
    * \brief dummy function for testing communication within the interlevel group
    *
    * The function simply performs a gather operation, where the root process gathers integers from the non-root
    * processes.
    *
    * COMMENT_HILMAR: Actually, I intended to make this function private and declare the functions
    *     ManagerCompCoord<space_dim_, world_dim_>::_test_communication()
    *   and
    *     ManagerCompNonCoord<space_dim_, world_dim_>::_test_communication()
    *   as friends. But since the InterlevelGroup class doesn't know the template parameters space_dim_ and world_dim_,
    *   this is not possible. Maybe there is a way to realise it, but I didn't have the time to find out how.
    *
    * \return flag whether test was succesful (0: succesful, >0: failed)
    */
    bool test_communication()
    {
      CONTEXT("InterlevelGroup::test_communication()");

      // the value to be gathered is the world rank of this process
      int world_rank = Process::rank;

      if(!is_root())
      {
        // non-root process

        // debug output
        Logger::log("Interlevel group " + stringify(_group_id) + ": Process " + stringify(_rank)
                     + " sends its world rank " + stringify(Process::rank) + " to the root process.\n");
        // non-root process sends it world rank to the root
        MPI_Gather(&world_rank, 1, MPI_INTEGER, nullptr, 0, MPI_DATATYPE_NULL, _rank_root, _comm);

        // additionally perform a reduce operation (sum up all world ranks)
        MPI_Reduce(&world_rank, nullptr, 1, MPI_INTEGER, MPI_SUM, _rank_root, _comm);

        // on non-root side there is nothing to check, so simply return 0 (test successful)
        return 0;
      }
      else
      {
        // root process

        // array on the root process for collecting the data from the other processes
        int* recv_ranks = new int[num_processes()];

        // root process collects all world ranks from the non-roots
        MPI_Gather(&world_rank, 1, MPI_INTEGER, recv_ranks, 1, MPI_INTEGER, _rank_root, _comm);

        // debugging output
        std::string s(stringify(recv_ranks[0]));
        int sum(recv_ranks[0]);
        for(unsigned int i(1) ; i < num_processes(); ++i)
        {
          s += " " + stringify(recv_ranks[i]);
          sum += recv_ranks[i];
        }

        int sum_reduce;
        // additionally perform a reduce operation (sum up all world ranks)
        MPI_Reduce(&world_rank, &sum_reduce, 1, MPI_INTEGER, MPI_SUM, _rank_root, _comm);

        Logger::log("Interlevel group " + stringify(_group_id) + ": Root process " + stringify(_rank)
                    + " gathered world ranks " + s + ", summing up to " + stringify(sum_reduce) + "\n");
        delete [] recv_ranks;

        // We cannot really do a clever check here, so we simply check wether sums are equal.
        // (The general problem is: When there is something wrong with the communication, then usually MPI
        // crashes completely. On the other hand, when communication is fine then usually correct values are sent.
        // So, it is very unlikely that communication works *AND* this test here returns false.)
        if(sum_reduce == sum)
        {
          return 0;
        }
        else
        {
          return 1;
        }
      }
    } // test_communication()
  };
} // namespace FEAST

#endif // KERNEL_INTERLEVEL_GROUP_HPP

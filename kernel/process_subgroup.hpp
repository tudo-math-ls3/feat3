#pragma once
#ifndef KERNEL_BM_PROCESS_SUBGROUP_HPP
#define KERNEL_BM_PROCESS_SUBGROUP_HPP 1

// includes, system
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <mpi.h>
#include <math.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/logger.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>
#include <kernel/work_group.hpp>
#include <kernel/graph.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief class describing a subgroup of a process group, consisting of some compute processes and eventually one
  *        extra coordinator process, sharing the same MPI communicator
  *
  * ProcessSubgroup objects are created by the manager. They consist of n compute processes and eventually one
  * extra process which is the coordinator of the parent process group. The latter is only the case if the coordinator
  * is not part of the compute processes anyway. In both cases (containing an extra coordinator process or not), a
  * ProcessSubgroup creates a WorkGroup object, which either consists of exactly the same processes as this
  * ProcessSubgroup (in the case there is no extra coordinator process) or of the n compute processes only (excluding
  * the extra coordinator process). Each ProcessSubgroup has its own MPI communicator. Since the coordinator of the
  * parent process group is part of this communicator, all necessary information (mesh, graph, ...) can be transferred
  * efficiently via collective communication routines. The extra WorkGroup object (excluding the eventual coordinator
  * process) with its own communicator is necessary to efficiently perform collective communication during the actual
  * computation (scalar products, norms, etc).
  *
  * Example:
  * The process group of a manager consists of six processes, the first one being the coordinator of the process
  * group. There is no dedicated load balancing process. The coordinator process (rank 0) reads the mesh and the solver
  * configuration and decides that the coarse grid problem is to be treated by two compute processes (process group
  * ranks 4 and 5) and the fine grid problems by six compute processes (ranks 0-5). Then two ProcessSubgroup objects are
  * created: The first one (for the coarse grid problem) consists of the two compute processes with ranks 4 and 5 and
  * the coordinator of the parent process group (rank 0). The other ProcessSubgroup (for the fine grid problem) object
  * consists of all six processes (ranks 0-5). Here, the coordinator of the parent process group (rank 0d) is also a
  * compute process. Both ProcessSubgroups then create a WorkGroup object. The coarse grid work group consists of the
  * two compute process (ranks 4 and 5) and thus differs from its ProcessSubgroup, while the fine grid work group
  * contains exactly the same six processes as its ProcessSubgroup.
  *
  *    processes parent process group:  0 1 2 3 4 5
  *    coarse grid process subgroup:    x       x x
  *    coarse grid work group:                  x x
  *    fine grid process subgroup:      x x x x x x
  *    fine grid work group:            x x x x x x
  *
  * If the parent process group contains seven processes with the seventh process being a dedicated load balancer
  * (rank 6), the process subgroups and work groups look like this:
  *
  *    processes parent process group:  0 1 2 3 4 5 6
  *    coarse grid process subgroup:    x       x x
  *    coarse grid work group:                  x x
  *    fine grid process subgroup:      x x x x x x
  *    fine grid work group:            x x x x x x
  *
  * I.e., the dedicated load balancing process does not belong to any of the subgroups or work groups.
  *
  * \author Hilmar Wobker
  * \author Dominik Goeddeke
  *
  */
  class ProcessSubgroup
    : public ProcessGroup
  {

  private:

    /* *****************
    * member variables *
    *******************/
    /// flag whether this group contains the coordinator of the parent process group as an extra process
    bool _contains_extra_coord;

    /**
    * \brief work group containing only the real compute processes of this ProcessSubgroup
    *
    * If #_contains_extra_coord == false, then the work group contains all processes of this ProcessGroup, otherwise it
    * it contains only the compute processes excluding the the extra coordinator process. The parent process group of
    * the work group is not this ProcessSubgroup, but the parent group of this ProcessSubgroup. That is why, when
    * creating the work group, the dummy calls of MPI_Comm_create(...) by the remaining processes of the parent process
    * group are performed in the manager class and not here.
    */
    WorkGroup* _work_group;


  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /**
    * \brief CTOR
    *
    * \param[in] num_processes
    * number of processes in this group
    *
    * \param[in] ranks_group_parent
    * ranks in the parent group
    *
    * \param[in] process_group_parent
    * parent group of processes
    *
    * \param[in] group_id
    * ID of this group
    *
    * \param[in] contains_extra_coord
    * flag whether this group contains the coordinator of the parent process group as an extra process
    */
    ProcessSubgroup(
      unsigned int const num_processes,
      int* const ranks_group_parent,
      ProcessGroup* const process_group_parent,
      unsigned int const group_id,
      bool contains_extra_coord)
      : ProcessGroup(num_processes, ranks_group_parent, process_group_parent, group_id),
        _contains_extra_coord(contains_extra_coord),
        _work_group(nullptr)
    {
      CONTEXT("ProcessSubgroup::ProcessSubgroup()");
      /* ******************************
      * test the logger functionality *
      ********************************/

//TODO: remove this debug output

      // debugging output
    // let the coordinator of the subgroup trigger some common messages
//      if(is_coordinator())
//      {
//        std::string s("I have COMM_WORLD rank " + stringify(Process::rank)
//                      + " and I am the coordinator of process subgroup " + stringify(_group_id));
//   //      Logger::log_master("Hello, master screen! " + s, Logger::SCREEN);
//        Logger::log_master("Hello, master file! " + s, Logger::FILE);
//   //      Logger::log_master("Hello, master screen and file! " + s, Logger::SCREEN_FILE);
//   //      Logger::log_master("Hello, default master screen and file! " + s);
//      }

      // write some individual messages to screen and file
      std::string s("I have COMM_WORLD rank " + stringify(Process::rank) + " and group rank " + stringify(_rank)
                    + " in process subgroup " + stringify(_group_id) + ".");
      if(is_coordinator())
      {
        s += " I am the coordinator!";
      }
      // vary the lengths of the messages a little bit
      for(int i(0) ; i < _rank+1 ; ++i)
      {
        s += " R";
      }
      for(unsigned int i(0) ; i < _group_id+1 ; ++i)
      {
        s += " G";
      }
      s += "\n";
//      log_indiv_master("Hello, master screen! " + s, Logger::SCREEN);
      log_indiv_master("Hello, master file! " + s, Logger::FILE);
//      log_indiv_master("Hello, master screen and file! " + s, Logger::SCREEN_FILE);
//      log_indiv_master("Hello, master default screen and file! " + s);

      // end of debugging output

      // Create work group consisting of the real compute processes only (the dummy calls of MPI_Comm_create(...) by the
      // remaining processes of the parent process group are performed in the manager class and not here). The group id
      // of the work group is the same as that of the process subgroup.
      if (!_contains_extra_coord)
      {
        // in the case there is no extra coordinator process, the work group of compute processes contains all processes
        // of this subgroup
        _work_group = new WorkGroup(_num_processes, _ranks_group_parent, _process_group_parent, _group_id);
      }
      else
      {
        // otherwise, the work group contains all processes of this subgroup except the coordinator process

        // We assume here that the coordinator of this process subgroup is also the coordinator of the parent group.
        // Since we always set the process with rank 0 to be the coordinator, this condition is fulfilled automatically.
        // However, when the latter is changed for some reason, the condition might be violated. In this case, we have
        // to provide the possibility to determine which process exactly should be the coordinator when creating
        // a process group (instead of unconditionally setting the first process to be the coordinator).
        if(is_coordinator())
        {
          ASSERT(_process_group_parent->is_coordinator(),
                 "Coordinator process of the subgroup has to be coordinator of the parent, too!");
        }

        if(!_process_group_parent->is_coordinator())
        {
          // copy all ranks of the parent group except its coordinator to an auxiliary array
          int* ranks_group_parent_aux = new int[_num_processes - 1];
          unsigned int shift(0);
          for(unsigned int i(0) ; i < _num_processes ; ++i)
          {
            if(_ranks_group_parent[i] != _process_group_parent->rank_coord())
            {
              ranks_group_parent_aux[i - shift] = _ranks_group_parent[i];
            }
            else
            {
              shift = 1;
            }
          }
          // assert that the coordinator of the parent group has actually been found
          ASSERT(shift == 1, "Coordinator in parent process group not found!");
          // create work group using the aux. array of parent ranks
          _work_group = new WorkGroup(_num_processes - 1, ranks_group_parent_aux, _process_group_parent, _group_id);
          // delete the aux array again
          delete [] ranks_group_parent_aux;
        }
        else
        {
          // *All* processes of the parent MPI group have to call the MPI_Comm_create() routine (otherwise the forking
          // will deadlock), so let the coordinator call the routine with dummy communicator and dummy group.  (The
          // dummy group is necessary here since the other MPI_Group object is hidden inside the WorkGroup CTOR above.)
          MPI_Comm dummy_comm;
          MPI_Group dummy_group;
          int rank_aux = _process_group_parent->rank();
          int mpi_error_code = MPI_Group_incl(_process_group_parent->group(), 1, &rank_aux, &dummy_group);
          validate_error_code_mpi(mpi_error_code, "MPI_Group_incl");
          mpi_error_code = MPI_Comm_create(_process_group_parent->comm(), dummy_group, &dummy_comm);
          validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");
          // COMMENT_HILMAR: First, I used this simpler version:
          //   mpi_error_code = MPI_Comm_create(_process_group_parent->comm(), MPI_GROUP_EMPTY, &dummy_comm);
          // It worked with OpenMPI 1.4.2 and MPICH2, but does not with OpenMPI 1.4.3. We are not quite sure yet, if that
          // is a bug in OpenMPI 1.4.3, or if this use of MPI_GROUP_EMPTY is incorrect.

          MPI_Comm_free(&dummy_comm);
          MPI_Group_free(&dummy_group);
          _work_group = nullptr;
        }
      }
    } // constructor


    /// DTOR (automatically virtual since DTOR of base class ProcessGroup is virtual)
    ~ProcessSubgroup()
    {
      CONTEXT("ProcessSubgroup::~ProcessSubgroup()");
      if(_work_group != nullptr)
      {
        delete _work_group;
        _work_group = nullptr;
      }
    }


    /**
    * \brief getter for the flag whether this subgroup contains an extra coordinator process
    *
    * \return pointer the flag #_contains_extra_coord
    */
    inline bool contains_extra_coord() const
    {
      CONTEXT("ProcessSubgroup::contains_extra_coord()");
      return _contains_extra_coord;
    }


    /**
    * \brief getter for the compute work group
    *
    * \return pointer to the work group #_work_group
    */
    inline WorkGroup* work_group() const
    {
      CONTEXT("ProcessSubgroup::work_group()");
      return _work_group;
    }
  };
} // namespace FEAST

#endif // guard KERNEL_BM_PROCESS_SUBGROUP_HPP

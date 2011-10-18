/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR below.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_WORK_GROUP_EXT_HPP
#define KERNEL_WORK_GROUP_EXT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/graph.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>
#include <kernel/work_group.hpp>

// includes, system
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <mpi.h>
#include <math.h>

namespace FEAST
{
  /**
  * \brief describes an extended work group, consisting of some compute processes (the actual work group) and
  *        eventually one extra coordinator process, sharing the same MPI communicator
  *
  * WorkGroupExt objects are created by the manager. They consist of n compute processes (the actual work group) and
  * eventually one extra process which is the coordinator of the parent process group. The latter is only the case if
  * the coordinator is not part of the compute processes anyway. In both cases (containing an extra coordinator process
  * or not), a WorkGroupExt creates a WorkGroup object, which either consists of exactly the same processes as this
  * WorkGroupExt (in the case there is no extra coordinator process) or of the n compute processes only (excluding
  * the extra coordinator process). Each WorkGroupExt has its own MPI communicator. Since the coordinator of the
  * parent process group is part of this communicator, all necessary information (mesh, graph, ...) can be transferred
  * efficiently via collective communication routines (which is actually the reason why this WorkGroupExt exists!).
  * The WorkGroup object (excluding the eventual coordinator process) with its own communicator is necessary to
  * efficiently perform collective communication during the actual computation (scalar products, norms, etc).
  *
  * Example:
  * The process group of a manager consists of six processes, the first one being the coordinator of the process
  * group. There is no dedicated load balancing process. The coordinator process (rank 0) reads the mesh and the solver
  * configuration and decides that the coarse grid problem is to be treated by two compute processes (process group
  * ranks 4 and 5) and the fine grid problems by six compute processes (ranks 0-5). Then two WorkGroupExt objects are
  * created: The first one (for the coarse grid problem) consists of the two compute processes with ranks 4 and 5 and
  * the coordinator of the parent process group (rank 0). The other WorkGroupExt (for the fine grid problem) object
  * consists of all six processes (ranks 0-5). Here, the coordinator of the parent process group (rank 0) is also a
  * compute process. Both WorkGroupExts then create a WorkGroup object. The coarse grid work group consists of the
  * two compute process (ranks 4 and 5) and thus differs from its extended work group, while the fine grid work group
  * contains exactly the same six processes as its extended work group.
  \verbatim
      processes parent process group:  0 1 2 3 4 5
      coarse grid ext. work group:     x       x x
      coarse grid work group:                  x x
      fine grid ext. work group:       x x x x x x
      fine grid work group:            x x x x x x
  \endverbatim
  * If the parent process group contains seven processes with the seventh process being a dedicated load balancer
  * (rank 6), the ext. work groups and work groups look like this:
  \verbatim
      processes parent process group:  0 1 2 3 4 5 6
      coarse grid ext. work group:     x       x x
      coarse grid work group:                  x x
      fine grid ext. work group:       x x x x x x
      fine grid work group:            x x x x x x
  \endverbatim
  * I.e., the dedicated load balancing process does not belong to any of the extended work group or work groups.
  *
  * \author Hilmar Wobker
  *
  */
  class WorkGroupExt
    : public ProcessGroup
  {

  private:

    /* *****************
    * member variables *
    *******************/
    /// flag whether this group contains the coordinator of the parent process group as an extra process
    bool _contains_extra_coordinator;

    /**
    * \brief work group containing only the real compute processes of this WorkGroupExt
    *
    * If #_contains_extra_coordinator == false, then the work group contains all processes of this ProcessGroup,
    * otherwise it contains only the compute processes excluding the the extra coordinator process. The parent process
    * group of the work group is not this WorkGroupExt, but the parent group of this WorkGroupExt. That is why, when
    * creating the work group, the dummy calls of MPI_Comm_create(...) by the remaining processes of the parent process
    * group are performed in the manager class and not here.
    */
    WorkGroup* _work_group;


  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /**
    * \brief CTOR which immediately constructs the work group belonging to this extended work group
    *
    * \param[in] comm
    * communicator shared by the group processes
    *
    * \param[in] group_id
    * ID of this group
    *
    * \param[in] contains_extra_coordinator
    * flag whether this group contains the coordinator of the parent process group as an extra process
    *
    * \param[in] is_coord_in_parent
    * flag whether this process is the coordinator of the parent process group (actually only needed for ASSERT)
    */
    WorkGroupExt(
      MPI_Comm comm,
      unsigned int const group_id,
      bool const contains_extra_coordinator,
      bool is_coord_in_parent)
      : ProcessGroup(comm, group_id),
        _contains_extra_coordinator(contains_extra_coordinator),
        _work_group(nullptr)
    {
      CONTEXT("WorkGroupExt::WorkGroupExt()");
      // Create work group consisting of the real compute processes only (the dummy calls of MPI_Comm_create(...) by the
      // remaining processes of the parent process group are performed in the manager class and not here). The group id
      // of the work group is the same as that of the extended work group.
      if (!_contains_extra_coordinator)
      {
        // in the case there is no extra coordinator process, the work group of compute processes contains all processes
        // of this extended work group
        MPI_Comm gr_comm;
        MPI_Comm_dup(_comm, &gr_comm);
        _work_group = new WorkGroup(gr_comm, _group_id);
        // COMMENT_HILMAR: One could also simply pass _comm itself to the WorkGroup CTOR instead of making a explicit
        //   copy of the communicator. But then we have problems when deleting the work group in the DTOR again: When
        //   'delete _work_group' is called, the communicator of the WorkGroupExt class would be deleted already which
        //   leads to a runtime error. Currently, I don't have a better idea, how to solve this problem.
      }
      else
      {
        // otherwise, the work group contains all processes of this extended work group except the coordinator process

        // We assume here that the coordinator of this extended work group is also the coordinator of the parent group.
        // Since the coordinator of the parent group is always part of this extended work group and we always set the
        // process with rank 0 to be the coordinator, this condition is fulfilled automatically.
        // However, when the latter is changed for some reason, the condition might be violated. In this case, we have
        // to provide the possibility to determine which process exactly should be the coordinator when creating
        // a process group (instead of unconditionally setting the first process to be the coordinator).
// COMMENT_HILMAR: A problem discovered by Dirk: If the code is compiled without FEAST_DEBUG_MODE, but with flags
//   -Wall -Wextra, then we get a warning (error when -Werror is used) that the paramater is_coord_in_parent is not used
//   (since the ASSERT call is removed). How should we deal with such a situation in general?
//   Regarding this special situation here: It's essential that the coordinator of the extended work group is also the
//   coordinator of the parent process group. Currently, this is guaranteed (see comment above). But I'm not sure if
//   this doesn't change in future. So, I would like to keep this inquiry here to be on the safe side!
        ASSERT(is_coord_in_parent == is_coordinator(),
               "Coordinator process of the extended work group has to be coordinator of the parent, too!");

        MPI_Group gr_without_coord;
        MPI_Comm gr_comm;
        int mpi_error_code = MPI_Group_excl(_group, 1, const_cast<int*>(&_rank_coord), &gr_without_coord);
        validate_error_code_mpi(mpi_error_code, "MPI_Group_excl");
        mpi_error_code = MPI_Comm_create(_comm, gr_without_coord, &gr_comm);
        validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");

        if(!is_coordinator())
        {
          _work_group = new WorkGroup(gr_comm, _group_id);
        }
      }
    } // constructor



    /// DTOR (automatically virtual since DTOR of base class ProcessGroup is virtual)
    ~WorkGroupExt()
    {
      CONTEXT("WorkGroupExt::~WorkGroupExt()");
      if(_work_group != nullptr)
      {
        delete _work_group;
        _work_group = nullptr;
      }
    }


    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the flag whether this extended work group contains an extra coordinator process
    *
    * \return pointer the flag #_contains_extra_coordinator
    */
    inline bool contains_extra_coordinator() const
    {
      CONTEXT("WorkGroupExt::contains_extra_coordinator()");
      return _contains_extra_coordinator;
    }


    /**
    * \brief getter for the compute work group
    *
    * \return pointer to the work group #_work_group
    */
    inline WorkGroup* work_group() const
    {
      CONTEXT("WorkGroupExt::work_group()");
      return _work_group;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief checks whether this process is the extra coordinator process of the extended work group
    *
    * \return true if this process is the coordinator of the process group and the process group contains an extra
    *         coordinator process, otherwise false
    */
    inline bool is_extra_coordinator() const
    {
      CONTEXT("ProcessGroup::is_extra_coord()");
      return is_coordinator() && contains_extra_coordinator();
    }


  };
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_WORK_GROUP_EXT_HPP

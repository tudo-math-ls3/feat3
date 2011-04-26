#pragma once
#ifndef KERNEL_MANAGER_COMP_HPP
#define KERNEL_MANAGER_COMP_HPP 1

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
#include <kernel/work_group_ext.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief defines manager of a compute process group (super class of coordinator and non-coordinator versions)
  *
  * \author Hilmar Wobker
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class ManagerComp
  {

  protected:

    /* *****************
    * member variables *
    *******************/
    /// process group the manager manages
    ProcessGroup* _process_group;

    /**
    * \brief rank of the load balancer process
    *
    * Coincides with the coordinator process, i.e. it has rank 0.
    */
    int _rank_load_balancer;

    /// vector of extended work groups the manager manages
    std::vector<WorkGroupExt*> _work_groups;

    /// number of work groups
    unsigned int _num_work_groups;

    /**
    * \brief array of number of workers in each extended work group (including the eventual extra coordinator process)
    *
    * Dimension: [#_num_work_groups]
    */
    unsigned int* _num_proc_in_ext_work_group;

    /**
    * \brief array indicating whether the extended work groups contain an extra process for the coordinator
    *
    * This coordinator will not be a compute process in the corresponding work group.
    * Usually a boolean would do the trick, but we cannot send boolean arrays via MPI. (C++ bindings are deprecated,
    * hence we cannot use MPI::BOOL. MPI_LOGICAL cannot be used, since this is not equivalent to C++'s bool.)
    * So, we use unsigned char here and treat is as if it were boolean.
    */
    unsigned char* _group_contains_extra_coord;

    /**
    * \brief 2-dim. array of process group ranks building the extended work groups
    *
    * Dimension: [#_num_work_groups][#_num_proc_in_ext_work_group[\a group_id]]
    */
    int** _ext_work_group_ranks;

    /**
    * \brief boolean array indicating to which work groups this process belongs
    *
    * Dimension: [#_num_work_groups]
    */
    bool* _belongs_to_group;


    /// number of inter level groups
    unsigned int _num_inter_level_groups;

    /**
    * \brief array of number of workers in each inter level group
    *
    * Dimension: [#_num_inter_level_groups]
    */
    unsigned int* _num_proc_in_inter_level_group;

    /**
    * \brief 2-dim. array of process group ranks building the inter level groups
    *
    * Dimension: [#_num_inter_level_groups][#_num_proc_in_inter_level_group[\a group_id]]
    */
    int** _inter_level_group_ranks;

    /* *****************
    * member functions *
    *******************/
    /**
    * \brief function that sets up (extended) work groups basing on the provided information
    *
    * This function is called on all processes of the process group.
    *
    * \author Hilmar Wobker
    */
    void _create_work_groups()
    {
      CONTEXT("ManagerComp::_create_work_groups()");

      /* *************************************************************************
      * send the 2D array _ext_work_group_ranks to the non-coordinator processes *
      * using a self-defined MPI datatype                                        *
      ***************************************************************************/

      // get base address
      MPI_Aint base;
      MPI_Get_address(_ext_work_group_ranks, &base);

      // calculate offsets of the subarray addresses w.r.t. base adress
      MPI_Aint* displacements = new MPI_Aint[_num_work_groups];
      for (unsigned int i(0) ; i < _num_work_groups ; ++i)
      {
        MPI_Get_address(_ext_work_group_ranks[i], &displacements[i]);
        displacements[i] -= base;
      }

      // create MPI datatype
      MPI_Datatype int_array_2d;
      MPI_Type_create_hindexed(_num_work_groups, reinterpret_cast<int*>(_num_proc_in_ext_work_group), displacements,
                               MPI_INTEGER, &int_array_2d);
      MPI_Type_commit(&int_array_2d);

      // let the coordinator send the array to all non-coordinator processes
      MPI_Bcast(_ext_work_group_ranks, 1, int_array_2d, _process_group->rank_coord(), _process_group->comm());

      // free the datatype definition again
      MPI_Type_free(&int_array_2d);
      // delete offset array
      delete [] displacements;


      /* ***********************************
      * now begin creating the work groups *
      *************************************/

      // boolean array indicating to which work groups this process belongs
      _belongs_to_group = new bool[_num_work_groups];
      for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)
      {
        // intialise with false
        _belongs_to_group[igroup] = false;
        for(unsigned int j(0) ; j < _num_proc_in_ext_work_group[igroup] ; ++j)
        {
          if(_process_group->rank() == _ext_work_group_ranks[igroup][j])
          {
            _belongs_to_group[igroup] = true;
          }
        }
      } // for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)

      // create WorkGroupExt objects including MPI groups and MPI communicators
      // It is not possible to set up all work groups in one call, since the processes building the work groups are
      // not necessarily disjunct. Hence, there are as many calls as there are work groups. All processes not belonging
      // to the work group currently created call the MPI_Comm_create() function with dummy comm. and dummy group.

      _work_groups.resize(_num_work_groups, nullptr);
      for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)
      {
        MPI_Group group_work;
        MPI_Comm group_comm;

        int mpi_error_code = MPI_Group_incl(_process_group->group(), _num_proc_in_ext_work_group[igroup],
                                            const_cast<int*>(_ext_work_group_ranks[igroup]), &group_work);
        validate_error_code_mpi(mpi_error_code, "MPI_Group_incl");

        // Create the group communicator for, among others, collective operations.
        // It is essential that *all* processes of the process group participate.
        mpi_error_code = MPI_Comm_create(_process_group->comm(), group_work, &group_comm);
        validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");

        if(belongs_to_group()[igroup])
        {
          _work_groups[igroup] = new WorkGroupExt(group_comm, igroup, (bool)_group_contains_extra_coord[igroup],
                                                  _process_group->is_coordinator());
        }
      } // for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)
    } // _create_work_groups()


  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /// CTOR
    ManagerComp(ProcessGroup* process_group)
      : _process_group(process_group),
        _rank_load_balancer(MPI_PROC_NULL),
        _num_work_groups(0),
        _num_proc_in_ext_work_group(nullptr),
        _group_contains_extra_coord(nullptr),
        _ext_work_group_ranks(nullptr),
        _belongs_to_group(nullptr)
    {
      CONTEXT("ManagerComp::ManagerComp()");
    }

    /// DTOR
    virtual ~ManagerComp()
    {
      CONTEXT("ManagerComp::~ManagerComp()");
      if (_ext_work_group_ranks != nullptr)
      {
        for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)
        {
          delete [] _ext_work_group_ranks[igroup];
        }
        delete [] _ext_work_group_ranks;
        _ext_work_group_ranks = nullptr;
      }

      if (_belongs_to_group != nullptr)
      {
        delete [] _belongs_to_group;
        _belongs_to_group = nullptr;
      }

      while (!_work_groups.empty())
      {
        // call the destructor of the last element in the vector
        delete _work_groups.back();
        // delete the pointer from the vector
        _work_groups.pop_back();
      }

      if (_num_proc_in_ext_work_group != nullptr)
      {
        delete [] _num_proc_in_ext_work_group;
      }

      if (_group_contains_extra_coord != nullptr)
      {
        delete [] _group_contains_extra_coord;
      }
    }


    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the array indicating to which work groups this process belongs
    *
    * \return array indicating to which work groups this process belongs
    */
    inline bool* belongs_to_group() const
    {
      CONTEXT("ManagerComp::belongs_to_group()");
      return _belongs_to_group;
    }



    /* *****************
    * member functions *
    *******************/
  };
} // namespace FEAST

#endif // guard KERNEL_MANAGER_COMP_HPP

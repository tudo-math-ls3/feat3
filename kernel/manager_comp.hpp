/* GENERAL_REMARK_BY_HILMAR:
 * See GENERAL_REMARK_BY_HILMAR in manager.hpp.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_MANAGER_COMP_HPP
#define KERNEL_MANAGER_COMP_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>
#include <kernel/work_group_ext.hpp>
#include <kernel/interlevel_group.hpp>

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

namespace FEAST
{
  /**
  * \brief defines manager of a compute process group (super class of coordinator and non-coordinator versions)
  *
  * See the description of class Manager.
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
    * \brief boolean array indicating to which (ext.) work groups this process belongs
    *
    * Dimension: [#_num_work_groups]
    */
    bool* _belongs_to_work_group;

    /// number of interlevel groups
    unsigned int _num_interlevel_groups;

    /**
    * \brief array of number of processes in each interlevel group
    *
    * Dimension: [#_num_interlevel_groups]
    */
    unsigned int* _num_proc_in_interlevel_group;

    /**
    * \brief array of positions of the root rank within the array _interlevel_group_ranks per interlevel group
    *
    * Dimension: [#_num_interlevel_groups]
    */
    unsigned int* _pos_of_root_in_interlevel_group;

    /**
    * \brief 2-dim. array of process group ranks building the interlevel groups
    *
    * Dimension: [#_num_interlevel_groups][#_num_proc_in_interlevel_group[\a group_id]]
    */
    int** _interlevel_group_ranks;

    /**
    * \brief boolean array indicating to which interlevel group this process belongs
    *
    * Dimension: [#_num_interlevel_groups]
    */
    bool* _belongs_to_interlevel_group;

    /// vector of interlevel groups the manager manages
    std::vector<InterlevelGroup*> _interlevel_groups;

    /* *****************
    * member functions *
    *******************/
    /**
    * \brief sets up (extended) work groups and interlevel groups basing on the provided information
    *
    * This function is called on all processes of the compute process group.
    *
    * \author Hilmar Wobker
    */
    void _create_work_groups_common()
    {
      CONTEXT("ManagerComp::_create_work_groups_common()");

      /* ************************************************************************
      * send the 2D arrays _ext_work_group_ranks and _interlevel_group_ranks to *
      * the non-coordinator processes using a self-defined MPI datatype         *
      **************************************************************************/

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


      if(_num_interlevel_groups > 0)
      {
        // now _interlevel_group_ranks
        MPI_Get_address(_interlevel_group_ranks, &base);

        // calculate offsets of the subarray addresses w.r.t. base adress
        displacements = new MPI_Aint[_num_interlevel_groups];
        for (unsigned int i(0) ; i < _num_interlevel_groups ; ++i)
        {
          MPI_Get_address(_interlevel_group_ranks[i], &displacements[i]);
          displacements[i] -= base;
        }

        // create MPI datatype
        MPI_Type_create_hindexed(_num_interlevel_groups, reinterpret_cast<int*>(_num_proc_in_interlevel_group),
                                 displacements, MPI_INTEGER, &int_array_2d);
        MPI_Type_commit(&int_array_2d);

        // let the coordinator send the array to all non-coordinator processes
        MPI_Bcast(_interlevel_group_ranks, 1, int_array_2d, _process_group->rank_coord(), _process_group->comm());

        // free the datatype definition again
        MPI_Type_free(&int_array_2d);
        // delete offset array
        delete [] displacements;
      }


      /* ***********************************
      * now begin creating the work groups *
      *************************************/

      // boolean array indicating to which work groups this process belongs
      _belongs_to_work_group = new bool[_num_work_groups];
      for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)
      {
        // intialise with false
        _belongs_to_work_group[igroup] = false;
        for(unsigned int j(0) ; j < _num_proc_in_ext_work_group[igroup] ; ++j)
        {
          if(_process_group->rank() == _ext_work_group_ranks[igroup][j])
          {
            _belongs_to_work_group[igroup] = true;
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
        MPI_Group group;
        MPI_Comm comm;

        int mpi_error_code = MPI_Group_incl(_process_group->group(), _num_proc_in_ext_work_group[igroup],
                                            const_cast<int*>(_ext_work_group_ranks[igroup]), &group);
        validate_error_code_mpi(mpi_error_code, "MPI_Group_incl");

        // Create the group communicator.
        // It is essential that *all* processes of the process group participate.
        mpi_error_code = MPI_Comm_create(_process_group->comm(), group, &comm);
        validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");

        if(belongs_to_work_group()[igroup])
        {
          _work_groups[igroup] = new WorkGroupExt(comm, igroup, _group_contains_extra_coord[igroup] != 0,
                                                  _process_group->is_coordinator());
        }
        // free the group again
        MPI_Group_free(&group);
      } // for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)


      /* *****************************************
      * now begin creating the interlevel groups *
      *******************************************/

      if(_num_interlevel_groups > 0)
      {
        // boolean array indicating to which interlevel groups this process belongs
        _belongs_to_interlevel_group = new bool[_num_interlevel_groups];
        for(unsigned int igroup(0) ; igroup < _num_interlevel_groups ; ++igroup)
        {
          // intialise with false
          _belongs_to_interlevel_group[igroup] = false;
          for(unsigned int j(0) ; j < _num_proc_in_interlevel_group[igroup] ; ++j)
          {
            if(_process_group->rank() == _interlevel_group_ranks[igroup][j])
            {
              _belongs_to_interlevel_group[igroup] = true;
            }
          }
        } // for(unsigned int igroup(0) ; igroup < _num_interlevel_groups ; ++igroup)

        // create InterlevelGroup objects including MPI groups and MPI communicators
        // Since interlevel groups are also not necessariyl disjunct, we proceed as for the work groups.

        _interlevel_groups.resize(_num_interlevel_groups, nullptr);
        for(unsigned int igroup(0) ; igroup < _num_interlevel_groups ; ++igroup)
        {
          MPI_Group group;
          MPI_Comm comm;

          int mpi_error_code = MPI_Group_incl(_process_group->group(), _num_proc_in_interlevel_group[igroup],
                                              const_cast<int*>(_interlevel_group_ranks[igroup]), &group);
          validate_error_code_mpi(mpi_error_code, "MPI_Group_incl");

          // Create the group communicator.
          // It is essential that *all* processes of the process group participate.
          mpi_error_code = MPI_Comm_create(_process_group->comm(), group, &comm);
          validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");

          if(belongs_to_interlevel_group()[igroup])
          {
            // since MPI_Group_incl(...) arranges the processes in the same order as they appear in the array
            // _interlevel_group_ranks[igroup], the local rank (w.r.t. to the newly created interlevel group) coincides
            // with the position in the array _interlevel_group_ranks[igroup]. So, we can pass the value
            // _pos_of_root_in_interlevel_group[igroup] as local root rank.
            _interlevel_groups[igroup] = new InterlevelGroup(comm, igroup, _pos_of_root_in_interlevel_group[igroup]);
          }
          // free the group again
          MPI_Group_free(&group);
        } // for(unsigned int igroup(0) ; igroup < _num_interlevel_groups ; ++igroup)
      } // if(_num_interlevel_groups > 0)
    } // _create_work_groups_common()


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
        _belongs_to_work_group(nullptr),
        _num_interlevel_groups(0),
        _num_proc_in_interlevel_group(nullptr),
        _pos_of_root_in_interlevel_group(nullptr),
        _interlevel_group_ranks(nullptr),
        _belongs_to_interlevel_group(nullptr)
    {
      CONTEXT("ManagerComp::ManagerComp()");
    }

    /// DTOR
    virtual ~ManagerComp()
    {
      CONTEXT("ManagerComp::~ManagerComp()");
      if (_num_proc_in_ext_work_group != nullptr)
      {
        delete [] _num_proc_in_ext_work_group;
      }

      if (_group_contains_extra_coord != nullptr)
      {
        delete [] _group_contains_extra_coord;
      }

      if (_ext_work_group_ranks != nullptr)
      {
        for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)
        {
          delete [] _ext_work_group_ranks[igroup];
        }
        delete [] _ext_work_group_ranks;
        _ext_work_group_ranks = nullptr;
      }

      if (_belongs_to_work_group != nullptr)
      {
        delete [] _belongs_to_work_group;
        _belongs_to_work_group = nullptr;
      }

      while (!_work_groups.empty())
      {
        // call the destructor of the last element in the vector
        delete _work_groups.back();
        // delete the pointer from the vector
        _work_groups.pop_back();
      }

      if (_num_proc_in_interlevel_group != nullptr)
      {
        delete [] _num_proc_in_interlevel_group;
      }

      if (_pos_of_root_in_interlevel_group != nullptr)
      {
        delete [] _pos_of_root_in_interlevel_group;
      }

      if (_interlevel_group_ranks != nullptr)
      {
        for(unsigned int igroup(0) ; igroup < _num_interlevel_groups ; ++igroup)
        {
          delete [] _interlevel_group_ranks[igroup];
        }
        delete [] _interlevel_group_ranks;
        _interlevel_group_ranks = nullptr;
      }

      if (_belongs_to_interlevel_group != nullptr)
      {
        delete [] _belongs_to_interlevel_group;
        _belongs_to_interlevel_group = nullptr;
      }

      while (!_interlevel_groups.empty())
      {
        // call the destructor of the last element in the vector
        delete _interlevel_groups.back();
        // delete the pointer from the vector
        _interlevel_groups.pop_back();
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
    inline bool* belongs_to_work_group() const
    {
      CONTEXT("ManagerComp::belongs_to_work_group()");
      return _belongs_to_work_group;
    }


    /**
    * \brief getter for the array indicating to which interlevel groups this process belongs
    *
    * \return array indicating to which interlevel groups this process belongs
    */
    inline bool* belongs_to_interlevel_group() const
    {
      CONTEXT("ManagerComp::belongs_to_interlevel_group()");
      return _belongs_to_interlevel_group;
    }
  }; //ManagerComp
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_MANAGER_COMP_HPP

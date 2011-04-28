#pragma once
#ifndef KERN_MANAG_COMP_NON_COORD_HPP
#define KERN_MANAG_COMP_NON_COORD_HPP 1

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
#include <kernel/manager_comp.hpp>

/// FEAST namespace
namespace FEAST
{

  /// forward declaration of class Manager (needed for friend declarations)
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class Manager;

  /**
  * \brief defines manager of a compute process group on non-coordinator processes
  *
  * \author Hilmar Wobker
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class ManagerCompNonCoord
    : public ManagerComp<space_dim_, world_dim_>
  {

    /// declare Manager function as friend to grant access to function _create_work_groups_prepare()
    friend void Manager<space_dim_, world_dim_>::create_work_groups();
    friend unsigned int Manager<space_dim_, world_dim_>::test_communication();

  private:

    /* **********
    * shortcuts *
    ************/
    /**
    * \brief getter for the vector of ext. work groups
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return reference to the vector of ext. work groups
    */
    inline const std::vector<WorkGroupExt*>& _work_groups() const
    {
      CONTEXT("ManagerCompNonCoord::work_groups()");
      return ManagerComp<space_dim_, world_dim_>::_work_groups;
    }


    /**
    * \brief getter for the number of ext. work groups
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return number of ext. work groups
    */
    inline unsigned int _num_work_groups() const
    {
      CONTEXT("ManagerCompNonCoord::_num_work_groups()");
      return ManagerComp<space_dim_, world_dim_>::_num_work_groups;
    }


    /**
    * \brief getter for the array indicating to which ext. work groups this process belongs
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return array indicating to which ext. work groups this process belongs
    */
    inline bool* _belongs_to_work_group() const
    {
      CONTEXT("ManagerCompNonCoord::_belongs_to_work_group()");
      return ManagerComp<space_dim_, world_dim_>::_belongs_to_work_group;
    }


    /**
    * \brief getter for the number of interlevel groups
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return number of interlevel groups
    */
    inline unsigned int _num_interlevel_groups() const
    {
      CONTEXT("ManagerCompNonCoord::_num_interlevel_groups()");
      return ManagerComp<space_dim_, world_dim_>::_num_interlevel_groups;
    }


    /**
    * \brief getter for the array indicating to which interlevel groups this process belongs
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return array indicating to which interlevel groups this process belongs
    */
    inline bool* _belongs_to_interlevel_group() const
    {
      CONTEXT("ManagerCompNonCoord::_belongs_to_interlevel_group()");
      return ManagerComp<space_dim_, world_dim_>::_belongs_to_interlevel_group;
    }


    /**
    * \brief getter for the process group
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return ProcessGroup pointer #_process_group
    */
    inline ProcessGroup* _process_group() const
    {
      CONTEXT("ManagerCompNonCoord::_process_group()");
      return ManagerComp<space_dim_, world_dim_>::_process_group;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief function that sets up (extended) work groups and interlevel groups basing on the provided information
    *
    * This function is called on all non-coordinators of the manager's compute process group.
    *
    * \author Hilmar Wobker
    */
    void _create_work_groups()
    {
      CONTEXT("ManagerCompNonCoord::_create_work_groups()");

      // now the non-coordinator processes receive data which is broadcasted by the coordinator:
      //   - _num_work_groups
      //   - _num_proc_in_ext_work_group
      //   - group_contains_extra_coord

      unsigned int temp;
      int mpi_error_code = MPI_Bcast(&temp, 1, MPI_UNSIGNED, _process_group()->rank_coord(), _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      ManagerComp<space_dim_, world_dim_>::_num_work_groups = temp;
      ASSERT(_num_work_groups() > 0, "Number of work groups must be greater than 0.");

      // allocate array
      unsigned int* num_proc_in_ext_wg = new unsigned int[_num_work_groups()];
      mpi_error_code = MPI_Bcast(num_proc_in_ext_wg, _num_work_groups(), MPI_UNSIGNED, _process_group()->rank_coord(),
                                 _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      ManagerComp<space_dim_, world_dim_>::_num_proc_in_ext_work_group = num_proc_in_ext_wg;

      // allocate array
      unsigned char* group_contains_extra_coord = new unsigned char[_num_work_groups()];
      mpi_error_code = MPI_Bcast(group_contains_extra_coord, _num_work_groups(), MPI_UNSIGNED_CHAR,
                                 _process_group()->rank_coord(), _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      ManagerComp<space_dim_, world_dim_>::_group_contains_extra_coord = group_contains_extra_coord;

      // let the non-coordinator processes allocate the array _ext_work_group_ranks for rank partitioning
      int** wg_ranks = new int*[_num_work_groups()];
      for (unsigned int i(0) ; i < _num_work_groups() ; ++i)
      {
        wg_ranks[i] = new int[ManagerComp<space_dim_, world_dim_>::_num_proc_in_ext_work_group[i]];
      }
      ManagerComp<space_dim_, world_dim_>::_ext_work_group_ranks = wg_ranks;

      // now the non-coordinator processes receive further data which is broadcasted by the coordinator:
      //   - _num_interlevel_groups
      //   - _num_proc_in_interlevel_group
      //   - _pos_of_root_in_interlevel_group

      mpi_error_code = MPI_Bcast(&temp, 1, MPI_UNSIGNED, _process_group()->rank_coord(), _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      ManagerComp<space_dim_, world_dim_>::_num_interlevel_groups = temp;

      if(_num_interlevel_groups() > 0)
      {
        // let the non-coordinator processes create the array _num_proc_in_interlevel_group and receive the data
        unsigned int* num_proc_in_interlev_gr = new unsigned int[_num_interlevel_groups()];
        mpi_error_code = MPI_Bcast(num_proc_in_interlev_gr, _num_interlevel_groups(), MPI_UNSIGNED,
                                   _process_group()->rank_coord(), _process_group()->comm());
        validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
        ManagerComp<space_dim_, world_dim_>::_num_proc_in_interlevel_group = num_proc_in_interlev_gr;

        // let the non-coordinator processes create the array _pos_of_root_in_interlevel_group and receive the data
        unsigned int* pos_of_root_in_interlev_gr = new unsigned int[_num_interlevel_groups()];
        mpi_error_code = MPI_Bcast(pos_of_root_in_interlev_gr, _num_interlevel_groups(), MPI_UNSIGNED,
                                   _process_group()->rank_coord(), _process_group()->comm());
        validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
        ManagerComp<space_dim_, world_dim_>::_pos_of_root_in_interlevel_group = pos_of_root_in_interlev_gr;

        // let the non-coordinator processes allocate the array _interlevel_group_ranks for rank partitioning
        int** interlev_gr_ranks = new int*[_num_interlevel_groups()];
        for (unsigned int i(0) ; i < _num_interlevel_groups() ; ++i)
        {
          interlev_gr_ranks[i] = new int[ManagerComp<space_dim_, world_dim_>::_num_proc_in_interlevel_group[i]];
        }
        ManagerComp<space_dim_, world_dim_>::_interlevel_group_ranks = interlev_gr_ranks;
      }

      // call routine for creating work groups (within this routine the arrays _ext_work_group_ranks and
      // _interlevel_group_ranks are transferred)
      ManagerComp<space_dim_, world_dim_>::_create_work_groups_common();
    } // _create_work_groups()


    /**
    * \brief receives the relevant parts of the global graph and sets distributed graphs in the work groups
    *
    * The local graphs tell the workers with which workers of their work group they have to communicate.
    * This function must be called on all non-coordinator processes of the process group, at the same time the
    * coordinator process must call the function transfer_graphs_to_workers(). Before this function can be used, the
    * function create_work_groups() must have been called.
    *
    * \sa transfer_graphs_to_workers(Graph**)
    *
    * \author Hilmar Wobker
    */
    void _receive_and_set_graphs()
    {
      CONTEXT("ManagerCompNonCoord::_receive_and_set_graphs()");

      for(unsigned int igroup(0) ; igroup < _num_work_groups() ; ++igroup)
      {
        if(_belongs_to_work_group()[igroup])
        {
          ASSERT(!_work_groups()[igroup]->is_coordinator(), "Routine must be called on non-coordinator processes.");
          // number of neighbours of the graph node corresponding to this process (use unsigned int datatype here
          // instead of index_glob_t since MPI routines expect it)
          unsigned int num_neighbours_local;
          index_glob_t* neighbours_local(nullptr);
          int rank_coord = _work_groups()[igroup]->rank_coord();

          // receive the number of neighbours from the coordinator process
          MPI_Scatter(nullptr, 0, MPI_DATATYPE_NULL, &num_neighbours_local, 1, MPI_INTEGER,
                      rank_coord, _work_groups()[igroup]->comm());

          // receive the neighbours
          neighbours_local = new index_glob_t[num_neighbours_local];
          MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DATATYPE_NULL, neighbours_local,
                       num_neighbours_local, MPIType<index_glob_t>::value(), rank_coord,
                       _work_groups()[igroup]->comm());

          // now create distributed graph structure within the compute work groups
          _work_groups()[igroup]->work_group()->set_graph_distributed(num_neighbours_local, neighbours_local);
          // The array neighbours_local is copied inside the constructor of the distributed graph object, hence it
          // can be deallocated again.
          delete [] neighbours_local;
        } // if(_belongs_to_work_group()[igroup])
      } // for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)

    } // _receive_and_set_graphs()


    /**
    * \brief dummy function for testing work group and interlevel group communication
    *
    * The function simply performs an exchange of integer arrays between the corresponding processes.
    * Note that this function and the corresponding version in ManagerCompCoord are slightly different, which is
    * why they are not implemented in ManagerComp.
    *
    * \return flag whether test was succesful (0: succesful, >0: failed)
    */
    bool _test_communication()
    {
      CONTEXT("ManagerCompNonCoord::_test_communication()");

      unsigned int work_group_test_ok = 0;
      unsigned int interlevel_group_test_ok = 0;

      // test local neighbourhood communication
      for(unsigned int igroup(0) ; igroup < _num_work_groups() ; ++igroup)
      {
        if (_belongs_to_work_group()[igroup])
        {
          work_group_test_ok = work_group_test_ok + _work_groups()[igroup]->work_group()->test_communication();
        }
      }

      // test communication within interlevel groups
      for(unsigned int igroup(0) ; igroup < _num_interlevel_groups() ; ++igroup)
      {
        if (_belongs_to_interlevel_group()[igroup])
        {
          interlevel_group_test_ok = interlevel_group_test_ok
            + ManagerComp<space_dim_, world_dim_>::_interlevel_groups[igroup]->test_communication();
        }
      }
      return work_group_test_ok + interlevel_group_test_ok;
    } // _test_communication



  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /// CTOR
    ManagerCompNonCoord(ProcessGroup* process_group)
      : ManagerComp<space_dim_, world_dim_>(process_group)
    {
      CONTEXT("ManagerCompNonCoord::ManagerCompNonCoord()");
    }

    /// DTOR
    ~ManagerCompNonCoord()
    {
      CONTEXT("ManagerCompNonCoord::~ManagerCompNonCoord()");
    }
  }; // class ManagerCompNonCoord
} // namespace FEAST

#endif // guard KERN_MANAG_COMP_NON_COORD_HPP

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

  private:

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



    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the vector of ext. work groups
    *
    * \return reference to the vector of ext. work groups
    */
    inline const std::vector<WorkGroupExt*>& work_groups() const
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
    inline unsigned int num_work_groups() const
    {
      CONTEXT("ManagerCompNonCoord::num_work_groups()");
      return ManagerComp<space_dim_, world_dim_>::_num_work_groups;
    }


    /**
    * \brief setter for the number of ext. work groups
    *
    * \param[in] number of ext. work groups
    */
    inline void set_num_work_groups(unsigned int num_work_groups)
    {
      CONTEXT("ManagerCompNonCoord::set_num_work_groups()");
      ManagerComp<space_dim_, world_dim_>::_num_work_groups = num_work_groups;
    }


    /**
    * \brief getter for the array of number of processes per ext. work group
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return pointer to array of number of processes per ext. work group
    */
    inline unsigned int* num_proc_in_work_group() const
    {
      CONTEXT("ManagerCompNonCoord::num_proc_in_work_group()");
      return ManagerComp<space_dim_, world_dim_>::_num_proc_in_work_group;
    }


    /**
    * \brief setter for the array of number of processes per ext. work group
    *
    * \param[in] pointer to array of number of processes per ext. work group
    */
    inline void set_num_proc_in_work_group(unsigned int* num_proc_in_work_group)
    {
      CONTEXT("ManagerCompNonCoord::set_num_proc_in_work_group()");
      ManagerComp<space_dim_, world_dim_>::_num_proc_in_work_group = num_proc_in_work_group;
    }


    /**
    * \brief getter for array indicating whether the ext. work groups contain an extra process for the coordinator
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return pointer to array indicating whether the ext. work groups contain an extra process for the coordinator
    */
    inline unsigned char* group_contains_extra_coord() const
    {
      CONTEXT("ManagerCompCoord::num_proc_in_work_group()");
      return ManagerComp<space_dim_, world_dim_>::_group_contains_extra_coord;
    }


    /**
    * \brief setter for the array indicating whether the ext. work groups contain an extra process for the coordinator
    *
    * \param[in] pointer to array indicating whether the ext. work groups contain an extra process for the coordinator
    */
    inline void set_group_contains_extra_coord(unsigned char* group_contains_extra_coord)
    {
      CONTEXT("ManagerCompCoord::set_group_contains_extra_coord()");
      ManagerComp<space_dim_, world_dim_>::_group_contains_extra_coord = group_contains_extra_coord;
    }


    /**
    * \brief getter for the 2D array of ext. work group ranks
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return pointer to 2D array of ext. work group ranks
    */
    inline int** work_group_ranks() const
    {
      CONTEXT("ManagerCompNonCoord::work_group_ranks()");
      return ManagerComp<space_dim_, world_dim_>::_work_group_ranks;
    }


    /**
    * \brief setter for the 2D array of ext. work group ranks
    *
    * \param[in] pointer to 2D array of ext. work group ranks
    */
    inline void set_work_group_ranks(int** work_group_ranks)
    {
      CONTEXT("ManagerCompNonCoord::set_work_group_ranks()");
      ManagerComp<space_dim_, world_dim_>::_work_group_ranks = work_group_ranks;
    }


    /**
    * \brief getter for the array indicating to which ext. work groups this process belongs
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return array indicating to which ext. work groups this process belongs
    */
    inline bool* belongs_to_group() const
    {
      CONTEXT("ManagerCompNonCoord::belongs_to_group()");
      return ManagerComp<space_dim_, world_dim_>::_belongs_to_group;
    }


    /**
    * \brief getter for the process group
    *
    * \return ProcessGroup pointer #_process_group
    */
    inline ProcessGroup* process_group() const
    {
      CONTEXT("ManagerCompNonCoord::process_group()");
      return ManagerComp<space_dim_, world_dim_>::_process_group;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief function that sets up ext. work groups and ext. work groups basing on the provided information
    *
    * This function is called on all processes of the process group.
    *
    * \author Hilmar Wobker
    */
    void create_work_groups()
    {
      CONTEXT("ManagerCompNonCoord::create_work_groups()");

      // now the non-coordinator processes receive data which is broadcasted by the coordinator:
      //   - _num_work_groups
      //   - _num_proc_in_work_group
      //   - group_contains_extra_coord

      unsigned int temp;
      int mpi_error_code = MPI_Bcast(&temp, 1, MPI_UNSIGNED, process_group()->rank_coord(), process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      set_num_work_groups(temp);
      ASSERT(num_work_groups() > 0, "Number of work groups must be greater than 0.");

      // allocate array
      unsigned int* num_proc_in_subgr = new unsigned int[num_work_groups()];
      mpi_error_code = MPI_Bcast(num_proc_in_subgr, num_work_groups(), MPI_UNSIGNED, process_group()->rank_coord(),
                                 process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      set_num_proc_in_work_group(num_proc_in_subgr);

      // allocate array
      unsigned char* group_contains_extra_coord = new unsigned char[num_work_groups()];
      mpi_error_code = MPI_Bcast(group_contains_extra_coord, num_work_groups(), MPI_UNSIGNED_CHAR,
                                 process_group()->rank_coord(), process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      set_group_contains_extra_coord(group_contains_extra_coord);

      // let the non-coordinator processes allocate the array _work_group_ranks for rank partitioning
      int** subgr_ranks = new int*[num_work_groups()];
      for (unsigned int i(0) ; i < num_work_groups() ; ++i)
      {
        subgr_ranks[i] = new int[num_proc_in_work_group()[i]];
      }
      set_work_group_ranks(subgr_ranks);

      // call routine for creating work groups (within this routine the array work_group_ranks is transferred)
      ManagerComp<space_dim_, world_dim_>::_create_work_groups();
    } // create_work_groups()


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
    void receive_and_set_graphs()
    {
      CONTEXT("ManagerCompNonCoord::receive_and_set_graphs()");

      for(unsigned int igroup(0) ; igroup < num_work_groups() ; ++igroup)
      {
        if(belongs_to_group()[igroup])
        {
          ASSERT(!work_groups()[igroup]->is_coordinator(), "Routine must be called on non-coordinator processes.");
          // number of neighbours of the graph node corresponding to this process (use unsigned int datatype here
          // instead of index_glob_t since MPI routines expect it)
          unsigned int num_neighbours_local;
          index_glob_t* neighbours_local(nullptr);
          int rank_coord = work_groups()[igroup]->rank_coord();

          // receive the number of neighbours from the coordinator process
          MPI_Scatter(nullptr, 0, MPI_DATATYPE_NULL, &num_neighbours_local, 1, MPI_INTEGER,
                      rank_coord, work_groups()[igroup]->comm());

          // receive the neighbours
          neighbours_local = new index_glob_t[num_neighbours_local];
          MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DATATYPE_NULL, neighbours_local,
                       num_neighbours_local, MPIType<index_glob_t>::value(), rank_coord,
                       work_groups()[igroup]->comm());

          // now create distributed graph structure within the compute work groups
          work_groups()[igroup]->work_group()->set_graph_distributed(num_neighbours_local, neighbours_local);
          // The array neighbours_local is copied inside the constructor of the distributed graph object, hence it
          // can be deallocated again.
          delete [] neighbours_local;
        } // if(belongs_to_group()[igroup])
      } // for(unsigned int igroup(0) ; igroup < _num_work_groups ; ++igroup)

// TODO: remove this test code here and call it from somewhere else
      // test local neighbourhood communication
      for(unsigned int igroup(0) ; igroup < num_work_groups() ; ++igroup)
      {
        if (belongs_to_group()[igroup])
        {
          work_groups()[igroup]->work_group()->do_exchange();
        }
      }
    } // receive_and_set_graphs()
  }; // class ManagerCompNonCoord
} // namespace FEAST

#endif // guard KERN_MANAG_COMP_NON_COORD_HPP

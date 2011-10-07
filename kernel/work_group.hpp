/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR below.
 * The link to interlevel groups is still missing and only rudimentarily prepared. Actually, I'm not really sure if
 * such a link is necessary here at all. It might turn out, that you will never have to access interlevel groups from
 * within a work group. Maybe everything will be managed from outside (i.e., in some class which knows the work groups
 * and the interlevel groups). At this stage I simply don't know yet, how and where exactly the work group / interlevel
 * group functionality is used, so I thought it is better to *not* fix this link already.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_WORK_GROUP_HPP
#define KERNEL_WORK_GROUP_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/logger.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>
//#include <kernel/interlevel_group.hpp>
#include <kernel/graph.hpp>

// includes, system
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <mpi.h>
#include <math.h>

namespace FEAST
{
  /**
  * \brief describes a work group, i.e. a set of worker processes sharing the same MPI communicator and performing
  *        the same compute task
  *
  * A WorkGroup object is part of a ProcessSubgroup. It describes a set of compute processes sharing the MPI
  * communicator. The member of the work group living on this process is called 'worker' (which is not realised via an
  * extra class). Communication between different work groups is done via the enclosing ProcessGroup communicator.
  * For an example concerning WorkGroup and ProcessSubgroup creation, see the description of class ProcessSubgroup.
  *
  * \author Hilmar Wobker
  *
  */
  class WorkGroup
    : public ProcessGroup
  {

  private:

    /* *****************
    * member variables *
    *******************/
//    /**
//    * \brief additional communicator representing an optimised process topology
//    *
//    * This communicator results from a call to MPI_Dist_graph_create(...), which remaps the processes of this work group
//    * to optimise communication patterns.
//    */
//   COMMENT_HILMAR: Everything that has to do with comm_opt is deactivated for the time being (see comment in routine
//      set_graph_distributed.
//    MPI_Comm _comm_opt;

    /**
    * \brief local portions of a distributed graph
    *
    * This object is sent to this process by the coordinator of the parent process group. It contains the parent process
    * group ranks of the neighbours of this worker. It can be used to create an MPI topology graph via
    * MPI_Dist_graph_create().
    */
    GraphDistributed* _graph_distributed;

    /// interlevel group this process builds with the processes of the finer grid work group
//    InterlevelGroup* _group_finer;

    /// interlevel group this process builds with the process of the coarser grid work group
//    InterlevelGroup* _group_coarser;

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
    */
    WorkGroup(
      MPI_Comm comm,
      unsigned int const group_id)
      : ProcessGroup(comm, group_id),
//        _comm_opt(MPI_COMM_NULL),
        _graph_distributed(nullptr)
    {
      CONTEXT("WorkGroup::WorkGroup()");

      // debugging output
      // write some individual messages to file
      String s("I have COMM_WORLD rank " + stringify(Process::rank) + " and group rank " + stringify(_rank)
                    + " in work group " + stringify(_group_id) + ".");
      if(is_coordinator())
      {
        s += " I am the coordinator!";
      }
      s += "\n";
      Logger::log(s);
      // end of debugging output
    } // constructor


    /// DTOR (automatically virtual since DTOR of base class ProcessGroup is virtual)
    ~WorkGroup()
    {
      CONTEXT("WorkGroup::~WorkGroup()");
      if(_graph_distributed != nullptr)
      {
        delete _graph_distributed;
        _graph_distributed = nullptr;
      }
    }


    /* ******************
    * getters & setters *
    ********************/
//  /**
//  * \brief getter for the optimised MPI communicator
//  *
//  * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
//  * cannot change the object.
//  *
//  * \return reference to MPI_Comm #_comm
//  */
//  inline const MPI_Comm& comm_opt() const
//  {
//    return _comm_opt;
//  }


    /**
    * \brief getter for the distributed graph object
    *
    * \return pointer to GraphDistributed #_graph_distributed
    */
    inline GraphDistributed* graph_distributed() const
    {
      CONTEXT("WorkGroup::graph_distributed()");
      return _graph_distributed;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief sets local graph portion of a distributed graph corresponding to one node (=process)
    *
    * \param[in] num_neighbours
    * number of neighbours of one graph node
    *
    * \param[in] neighbours
    * neighbours of this portion of the graph
    *
    * \note Currently, we do not distinguish face-, edge- and vertex-neighbours. The graph structure
    *       simply contains all neighbours.
    */
    void set_graph_distributed(
      Index const num_neighbours,
      Index * neighbours)
    {
      CONTEXT("WorkGroup::set_graph_distributed()");
      _graph_distributed = new GraphDistributed(num_neighbours, neighbours);

      // log into process-own log file
      Logger::log(_graph_distributed->print());

// COMMENT_HILMAR, 14.10.2010: The plan was to use MPI_Dist_graph_create(...) to create the MPI graph topology.
//   Unfortunately, this function has only been introduced with MPI-2.2 and is not yet implemented in OpenMPI yet.
//   The only alternative that already exists, MPI_Graph_create(...), requires that all processes know the complete
//   graph (and not only the coordinator process). IMHO, it doesn't make sense to further pursue this MPI process
//   topology topic in the moment. Let us stick with the standard MPI communicators and let us see what the ITMC
//   will say.
// COMMENT_HILMAR: By now, the weights and info objects are empty. Does it make sense to, e.g., use different weights
//   for edge neighbours on the one hand and diagonal neighbours on the other hand? Should diagonal neighbours appear
//   in the graph topology at all? Or should we only take the edge neighbours here?
//
//    unsigned int sources[1];
//    unsigned int degrees[1];
//    sources[0] = _rank;
//    degrees[0] = _graph_distributed->num_neighbours();
//    MPI_Dist_graph_create(_comm, 1, sources, degrees, _graph_distributed->neighbours(), nullptr,
//                          MPI_INFO_NULL, true, _comm_opt());
    }


    /**
    * \brief dummy function for testing local neighbour communication via non-blocking sends/receives
    *
    * The function simply performs an exchange of integer arrays between neighbouring processes. The aim is to
    * demonstrate that this is quite easy when using non-blocking send and receive calls.
    *
    * COMMENT_HILMAR: Actually, I intended to make this function private and declare the functions
    *     ManagerCompCoord<space_dim_, world_dim_>::_test_communication()
    *   and
    *     ManagerCompNonCoord<space_dim_, world_dim_>::_test_communication()
    *   as friends. But since the WorkGroup class doesn't know the template parameters space_dim_ and world_dim_,
    *   this is not possible. Maybe there is a way to realise it, but I didn't have the time to find out how.
    *
    * \return flag whether test was succesful (0: succesful, >0: failed)
    */
    unsigned int test_communication()
    {
      CONTEXT("WorkGroup::test_communication()");

      // flag whether test is succesful
      bool test_result = 0;

      // length of the integer arrays to be exchanged
      Index const n = 10;

      Index num_neighbours = _graph_distributed->num_neighbours();
      Index* neighbours = _graph_distributed->neighbours();

      // arrays for sending and receiving
      // (The MPI standard says that the send buffer given to MPI_Isend(...) should neither be overwritten nor read(!)
      // until the communication has been completed. Hence, we cannot send the same array to all neighbours, but we
      // have to send num_neighbours different arrays.)
      unsigned long** a = new unsigned long*[num_neighbours];
      unsigned long** a_recv = new unsigned long*[num_neighbours];

      // fill the array to be sent
      for(Index i(0) ; i < num_neighbours; ++i)
      {
        a[i] = new unsigned long[n];
        a_recv[i] = new unsigned long[n];
        for(Index j(0) ; j < n; ++j)
        {
          a[i][j] = 100000*(unsigned long)(_rank) + 100*(unsigned long)(neighbours[i]) + (unsigned long)(j);
        }
      }

      // debugging output
      String s;
      for(Index i(0) ; i < num_neighbours; ++i)
      {
        s = stringify(a[i][0]);
        for(Index j(1) ; j < n; ++j)
        {
          s +=  " " + stringify(a[i][j]);
        }
        Logger::log("Work group " + stringify(_group_id) + ": Process " + stringify(_rank) + " sends [" + s
                    + "] to neighbour " + stringify(neighbours[i]) + ".\n");
      }

      // request and status objects necessary for communication
      MPI_Request* requests = new MPI_Request[num_neighbours];
      MPI_Status* statuses = new MPI_Status[num_neighbours];

      // post sends to all neighbours
      for(Index i(0) ; i < num_neighbours; ++i)
      {
        MPI_Isend(a[i], n, MPI_UNSIGNED_LONG, int(neighbours[i]), 0, _comm, &requests[i]);
        // request objects are not needed
        MPI_Request_free(&requests[i]);
      }
      // post receives from all neighbours
      for(Index i(0) ; i < num_neighbours; ++i)
      {
        MPI_Irecv(a_recv[i], n, MPI_UNSIGNED_LONG, int(neighbours[i]), 0, _comm, &requests[i]);
      }
      // wait for all receives to complete
      MPI_Waitall(int(num_neighbours), requests, statuses);
      delete [] requests;
      delete [] statuses;

      // debugging output
      for(Index i(0) ; i < num_neighbours; ++i)
      {
        s = stringify(a_recv[i][0]);
        for(Index j(1) ; j < n; ++j)
        {
          s += " " + stringify(a_recv[i][j]);
          // We cannot really do a clever check here, so we simply check whether a_recv really contains the value we
          // expect. (The general problem is: When there is something wrong with the communication, then usually MPI
          // crashes completely. On the other hand, when communication is fine then usually correct values are sent.
          // So, it is very unlikely that communication works *AND* this test here returns false.)
          unsigned long expected_value = 100000*(unsigned long)(neighbours[i]) + 100*(unsigned long)(_rank)
            + (unsigned long)(j);
          if(expected_value != a_recv[i][j])
          {
            test_result = 1;
          }
        }
        Logger::log("Work group " + stringify(_group_id) + ": Process " + stringify(_rank) + " received [" + s
                    + "] from neighbour " + stringify(neighbours[i]) + ".\n");
      }
      for(unsigned int i(0) ; i < num_neighbours; ++i)
      {
        delete [] a[i];
        delete [] a_recv[i];
      }
      delete [] a;
      delete [] a_recv;

      return test_result;
    } // test_communication()
  };
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_WORK_GROUP_HPP

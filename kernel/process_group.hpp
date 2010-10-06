#pragma once
#ifndef KERNEL_PROCESS_GROUP_HPP
#define KERNEL_PROCESS_GROUP_HPP 1

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
#include <kernel/worker.hpp>

/**
* \brief group of processes sharing an MPI communicator
*
* \author Hilmar Wobker
*/
class ProcessGroup
{

private:

protected:
  /* *****************
  * member variables *
  *******************/
  /// MPI group representing the processes of the group
  MPI_Group _group;

  /// communicator shared by all processes of the group
  MPI_Comm _comm;

  /// number of processes in this group
  int _num_processes;

  /// rank of this process with respect to the ProcessGroup's communicator
  int _rank;

  /**
  * \brief array of ranks the processes of this group have in the parent group
  *
  * Dimension: [#_num_processes]
  */
  // COMMENT_HILMAR, 15.9.2010:
  // Currently, this is simply a pointer to the part of the corresponding array
  // _group_ranks_world[my_group] in universe.hpp. Maybe it makes sense to copy the array and delete the one in
  // universe.hpp.
  int* _ranks_group_parent;

  /// pointer to the process group from which this process group has been spawned
  ProcessGroup* _process_group_parent;

//    /// process groups spawned from this process group
// COMMENT_HILMAR, 15.9.2010: not sure yet if this is needed
//    ProcessGroup** _process_group_child;

  /// id of the process group
  int const _group_id;

  /// rank of the coordinator process within the group
  int const _rank_coord;

// COMMENT_HILMAR: every process group will certainly need its own MPI buffer... then activate this code.
//  /// buffer size
//  static int BUFFERSIZE_BYTES;
//
//  /// buffer for MPI communication
//  char* _buffer;


public:

  /* *************************
  * constructor & destructor *
  ***************************/
  /**
  * \brief constructor for the case the MPI_Group and the MPI_Comm already exist
  *
  * This constructor is intended for creating a process group object containing all COMM_WORLD processes. The
  * coordinator is set to rank 0.
  *
  * \param[in] comm
  * communicator shared by the group processes
  *
  * \param[in] num_processes
  * number of processes in the group
  */
  ProcessGroup(
    MPI_Comm comm,
    int num_processes)
    : _comm(comm),
      _num_processes(num_processes),
      _ranks_group_parent(nullptr),
      _process_group_parent(nullptr),
      _group_id(-1),
      _rank_coord(0)
  {
    // get MPI_Group object for the given communicator
    int mpi_error_code = MPI_Comm_group(_comm, &_group);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_group");

    // and finally look up the local rank of this process w.r.t. the group's communicator
    mpi_error_code = MPI_Group_rank(_group, &_rank);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_rank");

    // since this is the world group of processes, the local and the global rank should be equal
    assert(Process::rank == _rank);



// COMMENT_HILMAR: every process group will certainly need its own MPI buffer... then activate this code.
//    // in the world group buffers are not needed
//    _buffer = nullptr;
  }

  /**
  * \brief constructor for the case the MPI_Group and the corresponding communicator have to be created
  *
  * This constructor can be used for splitting the complete set of COMM_WORLD processes into subgroups for performing
  * two completely seperated tasks (e.g., for multiphysics).
  * The coordinator of this process group is set to rank 0.
  */
  ProcessGroup(
    int num_processes,
    int ranks_group_parent[],
    ProcessGroup* process_group_parent,
    int const group_id)
    : _num_processes(num_processes),
      _ranks_group_parent(ranks_group_parent),
      _process_group_parent(process_group_parent),
      _group_id(group_id),
      _rank_coord(0)
  {
    int mpi_error_code = MPI_Group_incl(_process_group_parent->_group, _num_processes,
                                        _ranks_group_parent, &_group);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_incl");

    // Create the group communicator for, among others, collective operations.
    // It is essential that *all* processes in MPI_COMM_WORLD participate since MPI_COMM_WORLD
    // is a parent sub-universe, i.e., something to spawn from
    mpi_error_code = MPI_Comm_create(_process_group_parent->_comm, _group, &_comm);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");

    // and finally look up the local rank of this process w.r.t. the group's communicator
    mpi_error_code = MPI_Group_rank(_group, &_rank);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_rank");

// COMMENT_HILMAR: every process group will certainly need its own MPI buffer... then activate this code.
//    // allocate communication buffers
//    int size_of_char;
//    mpi_error_code = MPI_Pack_size(1, MPI_CHAR, MPI_COMM_WORLD, &size_of_char);
//    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack_size");
//
//    if (size_of_char != 1)
//    {
//      // In the unlikely case that a char is not 1 byte, determine the size of the buffer. The -0.001 is a precaution
//      // for the case the division results in 42.00000000001 instead of 42.0.
//      int buffer_size = ceil(BUFFERSIZE_BYTES/size_of_char - 0.001);
// //    std::cout << "buffer size: " << StringUtils::stringify(buffer_size) <<  std::endl;
//      _buffer = new char[buffer_size];
//    }
//    else
//    {
//      // otherwise, array size (in bytes) and number of array entries are equal
//      _buffer = new char[BUFFERSIZE_BYTES];
//    }

  }

  /* ***********
  * destructor *
  *************/
  /// destructor
  ~ProcessGroup()
  {
// COMMENT_HILMAR: every process group will certainly need its own MPI buffer... then activate this code.
//    if (_buffer != nullptr)
//    {
//      delete [] _buffer;
//    }
  }
  /* ******************
  * getters & setters *
  ********************/
  /**
  * \brief getter for the MPI group
  *
  * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
  * cannot change the object.
  *
  * \return reference to MPI_Group #_group
  */
  inline const MPI_Group& group() const
  {
    return _group;
  }

  /**
  * \brief getter for the MPI communicator
  *
  * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
  * cannot change the object.
  *
  * \return reference to MPI_Comm #_comm
  */
  inline const MPI_Comm& comm() const
  {
    return _comm;
  }

  /**
  * \brief getter for the number of processes
  *
  * \return number of processes #_num_processes
  */
  inline int num_processes() const
  {
    return _num_processes;
  }

  /**
  * \brief getter for the pointer to the parent process group
  *
  * \return parent process group pointer #_process_group_parent
  */
  inline ProcessGroup* process_group_parent() const
  {
    return _process_group_parent;
  }

  /**
  * \brief getter for group ID
  *
  * \return group ID #_group_id
  */
  inline int group_id() const
  {
    return _group_id;
  }

  /**
  * \brief getter for the rank w.r.t. the group communicator this process belongs to
  *
  * \return rank #_rank w.r.t. process group communicator
  */
  inline int rank() const
  {
    return _rank;
  }

  /* *****************
  * member functions *
  *******************/
  /**
  * \brief checks whether this process is the coordinator of the process group
  *
  * \return true if this process is the coordinator of the process group, otherwise false
  */
  inline bool is_coordinator() const
  {
    return _rank == _rank_coord;
  }

  /**
  * \brief sending individual messages to the master
  *
  * With this function the processes of the process group can send individual log messages to the master. The
  * coordinator of the process group collects all these log messages in the order of process group ranks and sends
  * them as one MPI message to the master. Depending on the target the master then displays these messages on the
  * screen or writes them to the log file (one line per message). This function has to be called by \em all processes
  * of the process group! Note no further information will be appended to the messages (like rank of the process from
  * which the messages comes). Such infomration has to be coded directly into the message string.
  *
  * \param[in] message
  * string representing the log message
  *
  * \param[in] target
  * output target Logger::SCREEN, Logger::FILE or Logger::SCREEN_FILE (default: Logger::SCREEN_FILE)
  *
  * \sa Logger::log_master_array()
  *
  * \author Hilmar Wobker
  */
  void log_indiv_master(std::string message, Logger::target target = Logger::SCREEN_FILE)
  {

    // add 1 to the message length since string::c_str() adds the null termination symbol to the resulting char array
    int length = message.length() + 1;

    if (!is_coordinator())
    {
      /* ***************************************
      * code for all non-coordinator processes *
      *****************************************/

      // coordinator process gathers the lengths of the messages
      MPI_Gather(&length, 1, MPI_INT, nullptr, 0, MPI_DATATYPE_NULL, _rank_coord, _comm);

      // coordinator process gathers the messages
      // (Here, it is necessary to cast away the const'ness of string::c_str() since the MPI routine expects a
      // non-const send buffer.)
      MPI_Gatherv(const_cast<char *>(message.c_str()), length, MPI_CHAR, nullptr, nullptr, nullptr,
                  MPI_DATATYPE_NULL, _rank_coord, _comm);
// COMMENT_HILMAR: not sure, if it is actually legitimate to simply cast the const away or if one has to manually copy
// the string like this:
//    char bla[length];
//    strcpy(bla, message.c_str());
//    MPI_Gatherv(bla, length, MPI_CHAR, nullptr, nullptr, nullptr, MPI_DATATYPE_NULL, _rank_coord, _comm);
// With gcc 4.4.0 and intel 11.1 it works for tiny test problems. (If the function call is changed, then also
// change the call on coordinator side some lines below.)
    }
    else
    {
      /* ********************************
      * code on the coordinator process *
      **********************************/

      // receive buffer for storing the lengths of the messages
      int msg_lengths[_num_processes];
      // array for storing the start positions of the single messages in the receive buffer
      int msg_start_pos[_num_processes];

      // gather the lengths of the messages from the other processes
      MPI_Gather(&length, 1, MPI_INT, msg_lengths, 1, MPI_INT, _rank_coord, _comm);

      // set start positions of the single messages in the receive buffer
      msg_start_pos[0] = 0;
      for(int i(1) ; i < _num_processes ; ++i)
      {
        msg_start_pos[i] = msg_start_pos[i-1] + msg_lengths[i-1];
      }

      // determine total length of the messages and allocate receive buffer accordingly
      int total_length(msg_start_pos[_num_processes-1] + msg_lengths[_num_processes-1]);

      // receive buffer for (consecutively) storing the messages
      char messages[total_length];

//      // debug output
//      for(int i(0) ; i < _num_processes ; ++i)
//      {
//        std::cout << "local process " << StringUtils::stringify(i) << ": length of string = "
//                  << StringUtils::stringify(msg_lengths[i]) << std::endl;
//      }
//      std::cout << "Total length: " << StringUtils::stringify(total_length) << std::endl;

      // gather the messages from the other processes.
      MPI_Gatherv(const_cast<char *>(message.c_str()), length, MPI_CHAR, messages, msg_lengths, msg_start_pos,
                  MPI_CHAR, 0, _comm);

//      // debug output
//      // output the strings collected from all processes by using corresponding offsets in the
//      // char array (pointer arithmetic)
//      for(int i(0) ; i < _num_processes ; ++i)
//      {
//        std::cout << "Process " << Process::rank << " writes: "
//                  << std::string(messages + msg_start_pos[i], msg_lengths[i]) << std::endl;
//      }

      // now send everything to the master
      Logger::log_master_array(_num_processes, msg_lengths, total_length, messages, target);

    }
  } // log_indiv_master
}; // class ProcessGroup

// COMMENT_HILMAR: every process group will certainly need its own MPI buffer... then activate this code.
// // initialisation of static members
// COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
//int ProcessGroup::BUFFERSIZE_BYTES = 4194304;




/**
* \brief class describing a work group, i.e. a set of worker processes sharing the same MPI communicator and performing
*        the same task
*
* A WorkGroup object managing n processes typically consists of one Worker object (living on this process) and n-1
* RemoteProcess objects. WorkGroup objects are created by the load balancer. Example:
* The process group of a load balancer consists of 6 processes. One process (either that with local rank 0 or the
* dedicated one) reads the mesh and the solver configuration and decides that the coarse grid problem is to be
* treated by 2 processes (local ranks 0 and 1) and the fine grid problems by all 6 processes. Then 2 work
* groups are created: one consisting of the 2 processes with local ranks 0 and 1, and one consisting of all 6
* processes, each of the work groups having its own MPI communicator. In building the work groups, also the
* corresponding (Remote)Worker objects are created. Communication between different work groups or with the dedicated
* load balancer process is done via the enclosing ProcessGroup communicator.
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*
*/
class WorkGroup
  : public ProcessGroup
{

private:

  /* *****************
  * member variables *
  *******************/
//  /// worker object living on this process
//  Worker* _worker;

//  /**
//  * \brief vector of remote workers in the work group
//  *
//  * Here, RemoteWorker objects are used (instead of Worker objects) since they exist on remote processes.
//  */
//  std::vector<RemoteWorker*> _remote_workers;


public:

  /* *************************
  * constructor & destructor *
  ***************************/
  /// constructor
  WorkGroup(
    const int num_processes,
    int ranks_group_parent[],
    ProcessGroup* process_group_parent,
    const int group_id)
    : ProcessGroup(num_processes, ranks_group_parent, process_group_parent, group_id)
//      _worker(nullptr),
//      _remote_workers(nullptr)
  {
//    _worker = new Worker(_comm, _rank, process_group_parent->comm(), _process_group_parent->rank());


    /* ******************************
    * test the logger functionality *
    ********************************/

    // let the coordinator of the work group (rank = 0) trigger some common messages
    if(is_coordinator())
    {
      std::string s("I have COMM_WORLD rank " + StringUtils::stringify(Process::rank)
                    + " and I am the coordinator of work group " + StringUtils::stringify(_group_id));
      Logger::log_master("Hello, master screen! " + s, Logger::SCREEN);
      Logger::log_master("Hello, master file! " + s, Logger::FILE);
      Logger::log_master("Hello, master screen and file! " + s, Logger::SCREEN_FILE);
      Logger::log_master("Hello, default master screen and file! " + s);
    }

    // write some individual messages to screen and file
    std::string s("I have COMM_WORLD rank " + StringUtils::stringify(Process::rank)
                  + " and group rank " + StringUtils::stringify(_rank)
                  + " in work group " + StringUtils::stringify(_group_id));
    // vary the lengths of the messages a little bit
    for(int i(0) ; i < _rank+1 ; ++i)
    {
      s += " RANK";
    }
    for(int i(0) ; i < _group_id+1 ; ++i)
    {
      s += " GROUP";
    }
    log_indiv_master("Hello, master screen! " + s, Logger::SCREEN);
    log_indiv_master("Hello, master file! " + s, Logger::FILE);
    log_indiv_master("Hello, master screen and file! " + s, Logger::SCREEN_FILE);
    log_indiv_master("Hello, master default screen and file! " + s);
  }

  /// destructor
  ~WorkGroup()
  {
//    delete _worker;
  }
};

#endif // guard KERNEL_PROCESS_GROUP_HPP

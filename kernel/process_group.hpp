/* GENERAL_REMARK_BY_HILMAR:
 * If buffered communication is necessary (not sure...), each process group needs its own buffer. This has to be
 * organised somehow. There are some rudimentary preparation for this (commented out).
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_PROCESS_GROUP_HPP
#define KERNEL_PROCESS_GROUP_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/assertion.hpp>
#include <kernel/logger.hpp>
#include <kernel/process.hpp>
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
    unsigned int _num_processes;

    /// rank of this process with respect to the ProcessGroup's communicator
    int _rank;

    /**
    * \brief ID of the process group
    *
    * When several process groups are "spawned" from one parent process group, then these process groups can be
    * distinguished via this group ID. It has to be passed to the process group constructor, it is not incremented
    * automatically.
    */
    unsigned int const _group_id;

    /// rank of the coordinator process within the group
    int const _rank_coord;

//   COMMENT_HILMAR: every process group will eventually need its own MPI buffer... then activate this code.
//    /// buffer size
//    static unsigned int BUFFERSIZE_BYTES;
//
//    /// buffer for MPI communication
//    char* _buffer;


  public:

    /* **************************
    * constructors & destructor *
    ****************************/
    /**
    * \brief CTORs for the case the MPI_Group and the MPI_Comm already exist
    *
    * This constructor is, e.g., intended for creating a process group object containing all COMM_WORLD processes or
    * all COMM_WORLD processes except the master. The coordinator is set to rank 0, i.e., the first rank in the process
    * group. The group ID is set to the given group ID, or 0, if no ID is given.
    *
    * \param[in] comm
    * communicator shared by the group processes
    *
    * \param[in] group_id
    * ID of the process group (default value 0 when not provided)
    */
    ProcessGroup(
      MPI_Comm comm,
      unsigned int const group_id = 0)
      : _comm(comm),
        _num_processes(0),
        _group_id(group_id),
        _rank_coord(0)
    {
      CONTEXT("ProcessGroup::ProcessGroup()");
      // get MPI_Group object for the given communicator
      int mpi_error_code = MPI_Comm_group(_comm, &_group);
      validate_error_code_mpi(mpi_error_code, "MPI_Comm_group");

      // get the size of the group
      int num_proc;
      mpi_error_code = MPI_Comm_size(_comm, &num_proc);
      validate_error_code_mpi(mpi_error_code, "MPI_Comm_size");
      ASSERT(num_proc >= 1, "Number of processes must be at least 1.");
      _num_processes = (unsigned int) num_proc;

      // and finally look up the local rank of this process w.r.t. the group's communicator
      mpi_error_code = MPI_Group_rank(_group, &_rank);
      validate_error_code_mpi(mpi_error_code, "MPI_Group_rank");

// COMMENT_HILMAR: every process group will eventually need its own MPI buffer... then activate this code.
//      _buffer = nullptr;
    }


    /**
    * \brief CTOR for the case the parent group is split in disjunct process groups
    *
    * This constructor can be used for splitting the complete set of COMM_WORLD processes into subgroups for performing
    * two or more completely separated tasks (e.g., for multiphysics).
    * The coordinator of this process group is set to rank 0, i.e. the first rank in the process group.
    *
    * \param[in] process_group_parent
    * parent group of processes
    *
    * \param[in] group_id
    * ID of this group, used as 'color' code for the MPI_Comm_split routine
    */
    ProcessGroup(
      ProcessGroup const* process_group_parent,
      unsigned int const group_id)
      : _num_processes(0),
        _group_id(group_id),
        _rank_coord(0)
    {
      CONTEXT("ProcessGroup::ProcessGroup()");

      // Split process groups according to the group id which serves as 'color' code. The 0 is the key which determines
      // which process gets which rank in the new communicator. Since all get the same key (namely 0), the ranks
      // are ordered relative to the ranks of the parent group. (I.e., when the parent group ranks are 1,42,3,6, then
      // these four processes get the local ranks 0,3,1,2.)
      int mpi_error_code = MPI_Comm_split(process_group_parent->comm(), _group_id, 0, &_comm);
      validate_error_code_mpi(mpi_error_code, "MPI_Comm_split");

      // get MPI_Group object for the given communicator
      mpi_error_code = MPI_Comm_group(_comm, &_group);
      validate_error_code_mpi(mpi_error_code, "MPI_Comm_group");

      // get the size of the group
      int num_proc;
      mpi_error_code = MPI_Comm_size(_comm, &num_proc);
      validate_error_code_mpi(mpi_error_code, "MPI_Comm_size");
      ASSERT(num_proc >= 1, "Number of processes must be at least 1.");
      _num_processes = (unsigned int) num_proc;

      // and finally look up the local rank of this process w.r.t. the group's communicator
      mpi_error_code = MPI_Group_rank(_group, &_rank);
      validate_error_code_mpi(mpi_error_code, "MPI_Group_rank");

// COMMENT_HILMAR: every process group will eventually need its own MPI buffer... then activate this code.
//      // allocate communication buffers
//      int size_of_char;
//      mpi_error_code = MPI_Pack_size(1, MPI_CHAR, MPI_COMM_WORLD, &size_of_char);
//      validate_error_code_mpi(mpi_error_code, "MPI_Pack_size");
//
//      if (size_of_char != 1)
//      {
//        // In the unlikely case that a char is not 1 byte, determine the size of the buffer. The -0.001 is a precaution
//        // for the case the division results in 42.00000000001 instead of 42.0.
//        int buffer_size = ceil(BUFFERSIZE_BYTES/size_of_char - 0.001);
//        _buffer = new char[buffer_size];
//      }
//      else
//      {
//        // otherwise, array size (in bytes) and number of array entries are equal
//        _buffer = new char[BUFFERSIZE_BYTES];
//      }
    }


    /// DTOR
    virtual ~ProcessGroup()
    {
      CONTEXT("ProcessGroup::~ProcessGroup()");
      if(_comm != MPI_COMM_WORLD && _comm != MPI_COMM_NULL)
      {
        MPI_Comm_free(&_comm);
        MPI_Group_free(&_group);
      }
// COMMENT_HILMAR: every process group will eventually need its own MPI buffer... then activate this code.
//      if (_buffer != nullptr)
//      {
//        delete [] _buffer;
//      }
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
      CONTEXT("ProcessGroup::group()");
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
      CONTEXT("ProcessGroup::comm()");
      return _comm;
    }

    /**
    * \brief getter for the number of processes
    *
    * \return number of processes #_num_processes
    */
    inline unsigned int num_processes() const
    {
      CONTEXT("ProcessGroup::num_processes()");
      return _num_processes;
    }

    /**
    * \brief getter for group ID
    *
    * \return group ID #_group_id
    */
    inline unsigned int group_id() const
    {
      CONTEXT("ProcessGroup::group_id()");
      return _group_id;
    }

    /**
    * \brief getter for the rank w.r.t. the group communicator this process belongs to
    *
    * \return rank #_rank w.r.t. process group communicator
    */
    inline int rank() const
    {
      CONTEXT("ProcessGroup::rank()");
      return _rank;
    }

    /**
    * \brief getter for the rank w.r.t. the group communicator this process belongs to
    *
    * \return rank #_rank w.r.t. process group communicator
    */
    inline int rank_coord() const
    {
      CONTEXT("ProcessGroup::rank_coord()");
      return _rank_coord;
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
      CONTEXT("ProcessGroup::is_coordinator()");
      return _rank == _rank_coord;
    }


    /**
    * \brief sending individual messages to the master
    *
    * With this function the processes of the process group can send individual log messages to the master. The
    * coordinator of the process group collects all these log messages in the order of process group ranks and sends
    * them as one MPI message to the master. Depending on the target the master then displays these messages on the
    * screen or writes them to the log file (one line per message). This function has to be called by \em all processes
    * of the process group! Note that no further information will be appended to the messages (like rank of the process
    * from which the messages comes). Such infomration has to be coded directly into the message string.
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
      CONTEXT("ProcessGroup::log_indiv_master()");
      // add 1 to the message length since string::c_str() adds the null termination symbol to the resulting char array
      unsigned int length = unsigned int(message.length()) + 1;

      if (!is_coordinator())
      {
        /* ***************************************
        * code for all non-coordinator processes *
        *****************************************/

        // coordinator process gathers the lengths of the messages
        MPI_Gather(&length, 1, MPI_UNSIGNED, nullptr, 0, MPI_DATATYPE_NULL, _rank_coord, _comm);

        // coordinator process gathers the messages
        // (Here, it is necessary to cast away the const'ness of string::c_str() since the MPI routine expects a
        // non-const send buffer.)
        MPI_Gatherv(const_cast<char *>(message.c_str()), length, MPI_CHAR, nullptr, nullptr, nullptr,
                    MPI_DATATYPE_NULL, _rank_coord, _comm);
      }
      else
      {
        /* ********************************
        * code for the coordinator process *
        **********************************/

        // receive buffer for storing the lengths of the messages
        unsigned int* msg_lengths = new unsigned int[_num_processes];
        // array for storing the start positions of the single messages in the receive buffer
        unsigned int* msg_start_pos = new unsigned int[_num_processes];

        // gather the lengths of the messages from the other processes
        MPI_Gather(&length, 1, MPI_UNSIGNED, msg_lengths, 1, MPI_UNSIGNED, _rank_coord, _comm);

        // set start positions of the single messages in the receive buffer
        msg_start_pos[0] = 0;
        for(unsigned int i(1) ; i < _num_processes ; ++i)
        {
          msg_start_pos[i] = msg_start_pos[i-1] + msg_lengths[i-1];
        }

        // determine total length of the messages and allocate receive buffer accordingly
        unsigned int total_length(msg_start_pos[_num_processes-1] + msg_lengths[_num_processes-1]);

        // receive buffer for (consecutively) storing the messages
        char* messages = new char[total_length];

//        // debug output
//        for(unsigned int i(0) ; i < _num_processes ; ++i)
//        {
//          std::cout << "local process " << stringify(i) << ": length of string = "
//                    << stringify(msg_lengths[i]) << std::endl;
//        }
//        std::cout << "Total length: " << stringify(total_length) << std::endl;

        // gather the messages from the other processes.
        MPI_Gatherv(const_cast<char *>(message.c_str()), length, MPI_CHAR, messages,
                    reinterpret_cast<int*>(msg_lengths), reinterpret_cast<int*>(msg_start_pos), MPI_CHAR, _rank_coord,
                    _comm);

//        // debug output
//        // output the strings collected from all processes by using corresponding offsets in the
//        // char array (pointer arithmetic)
//        for(unsigned int i(0) ; i < _num_processes ; ++i)
//        {
//          std::cout << "Process " << Process::rank << " writes: "
//                    << std::string(messages + msg_start_pos[i], msg_lengths[i]) << std::endl;
//        }

        // now send everything to the master
        Logger::log_master_array(_num_processes, msg_lengths, total_length, messages, target);
        delete [] msg_lengths;
        delete [] msg_start_pos;
        delete [] messages;
      }
    } // log_indiv_master
  }; // class ProcessGroup

// COMMENT_HILMAR: every process group will certainly need its own MPI buffer... then activate this code.
// // initialisation of static members
// COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
//int ProcessGroup::BUFFERSIZE_BYTES = 4194304;
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_BM_PROCESS_GROUP_HPP

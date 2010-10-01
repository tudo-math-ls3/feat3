#pragma once
#ifndef KERNEL_LOGGER_HPP
#define KERNEL_LOGGER_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/service_ids.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/communication.hpp>
#include <kernel/process_group.hpp>


/**
* \brief class providing logging mechanisms
*
* <ul>
*   <li>
*   Each process is connected to a single log file. The log directory \c logdir, as well as the base name \c basename
*   of the log files is set in some basic configuration file. The name of the log file is '\c logdir/basename\<n\>\c .log'
*   where \c n is the MPI_COMM_WORLD rank of the process. \c n is displayed with at least three digits using leading
*   zeros. When 1000 or more MPI processes are used, the number of digits is automatically increased. The basename is
*   empty by default, and the log directory is \c ./log by default.
*   <li>
*   Only the master process is allowed to produce screen log.
*   <li>
*   The master file log and the screen log are basically independent of each other, i.e. there are three functions
*     \code log_master(..., LOGGER::SCREEN)       (--> message only appears on the screen)\endcode
*     \code log_master(..., LOGGER::FILE)         (--> message only appears in the master log file)\endcode
*     \code log_master(..., LOGGER::SCREEN_FILE)  (--> message appears on the screen and in the master log file)\endcode
*   which can be called by every process. When a non-master process calls them, communication is involved (sending the
*   message to the master).
*   <li>
*   When a group of processes wants to trigger individual messages (see logging scenario 6), there are the three
*   variants
*     \code log_indiv_master(..., LOGGER::SCREEN) \endcode
*     \code log_indiv_master(..., LOGGER::FILE) \endcode
*     \code log_indiv_master(..., LOGGER::SCREEN_FILE) \endcode
*   Since this can only be done within process groups, these functions are part of the class ProcessGroup.
*   <li>
*   The user has the option to \em globally synchronise master file log and screen log in four variants:
*   <ol>
*     <li> All messages sent to the master file log automatically appear on the screen.
*     <li> All messages sent to the screen automatically appear in the master file log.
*     <li> 1. + 2.
*     <li> no synchronisation
*   </ol>
*   The default is 4.
*   COMMENT_HILMAR: This feature is not implemented yet.
*   <li>
*   Each process can write messages to its log file via the function
*     \c log(...)
*   (i.e., on the master process the functions log(...) and log_master(..., LOGGER::FILE) are equivalent)
*   <li>
*   What sort of logging scenarios are there?
*   <ol>
*     <li>
*     <em> single process writes a message to its log file:</em>\n
*       call (on the writing process):
*       \code log("point (x,y) found in element 17"); \endcode
*       result in 002.log:
*       \code point (x,y) found in element 17 \endcode
*     <li>
*     <em>single process triggers a message on screen and/or in master log file:</em>\n
*       call (on the triggering process):
*       \code log_master("point (x,y) found in element 17"); \endcode
*       result in 010.log (assuming that the master has rank 10) and on screen:
*       \code point (x,y) found in element 17 \endcode
*     <li>
*     <em>group of processes writes a common message to the log files:</em>\n
*       call (on all writing processes):
*       \code log("global solver: starting iter 23"); \endcode
*       output in 000.log, 001.log, ..., 009.log:
*       \code  global solver: starting iter 23 \endcode
*       \code  global solver: starting iter 23 \endcode
*       \code    ... \endcode
*       \code  global solver: starting iter 23 \endcode
*     <li>
*     <em>group of processes triggers a common message on screen and/or in master log file:</em>\n
*       call (on all writing processes):
*         \code if(i_am_coordinator) \endcode
*         \code { \endcode
*         \code   log_master("global solver: starting iter 23"); \endcode
*         \code } \endcode
*       output in 010.log:
*       \code global solver: starting iter 23 \endcode
*       output on screen:
*       \code global solver: starting iter 23 \endcode
*       In this case, not every process sends its message but only one 'coordinator' process. The logic whether
*       this is one common message from a group of processes has to be provided by the user (see if clause). There is
*       no special routine for this case.
*     <li>
*     <em>group of processes writes individual messages to their log files:</em>\n
*       call on process 0, process 1, ..., process 9:
*       \code log("local solver: conv. rate: 0.042"); \endcode
*       \code log("local solver: conv. rate: 0.023"); \endcode
*       \code    ...\endcode
*       \code log("local solver: conv. rate: 0.666"); \endcode
*       output in 000.log, 001.log, ..., 009.log:
*       \code local solver: conv. rate: 0.042 \endcode
*       \code local solver: conv. rate: 0.023 \endcode
*       \code    ...\endcode
*       \code local solver: conv. rate: 0.666 \endcode
*     <li>
*     <em>group of processes triggers individual messages on screen and/or in master log file:</em>\n
*       call on process 0, process 1, ..., process 9:
*       \code ProcessGroup::log_indiv_master("process 0: local solver: conv. rate: 0.042"); \endcode
*       \code ProcessGroup::log_indiv_master("process 1: local solver: conv. rate: 0.023"); \endcode
*       \code    ...\endcode
*       \code ProcessGroup::log_indiv_master("process 9: local solver: conv. rate: 0.666");\endcode
*       output in 010.log:
*       \code process 0: local solver: conv. rate: 0.042 \endcode
*       \code process 1: local solver: conv. rate: 0.023 \endcode
*       \code    ...\endcode
*       \code process 9: local solver: conv. rate: 0.666 \endcode
*       output on screen:
*       \code process 0: local solver: conv. rate: 0.042 \endcode
*       \code process 1: local solver: conv. rate: 0.023 \endcode
*       \code    ...\endcode
*       \code process 9: local solver: conv. rate: 0.666 \endcode
*      Since this can only be done within process groups, these functions are part of the class ProcessGroup.
*      Two possibilities how to transfer the messages to the master:
*      <ol>
*        <li> All processes of the group send their message to the coordinator process of the group (via MPI_Gather),
*             which collects them in some array structure and sends them to the master.
*        <li> All processes of the group send their message directly to the master, which collects them in some
*             array structure and displays/writes them.
*      </ol>
*      COMMENT_HILMAR: I'm not sure how to implement this version... When two process groups
*      send indiv. messages then the master should group them accordingly. The master can only
*      communicate via the COMM_WORLD communicator to the two process groups, so one has to distinguish them
*      via MPI tags. Imagine the first process group contains 3 processes using tag 0 and the second process
*      group 5 processes using tag 1. When the request of the first group arrives first, then the master creates
*      an array of length 3 to store three strings and starts a loop calling two further receives of strings
*      (the first one it alreday got) listening to tag 0 only. But what happens now with the requests of
*      the second process group, which can (and must not) accepted in that loop (due to the wrong tag)? Are
*      they automatically postponed? Or does the process deadlock?
*      Due to the uncertainties with the second approach, I implemented the first one.
*  </ol>
* </ul>
* \todo everything concerning file output still has to be implemented
*
* \author Hilmar Wobker
*/
class Logger
{

private:

public:

  /**
  * \brief log targets used in various logging routines
  *
  * Log targets are either the screen (target SCREEN), a log file (target FILE), or both (target SCREEN_FILE).
  * Only the master process is allowed to write to the screen.
  */
  enum target
  {
    /// only write to screen
    SCREEN,
    /// only write to the log file
    FILE,
    /// write to screen and log file
    SCREEN_FILE
  };


  /**
  * \brief triggers logging of a message (given as string) on the master process
  *
  * This function receives a string representing a log message. It triggers the function receive() in the master's
  * service loop and sends the message to the master, which writes them to screen and/or log file.
  *
  * \param[in] message
  * string representing the log message
  *
  * \param[in] targ
  * output target SCREEN, FILE or SCREEN_FILE (default: SCREEN_FILE)
  *
  * \sa receive()
  *
  * \author Hilmar Wobker
  */
  static void log_master(
    std::string const message,
    target targ = SCREEN_FILE)
  {

// TODO: write some simple wrapper routine for the following MPI calls (like in FEAST1): init_msg, write_msg, send_msg
    // reset buffer
    Comm::MCW_buffer_pos = 0;
    // write the message id to the buffer
    int id = ServiceIDs::LOG_RECEIVE;
    int mpi_error_code = MPI_Pack(&id, 1, MPI_INTEGER, Comm::MCW_buffer,
                                  Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

// TODO: define specific MPI_Datatype for the following data (see example on p. 125f in MPI2.2 standard)

    // add 1 to the message length
    int msg_length(message.size()+1);

    // write length of the message to the buffer
    mpi_error_code = MPI_Pack(&msg_length, 1, MPI_INTEGER, Comm::MCW_buffer,
                              Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

    // write char array itself to the buffer
    mpi_error_code = MPI_Pack(const_cast<char *>(message.c_str()), msg_length, MPI_CHARACTER, Comm::MCW_buffer,
                              Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
// COMMENT_HILMAR: Legitimate to simply cast the const away? (see comment in process_group.hpp)

    // write log target to the buffer
    mpi_error_code = MPI_Pack(&targ, 1, MPI_INTEGER, Comm::MCW_buffer,
                              Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

    // send message
    mpi_error_code = MPI_Send(Comm::MCW_buffer, Comm::MCW_buffer_pos, MPI_PACKED, Process::rank_master, 0,
                              MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Send");
  }


  /**
  * \brief writing a messages to the screen and/or to the master's log file
  *
  * This function runs on the master and is triggered by the function log_master(). It receives one MPI message
  * consisting of a char array representing one log message. Depending on the sent output target the function then
  * writes the message to the screen and/or to the log file.
  *
  * \note Of course, one could realise this via the function receive_array(). But since single messages are sent much
  * more often then arrays of messages, it is justified to use a specially tailored function for this.
  *
  * \note There is no way to know in which order the master receives messages concurrently sent from different
  * processes. But when the master reacts to one request, then it will complete it before doing something else.
  * Hence, these concurrent requests should not be problematic.
  *
  * \sa log_master
  *
  * \author Hilmar Wobker
  */
  static void receive()
  {
    // read length of the messages from the buffer
    int msg_length;
    int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, &msg_length,
                                1, MPI_INTEGER, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");

//    // debug output
//    std::cout << "read length of message: " << msg_length << std::endl;

    // char array for storing the message
    char message[msg_length];
    // read char array from the buffer
    mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, message,
                                msg_length, MPI_CHAR, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");

    // read log target from the buffer
    int target(0);
    mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, &target, 1,
                                MPI_INTEGER, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");

    // display message on screen if requested
    if (target == SCREEN || target == SCREEN_FILE)
    {
      // use corresponding offsets in the char array (pointer arithmetic)
      std::cout << message << std::endl;
    }

    // write messages to master's log file if requested
    if (target == FILE || target == SCREEN_FILE)
    {
      // use corresponding offsets in the char array (pointer arithmetic)
      //bla_write_to_file_blub(std::string(message, msg_length));

      // as long as not file output is implemented, write messages to the screen
      std::cout << "FILE OUTPUT: " << message << std::endl;
    }
  }


  /**
  * \brief triggers logging of distinct messages (given as one char array) on the master process
  *
  * This function receives one large char array which contains a number of distinct messages. It triggers the function
  * receive_array() in the master's service loop and sends the messages to the master, which writes them to
  * screen and/or log file.
  *
  * \param[in] num_messages
  * the number of messages the char array \a messages contains
  *
  * \param[in] msg_lengths
  * array of lengths of the single messages, dimension: [\a num_messages]
  *
  * \param[in] total_length
  * total length of all messages (could also be computed in this function, but is often already available outside)
  *
  * \param[in] messages
  * char array containing the messages (each message terminated by null symbol), dimension: [\a total_length]
  *
  * \param[in] targ
  * output target SCREEN, FILE or SCREEN_FILE (default: SCREEN_FILE)
  *
  * \sa log_master_array(std::vector<std::string> const,target), receive_array
  *
  * \author Hilmar Wobker
  */
  static void log_master_array(
    int num_messages,
    int* const msg_lengths,
    int const total_length,
    char* const messages,
    target targ = SCREEN_FILE)
  {
// TODO: write some simple wrapper routine for the following MPI calls (like in FEAST1): init_msg, write_msg, send_msg
    // reset buffer
    Comm::MCW_buffer_pos = 0;
    // write the message id to the buffer
    int id = ServiceIDs::LOG_RECEIVE_ARRAY;
    int mpi_error_code = MPI_Pack(&id, 1, MPI_INTEGER, Comm::MCW_buffer,
                                  Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

// TODO: define specific MPI_Datatype for the following data (see example on p. 125f in MPI2.2 standard)

    // write number of messages the char array consists of to the buffer
    mpi_error_code = MPI_Pack(&num_messages, 1, MPI_INTEGER, Comm::MCW_buffer,
                              Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

    // write lengths of the messages to the buffer
    mpi_error_code = MPI_Pack(msg_lengths, num_messages, MPI_INTEGER, Comm::MCW_buffer,
                              Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

    // write char array itself to the buffer
    mpi_error_code = MPI_Pack(messages, total_length, MPI_CHARACTER, Comm::MCW_buffer,
                              Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

    // write log target to the buffer
    mpi_error_code = MPI_Pack(&targ, 1, MPI_INTEGER, Comm::MCW_buffer,
                              Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");

    // send message
    mpi_error_code = MPI_Send(Comm::MCW_buffer, Comm::MCW_buffer_pos, MPI_PACKED, Process::rank_master, 0,
                              MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Send");
  }


  /**
  * \brief triggers logging of distinct messages (given as vector of strings) on the master process
  *
  * This function receives a vector of strings representing distinct messages. It converts the messages to one large
  * char array and calls the first version of the routine log_master_array(), which then triggers logging of these
  * messages on the master process.
  *
  * \param[in] messages
  * vector of strings representing the distinct messages
  *
  * \param[in] targ
  * output target SCREEN, FILE or SCREEN_FILE (default: SCREEN_FILE)
  *
  * \sa log_master_array(int, int*, int, char*, target), receive_array
  *
  * \author Hilmar Wobker
  */
  static void log_master_array(
    std::vector<std::string> const &messages,
    target targ = SCREEN_FILE)
  {

// COMMENT_HILMAR: Maybe it is more efficient to *not* realise this via the first version of log_master_array(), but
// to trigger an extra receive routine on master side...

    // convert the vector of strings into one long char array and determine further information needed by
    // the other version of the function log_master_array(...).
    int num_messages(messages.size());
    int msg_lengths[num_messages];
    int total_length(0);
    for(int i(0) ; i < num_messages ; ++i)
    {
      // add 1 due to the null termination symbol added by string::c_str()
      msg_lengths[i] = messages[i].size()+1;
      total_length += msg_lengths[i];
    }
    char messages_char[total_length];
    // set pointer to beginning of message
    char* pos = messages_char;
    for(int i(0) ; i < num_messages ; ++i)
    {
      // copy message to corresponding part of the char array
      strcpy(pos, messages[i].c_str());
      // move pointer
      pos += msg_lengths[i];
    }
    // now call the other version of the function log_master_array()
    log_master_array(num_messages, msg_lengths, total_length, messages_char, targ);
  } // log_master_array()


  /**
  * \brief writing a number of distinct messages to the screen and/or to the master's log file
  *
  * This function runs on the master and is triggered by the (overloaded) function log_master_array(). It receives one
  * MPI message consisting of one long char array representing an array of distinct log messages, plus further
  * information how to partition this char array into the single messages. Depending on the sent output target the
  * function then writes one line per message to the screen and/or to the log file.
  *
  * \note There is no way to know in which order the master receives messages concurrently sent from different
  * processes. But when the master reacts to one request, then it will complete it before doing something else.
  * Hence, these concurrent requests should not be problematic.
  *
  * \sa log_master_array(int, int*, int, char*, target), log_master_array(std::vector<std::string> const,target)
  *
  * \author Hilmar Wobker
  */
  static void receive_array()
  {

// COMMENT_HILMAR:
// This function allocates and deallocates three arrays. Maybe it is more efficient to use and reuse some
// "pre-allocated" storage...

    // read number of messages the char array consists of from to the buffer
    int num_messages(0);

    int mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, &num_messages,
                                    1, MPI_INTEGER, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");

//    // debug output
//    std::cout << "read num_messages: " << num_messages << std::endl;

    // allocate array for sotring message lengths
    int msg_lengths[num_messages];
    // read message lengths from the buffer
    mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, msg_lengths,
                                num_messages, MPI_INTEGER, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");

//    // debug output
//    for (int i(0) ; i < num_messages ; ++i)
//    {
//      std::cout << "read length of message " << i << ": " << msg_lengths[i] << std::endl;
//    }

    // allocate array for storing the start positions of the single messages in the receive buffer
    int msg_start_pos[num_messages];

    // set start positions of the single messages in the receive buffer
    msg_start_pos[0] = 0;
    for(int i(1) ; i < num_messages ; ++i)
    {
      msg_start_pos[i] = msg_start_pos[i-1] + msg_lengths[i-1];
    }
    int total_length(msg_start_pos[num_messages-1] + msg_lengths[num_messages-1]);

//    // debug output
//    std::cout << "total length of message " << total_length << std::endl;

    // allocate char array for (consecutively) storing the messages
    char messages[total_length];
    // read char array itself from the buffer
    mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, messages,
                                total_length, MPI_CHAR, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");

    // read log target from the buffer
    int target(0);
    mpi_error_code = MPI_Unpack(Comm::MCW_buffer, Comm::MCW_received_bytes, &Comm::MCW_buffer_pos, &target,
                                    1, MPI_INTEGER, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Unpack");

    // display messages on screen if requested
    if (target == SCREEN || target == SCREEN_FILE)
    {
      for(int i(0) ; i < num_messages ; ++i)
      {
        // use corresponding offsets in the char array (pointer arithmetic)
        std::cout << messages + msg_start_pos[i] << std::endl;
      }
    }

    // write messages to master's log file if requested
    if (target == FILE || target == SCREEN_FILE)
    {
      for(int i(0) ; i < num_messages ; ++i)
      {
        // use corresponding offsets in the char array (pointer arithmetic)
        //bla_write_to_file_blub(messages + msg_start_pos[i]));

        // as long as not file output is implemented, write messages to the screen
        std::cout << "FILE OUTPUT: " << messages + msg_start_pos[i] << std::endl;
      }
    }
  }
}; // class Logger

#endif // guard KERNEL_LOGGER_HPP

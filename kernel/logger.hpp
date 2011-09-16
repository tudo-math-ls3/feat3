/* GENERAL_REMARK_BY_HILMAR:
 * See issue_00023.txt and COMMENT_HILMAR in this file.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_LOGGER_HPP
#define KERNEL_LOGGER_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/service_ids.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/comm.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

namespace FEAST
{
  /**
  * \brief class providing logging mechanisms
  *
  * <ul>
  *   <li>
  *   Each process is connected to a single log file. The log directory \c logdir, as well as the base name \c basename
  *   of the log files is set in some basic configuration file. The name of the log file is
  *   '\c logdir/basename\<n\>\c .log' where \c n is the MPI_COMM_WORLD rank of the process. \c n is displayed with at
  *   least three  digits using leading zeros. When 1000 or more MPI processes are used, the number of digits is
  *   automatically increased. The basename is empty by default, and the log directory is \c ./log by default.
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
  * COMMENT_HILMAR: When two process groups send indiv. messages then the master has to group them
  * accordingly. The master can only communicate via the COMM_WORLD communicator to the two process groups,
  * so one has to distinguish them via MPI tags. Imagine the first process group contains 3 processes using
  * tag 0 and the second process group 5 processes using tag 1. When the request of the first group arrives
  * first, then the master creates an array of length 3 to store three strings and starts a loop calling
  * two further receives of strings (the first one it alreday got) listening to tag 0 only. The requests of
  * the second process group are postponed by the master automatically.
  * Currently, the first version is implemented. It might be more efficient, however, to use the second version
  * (which doesn't need synchronisation between the members of the process group).
  *      </ol>
  *  </ol>
  * </ul>
  *
  * \author Hilmar Wobker
  */
  class Logger
  {

  public:

    /// variable storing the base name of the log file
    static std::string file_base_name;

    /// variable storing the complete name of the log file
    static std::string file_name;

    /// log file stream
    static std::ofstream file;

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
    * \brief opens the log file
    *
    * This function opens a log file. It has to be called on each process. The filename is assembled as
    * '<base_name><world_rank>.log', i.e. each MPI_COMM_WORLD process has his own log file.
    *
    * \param[in] base_name
    * The base name of the log file
    *
    * \todo Currently, the file is simply created in the current directory. E.g., when starting the automatic test
    * system ('ctest test'), all files are created in this root directory. When starting the program manually, e.g.,
    * in the kernel folder, the files are created there. So, the system has to be enhanced so that the user can specify
    * the absolute path to the desired log file folder (Windows and Unix).
    *
    * \author Hilmar Wobker
    */
    static void open_log_file(std::string const base_name)
    {
      CONTEXT("Logger::open_log_file()");
      // set base name of the log file
      file_base_name = base_name;
      // add the MPI_COMM_WORLD rank of this process to the base name
      // inquire how many digits are needed (depending on the number of processes)
      int num_digits((int)log10((double)Process::num_processes) + 1);
      // use at least 3 digits
      num_digits = std::max(3, num_digits);
      // +1 needed for termination symbol
      char* s = new char[num_digits+1];
      // use format with left-padded zeros
      sprintf(s, "%0*i", num_digits, Process::rank);
      // set full filename
      file_name = file_base_name + s + ".log";
      delete [] s;

      if (!file.is_open())
      {
        file.open(file_name.c_str());
        // catch file opening error
        if (file.fail())
        {
          abort_mpi("Error! Could not open log file " + file_name + ".");
        }
      }
      else
      {
        abort_mpi("Error! Log file " + file_name + " is already opened!");
      }
    }


    /**
    * \brief closes the log file
    *
    * \author Hilmar Wobker
    */
    static void close_log_file()
    {
      CONTEXT("Logger::close_log_file()");
      if (file.is_open())
      {
        file.close();
      }
      else
      {
        abort_mpi("Error! Log file " + file_name + " cannot be closed since it does not seems to have been "
                  + "openend before...!");
      }
    }


    /**
    * \brief writes a message to the log file of this process
    *
    * This function receives a string representing a log message and writes it to the log file attached to this process.
    * Note that no line break is performed.
    *
    * \param[in] message
    * string representing the log message
    *
    * \author Hilmar Wobker
    */
    static void log(std::string const message)
    {
      CONTEXT("Logger::log()");
      file << message;
    }


// COMMENT_HILMAR: This version seems not to be necessary. The above function log(std::string const &message) is also
// able to process char pointers (see for example the call in Logger::receive_array())
// Using this version here additionally triggers the gcc compiler to throw a warning when passing a constant string
// to the log(...) function, e.g. 'log("BRAL")'.
//   ...: warning: deprecated conversion from string constant to 'char*'
//  static void log(char message[])
//  {
//    file << message;
//  }


    /**
    * \brief writes a number of messages to the log file of this process
    *
    * This function receives a vector of strings representing a number of log message and writes them to the log file
    * attached to this process. Note that no line break is performed.
    *
    * \param[in] messages
    * vector of strings representing the messages
    *
    * \author Hilmar Wobker
    */
    static void log(std::vector<std::string> const &messages)
    {
      CONTEXT("Logger::log()");
      for(unsigned int i(0) ; i<messages.size() ; ++i)
      {
        file << messages[i];
      }
    }


    /**
    * \brief triggers logging of a message (given as string) on the master process
    *
    * This function receives a string representing a log message. It triggers the function receive() in the master's
    * service loop and sends the message to the master, which writes them to screen and/or log file. Note that no line
    * break is performed.
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
      std::string const &message,
      target targ = SCREEN_FILE)
    {
      CONTEXT("Logger::log_master()");
      // init a new message with corresponding ID
      Comm::init(ServiceIDs::LOG_RECEIVE);

// TODO: define specific MPI_Datatype for the following data? (see example on p. 125f in MPI2.2 standard)

      // write length of the log message to the buffer (add 1 to the length since string::c_str() adds null termination
      // symbol)
      Comm::write((unsigned int)message.size()+1);

      // write string itself to the buffer
      Comm::write(message);

      // write log target to the buffer
      Comm::write((int)targ);

      // send message
      Comm::send();
    }


    /**
    * \brief receiving a log message and writing it to the screen and/or to the master's log file
    *
    * This function runs on the master and is triggered by the function log_master(). It receives one MPI message
    * consisting of a char array representing one log message. Depending on the sent output target the function then
    * writes the message to the screen and/or to the log file. Note that no line break is performed.
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
      CONTEXT("Logger::receive()");
      // read length of the messages from the buffer
      unsigned int msg_length;
      Comm::read(msg_length);

      // char array for storing the message
      char* message = new char[msg_length];
      // read char array from the buffer
      Comm::read(msg_length, message);

      // read log target from the buffer
      unsigned int target;
      Comm::read(target);

      // display message on screen if requested
      if (target == SCREEN || target == SCREEN_FILE)
      {
        std::cout << message;
      }

      // write messages to master's log file if requested
      if (target == FILE || target == SCREEN_FILE)
      {
        log(message);
      }
      delete [] message;
    }


    /**
    * \brief triggers logging of distinct messages (given as one char array) on the master process
    *
    * This function receives one large char array which contains a number of distinct messages. It triggers the function
    * receive_array() in the master's service loop and sends the messages to the master, which writes them to
    * screen and/or log file. Note that no line break is performed.
    *
    * \param[in] num_messages
    * the number of messages the char array \a messages contains
    *
    * \param[in] msg_lengths
    * array of lengths of the single messages (including the null termination symbol), dimension: [\a num_messages]
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
      unsigned int const num_messages,
      unsigned int const* msg_lengths,
      unsigned int const total_length,
      char const* messages,
      target targ = SCREEN_FILE)
    {
      CONTEXT("Logger::log_master_array()");
      // init a new message with corresponding ID
      Comm::init(ServiceIDs::LOG_RECEIVE_ARRAY);

// TODO: define specific MPI_Datatype for the following data? (see example on p. 125f in MPI2.2 standard)

      // write number of log messages the char array consists of to the buffer
      Comm::write(num_messages);

      // write lengths of the log messages to the buffer
      Comm::write(msg_lengths, num_messages);

      // write char array itself to the buffer
      Comm::write(messages, total_length);

      // write log target to the buffer
      Comm::write((int)targ);

      // send message
      Comm::send();
    }


    /**
    * \brief triggers logging of distinct messages (given as vector of strings) on the master process
    *
    * This function receives a vector of strings representing distinct messages. It converts the messages to one long
    * string and calls the first version of the routine log_master_array(), which then triggers logging of these
    * messages on the master process. Note that no line break is performed.
    *
    * \param[in] messages
    * vector of strings representing the distinct messages
    *
    * \param[in] targ
    * output target SCREEN, FILE or SCREEN_FILE (default: SCREEN_FILE)
    *
    * \sa log_master_array(unsigned int, unsigned int*, unsigned int, char*, target), receive_array
    *
    * \author Hilmar Wobker
    */
    static void log_master_array(
      std::vector<std::string> const &messages,
      target targ = SCREEN_FILE)
    {
      CONTEXT("Logger::log_master_array()");

// COMMENT_HILMAR: Maybe it is more efficient to *not* realise this via the first version of log_master_array(), but
// to trigger an extra service receive routine on master side...

      // get number of messages
      unsigned int num_messages(messages.size());
      ASSERT(num_messages > 0, "There must be at least one message in the message vector!");

      // store lengths of the single strings
      unsigned int* msg_lengths = new unsigned int[num_messages];
      for(unsigned int i(0) ; i < num_messages ; ++i)
      {
        // add 1 due to the null termination symbol which we insert manually
        msg_lengths[i] = messages[i].size() + 1;
      }

      // Convert the vector of strings into one long string, where the single strings are separated by null termination
      // symbols '\0' (which eases output of the resulting char array on receiver side). Do not add '\0' to the last
      // string since this is added automatically by the c_str() method used below.

      // initialise the string with the first message
      std::string msgs_as_one_string(messages[0]);
      for(unsigned int i(1) ; i < num_messages ; ++i)
      {
        // append null termination symbol
        msgs_as_one_string.append(1,'\0');
        // append next string
        msgs_as_one_string.append(messages[i]);
      }

      // now call the other version of the function log_master_array(), pass char array representation of the string
      log_master_array(num_messages, msg_lengths, msgs_as_one_string.size()+1, msgs_as_one_string.c_str(), targ);

      // delete aux. array again
      delete [] msg_lengths;
    } // log_master_array()


    /**
    * \brief receiving a number of distinct messages and writing them to the screen and/or to the master's log file
    *
    * This function runs on the master and is triggered by the (overloaded) function log_master_array(). It receives one
    * MPI message consisting of one long char array representing an array of distinct log messages, plus further
    * information how to partition this char array into the single messages. Depending on the sent output target the
    * function then writes one line per message to the screen and/or to the log file. Note that no line break is
    * performed.
    *
    * \note There is no way to know in which order the master receives messages concurrently sent from different
    * processes. But when the master reacts to one request, then it will complete it before doing something else.
    * Hence, these concurrent requests should not be problematic.
    *
    * \sa log_master_array(unsigned int, unsigned int*, unsigned int, char*, target),
    *     log_master_array(std::vector<std::string> const,target)
    *
    * \author Hilmar Wobker
    */
    static void receive_array()
    {
      CONTEXT("Logger::receive_array()");

// COMMENT_HILMAR:
// This function allocates and deallocates three arrays. Maybe it is more efficient to use and reuse some
// "pre-allocated" storage...

      // read number of messages the char array consists of from to the buffer
      unsigned int num_messages;
      Comm::read(num_messages);

      // allocate array for sotring message lengths
      unsigned int* msg_lengths = new unsigned int[num_messages];
      // read message lengths from the buffer
      Comm::read(num_messages, msg_lengths);

      // allocate array for storing the start positions of the single messages in the receive buffer
      unsigned int* msg_start_pos = new unsigned int[num_messages];

      // set start positions of the single messages in the receive buffer
      msg_start_pos[0] = 0;
      for(unsigned int i(1) ; i < num_messages ; ++i)
      {
        msg_start_pos[i] = msg_start_pos[i-1] + msg_lengths[i-1];
      }
      unsigned int total_length(msg_start_pos[num_messages-1] + msg_lengths[num_messages-1]);
      delete [] msg_lengths;

      // allocate char array for (consecutively) storing the messages
      char* messages = new char[total_length];
      // read char array itself from the buffer
      Comm::read(total_length, messages);

      // read log target from the buffer
      unsigned int target;
      Comm::read(target);

      // display messages on screen if requested
      if (target == SCREEN || target == SCREEN_FILE)
      {
        for(unsigned int i(0) ; i < num_messages ; ++i)
        {
          // use corresponding offsets in the char array (pointer arithmetic)
          std::cout << messages + msg_start_pos[i];
        }
      }

      // write messages to master's log file if requested
      if (target == FILE || target == SCREEN_FILE)
      {
        for(unsigned int i(0) ; i < num_messages ; ++i)
        {
          // use corresponding offsets in the char array (pointer arithmetic)
          log(messages + msg_start_pos[i]);
        }
      }
      delete [] msg_start_pos;
      delete [] messages;
    } // receive_array()
  }; // class Logger

  std::string Logger::file_name;
  std::string Logger::file_base_name;
  std::ofstream Logger::file;

} // namespace FEAST

#endif // KERNEL_LOGGER_HPP

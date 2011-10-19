/* GENERAL_REMARK_BY_HILMAR:
 * IMHO a master process is not really necessary. However, we decided that we will stick to it for the time being.
 * The advantage of using a master process is that one can (at least up to a certain level)  produce coherent screen
 * output. E.g., two process groups are triggering screen output at the same time. Without a master process, one
 * cannot control how the two outputs appear on the screen, i.e., they can completely intermix. However, when really
 * large parallel jobs are performed, screen output doesn't make sense anyway, so this is quite a weak argument. Getting
 * individual process group output can be achieved much better by using keywords and grep commands after the
 * computation. And if screen output is used to actively follow a computation while it is running, then one can also
 * live with a cluttered output.
 * If it comes to multiphysics computations or similar, i.e. where two completely independent main process groups (and
 * load balancers) are created, and the user wants to coordinate these two computations, then this is also possible
 * without a master process.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_FOUNDATION_MASTER_HPP
#define KERNEL_FOUNDATION_MASTER_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/mpi_utils.hpp>
#include <kernel/foundation/comm.hpp>
#include <kernel/foundation/service_ids.hpp>
#include <kernel/foundation/process_group.hpp>
#include <kernel/foundation/logger.hpp>

// includes, system
#include <iostream>
#include <stdlib.h>
#include <mpi.h>

namespace FEAST
{
  /**
  * \brief class defining the master
  *
  * \author Hilmar Wobker
  * \author Dominik Goeddeke
  */
  class Master
  {

  private:

    /* *****************
    * member variables *
    *******************/
    /// flag for finishing the infinite service loop
    bool _finish_service;

#ifndef NDEBUG
    /// MPI status object needed in MPI_Recv call
    MPI_Status status;
#endif

  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /// constructor
    Master()
     : _finish_service(false)
    {
    }

    /* *****************
    * member functions *
    *******************/
    /// infinite service loop waiting for messages
    void service()
    {
      CONTEXT("Master::service()");
      std::cout << "Master at your service... " << std::endl;

      // service ID read from the sent messages
      int id;
      while(true)
      {
        // receive messages with any tag and from any source
        // (Warning! Don't replace the status object by MPI_STATUS_IGNORE! It is needed in the following call
        // of MPI_get_count().)
        MPI_Recv(Comm::MCW_buffer, Comm::MCW_BUFFERSIZE, MPI_PACKED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_MACRO);

#ifndef NDEBUG
        // Read the size of the received message (in bytes) from the status object.
        // The actual size of the message is not really needed, so we only read it when not in NDEBUG mode.
        // If the correct size is available, then the MPI_Unpack routine automatically checks whether the buffer position
        // exceeds the length of the sent message (Comm::MCW_buffer_pos > Comm::MCW_received_bytes) and throws an error
        // like this:
        //   *** An error occurred in MPI_Unpack
        //   *** on communicator MPI_COMM_WORLD
        //   *** MPI_ERR_TRUNCATE: message truncated
        MPI_Get_count(&status, MPI_PACKED, &Comm::MCW_received_bytes);
#else
        // In NDEBUG mode the actual size of the sent message is not read. Instead, the size of the buffer array is
        // used. Hence, it may happen that the MPI_Unpack routine reads rubbish from the buffer without throwing an
        // truncation error.
        Comm::MCW_received_bytes = Comm::MCW_BUFFERSIZE;
#endif

        // reset buffer content pointer
        Comm::MCW_buffer_pos = 0;

        // read first integer which is the ID of the message
        Comm::read(id);

//        // debug output
//        std::cout << "Master " << Process::rank << " received message with ID "
//                  << id << " and size " << Comm::MCW_received_bytes << "." << std::endl;

        switch(id)
        {
          case(ServiceIDs::LOG_RECEIVE):
          {
            Logger::receive();
            break;
          }
          /*case(ServiceIDs::LOG_RECEIVE_ARRAY):
          {
            Logger::receive_array();
            break;
          }*/
          case(ServiceIDs::MASTER_FINISH_SERVICE):
          {
            _finish_service = true;
            break;
          }
          default:
          {
            abort_mpi("ServiceID not found!");
          }
        }

        // exit the loop if some stop signal occured
        if(_finish_service)
        {
          break;
        }
      }
      std::cout << "Master finished service." << std::endl;
    }
  }; // class Master
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_FOUNDATION_MASTER_HPP

#pragma once
#ifndef KERNEL_MASTER_HPP
#define KERNEL_MASTER_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/service_ids.hpp>
#include <kernel/process_group.hpp>
#include <kernel/logger.hpp>


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

  MPI_Status status;

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
  /// dummy wait function
  void wait()
  {
    for(int i(0) ; i<1 ; ++i)
    {
      sleep(1.0);
      std::cout << "Master process with world rank " << Process::rank <<" is waiting..." << std::endl;
    }
  }

  /// infinite service loop waiting for messages
  void service()
  {
    std::cout << "Master at your service... " << std::endl;
    sleep(2.0);

    ServiceIDs::service_id id;
    while(true)
    {
      // receive messages with any tag and from any source
      // (Warning! Don't replace the status object by MPI_STATUS_IGNORE! It is needed in the following call
      // of MPI_get_count().)
      MPI_Recv(MPIUtils::buffer, MPIUtils::BUFFERSIZE_BYTES, MPI_PACKED, MPI_ANY_SOURCE,
               MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      // read the size of the received message (in bytes) from the status object
      MPI_Get_count(&status, MPI_PACKED, &MPIUtils::received_bytes);
// COMMENT_HILMAR:
// Is the actual size of the message really needed? Or can it be replaced by MPIUtils::BUFFERSIZE_BYTES in the
// following calls of MPI_Unpack? I think we only lose the functionality that the MPI_Unpack routine automatically
// checks whether the buffer position exceeds the length of the sent message
// (MPIUtils::buffer_pos > MPIUtils::received_bytes). But this error would be so critical that it will crash the
// program anyway, hence this check is not really necessary, is it? (BTW: Does the routine *really* check? I think so,
// but actually I'm not sure...)

      // reset buffer content pointer
      MPIUtils::buffer_pos = 0;

      // read first integer which is the ID of the message
      MPI_Unpack(MPIUtils::buffer, MPIUtils::received_bytes, &MPIUtils::buffer_pos, &id, 1,
                 MPI_INTEGER, MPI_COMM_WORLD);

//      // debug output
//      std::cout << "Master " << Process::rank << " received message with ID "
//                << id << " and size " << MPIUtils::received_bytes << "." << std::endl;

      switch(id)
      {
        case(ServiceIDs::LOG_RECEIVE):
        {
          Logger::receive();
          break;
        }
        case(ServiceIDs::LOG_RECEIVE_ARRAY):
        {
          Logger::receive_array();
          break;
        }
        case(ServiceIDs::MASTER_FINISH_SERVICE):
        {
          _finish_service = true;
          break;
        }
        default:
        {
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

#endif // guard KERNEL_MASTER_HPP

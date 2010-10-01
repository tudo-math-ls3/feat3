#pragma once
#ifndef KERNEL_COMM_HPP
#define KERNEL_COMM_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/service_ids.hpp>


/**
* \brief class providing basic communication functions
*
* \author Hilmar Wobker
*/
class Comm
{

private:

  /// MPI status object needed in MPI_Recv call
  MPI_Status status;

public:

  /* *****************
  * member variables *
  *******************/
  /// global buffer used for MPI_COMM_WORLD (MCW) communication
  static char* MCW_buffer;
  /// size of buffer Comm::MCW_buffer in bytes
  static int MCW_BUFFERSIZE;
  /// current position in buffer Comm::MCW_buffer
  static int MCW_buffer_pos;
  /// current size of buffer Comm::MCW_buffer
  static int MCW_received_bytes;

  /* *************************
  * constructor & destructor *
  ***************************/

  /* *****************
  * member functions *
  *******************/

  /**
  * \brief init a new MPI_COMM_WORLD message
  *
  * This function resets the MPI_COMM_WORLD buffer #MCW_buffer and writes the message ID to its first position.
  *
  * \param[in] service_id
  * ID of the message, where only values of the enumeration ServiceIDs::service_id are accepted.
  */
  static void init_msg(ServiceIDs::service_id service_id)
  {
    // reset buffer
    Comm::MCW_buffer_pos = 0;
    // write the message id to the buffer
    int mpi_error_code = MPI_Pack(&service_id, 1, MPI_INTEGER, Comm::MCW_buffer,
                                  Comm::MCW_BUFFERSIZE, &Comm::MCW_buffer_pos, MPI_COMM_WORLD);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
  }

}; // class Comm

// COMMENT_HILMAR: JUST TEMPORARILY
// initialisation of static members
// COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
int Comm::MCW_BUFFERSIZE = 4194304;
char* Comm::MCW_buffer = new char[MCW_BUFFERSIZE];
int Comm::MCW_buffer_pos = 0;
int Comm::MCW_received_bytes = 0;
// COMMENT_HILMAR: JUST TEMPORARILY

#endif // guard KERNEL_COMM_HPP

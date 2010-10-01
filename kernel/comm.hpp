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
  /// size of buffer #buffer in bytes
  static int MCW_BUFFERSIZE;
  /// current position in buffer #buffer
  static int MCW_buffer_pos;
  /// current size of buffer #buffer
  static int MCW_received_bytes;

  /* *************************
  * constructor & destructor *
  ***************************/
  /* *****************
  * member functions *
  *******************/

  // init message buffer
//  static void init_msg(int id, MPI_Comm& comm, void* buffer, int buffersize_bytes, int pos_bytes)
//  {
//    int mpi_error_code = mpi_pack(id, 1, MPI_INTEGER, buffer, buffersize_bytes, pos_bytes, comm)
//  }

}; // class Communication

// COMMENT_HILMAR: JUST TEMPORARILY
// initialisation of static members
// COMMENT_HILMAR: Use some arbitrary size for the time being. This has to be parameterised somehow...
int Comm::MCW_BUFFERSIZE = 4194304;
char* Comm::MCW_buffer = new char[MCW_BUFFERSIZE];
int Comm::MCW_buffer_pos = 0;
int Comm::MCW_received_bytes = 0;
// COMMENT_HILMAR: JUST TEMPORARILY

#endif // guard KERNEL_COMM_HPP

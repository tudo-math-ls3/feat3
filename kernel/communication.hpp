#pragma once
#ifndef KERNEL_COMMUNICATION_HPP
#define KERNEL_COMMUNICATION_HPP 1

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
class Communication
{

private:

  /* *****************
  * member variables *
  *******************/
  /// MPI status object needed in MPI_Recv call
  MPI_Status _status;

public:

  /* *************************
  * constructor & destructor *
  ***************************/
  /* *****************
  * member functions *
  *******************/
}; // class Communication

#endif // guard KERNEL_COMMUNICATION_HPP

#pragma once
#ifndef KERNEL_MASTER_HPP
#define KERNEL_MASTER_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>


/**
* \brief class defining the master
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*/
class Master
{

private:
  /* *****************
  * member variables *
  *******************/
  /**
  * \brief constructor
  */
  bool _service_finish;

public:

  /* *************
  * constructors *
  ***************/
  /**
  * \brief constructor
  */
  Master()
   : _service_finish(false);
  {
  }

  /* *****************
  * member functions *
  *******************/
  /**
  * \brief dummy wait function
  */
  void wait()
  {
    for (int i(0) ; i<1 ; ++i)
    {
      sleep(1.0);
      std::cout << "Master process with world rank " << Process::rank <<" is waiting..." << std::endl;
    }
  }

  /**
  * \brief infinite service loop waiting for messages
  */
  void service()
  {
    while (true)
    {
      mpi_recv(buf, PAR_BUFFERSIZE_BYTES, MPI_PACKED, MPI_ANY_SOURCE, &
               MPI_ANY_TAG, MPI_COMM_WORLD, statuss)

      // exit the loop if there is some stop signal
      if (_bmasterFinish)
      {
        break;
      }
    }
  }
}; // class Master

#endif // guard KERNEL_MASTER_HPP

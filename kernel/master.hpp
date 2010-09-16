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
  /* *****************
   * private members *
   *******************/
  private:

  /* ****************
   * public members *
   ******************/
  public:
    /* **************
     * constructors *
     ****************/
    /**
    * \brief constructor
    */
    Master()
    {
    }

    /* ******************
     * member functions *
     ********************/
    // dummy function
    void wait()
    {
      for (int i(0) ; i<1 ; ++i)
      {
        sleep(1.0);
        std::cout << "Master process with world rank " << Process::rank <<" is waiting..." << std::endl;
      }
    }
}; // class Master

#endif // guard KERNEL_MASTER_HPP

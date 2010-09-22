#pragma once
#ifndef KERNEL_PROCESS_HPP
#define KERNEL_PROCESS_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

// includes, Feast
#include <kernel/base_header.hpp>


/**
* \brief base class encapsulating an MPI process
*
* For each MPI process spawned at program start, one Process object will created on this MPI process. It will live
* throughout the program. Several classes store a pointer to this process object, i.e. an object of such a class has
* the property to live on this process.
* COMMENT_HILMAR: Maybe it makes sense to define this class as static or as a Singleton... we'll see.
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*/
class Process
{

private:


public:

  /* *****************
  * member variables *
  *******************/
  /**
  * \brief rank of the process within MPI_COMM_WORLD, set once via constructor and never changed again
  */
  static int rank;

  /**
  * \brief rank of the master process within MPI_COMM_WORLD, set once via constructor and never changed again
  *
  * Every process has to know the rank of the master process in order to trigger screen output (although this will
  * be mainly done by some coordinator processes).
  */
  static int rank_master;

  /**
  * \brief flag whether this process is the master process
  */
  static bool is_master;

  /* *************
  * constructors *
  ***************/
  /**
  * \brief empty constructor
  */
  Process();

}; // class Process

// define static member variables
int Process::rank = MPI_PROC_NULL;
int Process::rank_master = MPI_PROC_NULL;
bool Process::is_master = false;



///**
//* \brief base class encapsulating an MPI process
//*
//* For each MPI process spawned at program start, one Process object will created on this MPI process. It will live
//* throughout the program. Several classes store a pointer to this process object, i.e. an object of such a class has
//* the property to live on this process.
//* COMMENT_HILMAR: Maybe it makes sense to define this class as static or as a Singleton... we'll see.
//*
//* @author Hilmar Wobker
//* @author Dominik Goeddeke
//*/
//class Process
//{
//  /* *****************
//   * private members *
//   *******************/
//  private:
//    /* ******************
//     * member variables *
//     ********************/
//    /**
//    * \brief rank of the process within MPI_COMM_WORLD, set once via constructor and never changed again
//    */
//    const int _rank;
//    /**
//    * \brief rank of the master process within MPI_COMM_WORLD, set once via constructor and never changed again
//    *
//    * Every process has to know the rank of the master process in order to trigger screen output (although this will
//    * be mainly done by some coordinator processes).
//    */
//    const int _rank_master;
//
//  /* ****************
//   * public members *
//   ******************/
//  public:
//    /* **************
//     * constructors *
//     ****************/
//    /**
//    * \brief constructor requiring two parameters
//    */
//    Process(
//      const int rank,
//      const int rank_master)
//      : _rank(rank),
//        _rank_master(rank_master)
//    {
//    }
//
//    /* *******************
//     * getters & setters *
//     *********************/
//    /**
//    * \brief getter for the MPI_COMM_WORLD rank
//    */
//    inline int rank() const
//    {
//      return _rank;
//    }
//
//    /**
//    * \brief getter for the MPI_COMM_WORLD rank of the master process
//    */
//    inline int rank_master() const
//    {
//      return _rank_master;
//    }
//}; // class Process


///**
//* \brief class defining a remote process
//*
//* @author Hilmar Wobker
//*/
//class RemoteProcess
//  : public Process
//{
//  private:
//    /**
//    * \brief pointer to the MPI communicator shared by the RemoteProcess and this process
//    */
//    MPI_Comm* const _comm;
//
//    /**
//    * \brief rank of the remote process with respect to the shared communicator
//    */
//    const int _rank;
//
//  /* ****************
//   * public members *
//   ******************/
//  public:
//    /* **************
//     * constructors *
//     ****************/
//    /**
//    * \brief constructor
//    */
//    RemoteProcess(
//      const int rank_world,
//      const int rank_master,
//      MPI_Comm* const comm,
//      const int rank)
//      : Process(rank_world, rank_master),
//        _comm(comm),
//        _rank(rank)
//    {
//    }
//}; // class RemoteProcess

#endif // guard KERNEL_PROCESS_HPP

// Hilmar: BRAL: In dieser Datei werden erstmal einige Klassen gesammelt. Es macht spaeter bestimmt Sinn, einige davon
// auszulagern.
// Dom: FORTYTWO: zustimm

#pragma once
#ifndef KERNEL_PROCESS_HPP
#define KERNEL_PROCESS_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>


/**
* \brief Abstract base class encapsulating an MPI process
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*/
class Process
{
  /* *******************
   * protected members *
   *********************/
  protected:
    /* ******************
     * member variables *
     ********************/
    /**
    * \brief rank of the process within MPI_COMM_WORLD, set once via constructor and never changed again
    */
    const int _rank_world;
    /**
    * \brief rank of the master process within MPI_COMM_WORLD, set once via constructor and never changed again
    *
    * Every process has to know the rank of the master process in order to trigger screen output (although this will
    * be mainly done by some coordinator processes).
    */
    const int _rank_master;

  /* ****************
   * public members *
   ******************/
  public:
    /* **************
     * constructors *
     ****************/
    /**
    * \brief constructor requiring two parameters
    */
    Process(
      const int rank_world,
      const int rank_master)
      : _rank_world(rank_world),
        _rank_master(rank_master)
    {
    }

    /* *******************
     * getters & setters *
     *********************/
    /**
    * \brief getter for the MPI_COMM_WORLD rank
    */
    inline int rank_world() const
    {
      return _rank_world;
    }
}; // class Process


/**
* \brief class defining the master process
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*/
class Master
  : public Process
{
  /* ****************
   * public members *
   ******************/
  public:
    /* **************
     * constructors *
     ****************/
    /**
    * \brief constructor requiring one parameter
    */
    Master(const int rank_world)
      : Process(rank_world, rank_world)
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
        std::cout << "Master process with world rank " << _rank_world <<" is waiting..." << std::endl;
      }
    }
}; // class Master


/**
* \brief class defining defining a remote worker process
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*/
class RemoteWorker
  : public Process
{
  private:
    // rank of the remote worker within the communicator of the corresponding work group
    const int _rank_work_group;

  /* ****************
   * public members *
   ******************/
  public:
    /* **************
     * constructors *
     ****************/
    /**
    * \brief constructor requiring four parameters
    */
    RemoteWorker(
      const int rank_world,
      const int rank_master,
      const int rank_work_group)
      : Process(rank_world, rank_master),
        _rank_work_group(rank_work_group)
    {
    }
}; // class RemoteWorker


/**
* \brief class defining a worker process
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*/
class Worker
  : public Process
{
  private:
    // rank of this process within the local communicator
    const int _rank_work_group;

    // communicator used by the work group this worker belongs to
    const MPI_Comm _comm_work_group;

    // array of workers treating the neighbour subdomains
    RemoteWorker** _neighbours;

    // array of workers in the "next finer" work group this worker has to communicate with
    RemoteWorker** _finer; // besseren Namen finden!

    // worker in the "next coarser" work group this worker has to communicate with (this is always only one!)
    RemoteWorker* _coarser; // besseren Namen finden!

/*
 * BRAL: zu _finer und _coarser: Dasselbe Beispiel wie in load_balancer.hpp:

 * distribution of submeshes to workers A-G on different levels:
 *
 *     ---------------      ---------------      ---------------
 *     |             |      |      |      |      |      |      |
 *     |             |      |      |      |      |  D   |  G   |
 *     |             |      |      |      |      |      |      |
 *     |      A      |      |  B   |  C   |      ---------------
 *     |             |      |      |      |      |      |      |
 *     |             |      |      |      |      |  E   |  F   |
 *     |             |      |      |      |      |      |      |
 *     ---------------      ---------------      ---------------
 *       level 0               level 1              levels 2-L
 *
 *     * case a:
 *       process group rank:  0  1  2  3
 *             work group 2:  D  E  F  G         (four workers for the problems on level 2-L)
 *             work group 1:  B     C            (two workers for the problem on level 1)
 *             work group 0:  A                  (one worker for the coarse mesh problem on level 0)
 *
 *       Communications:
 *       A <--> B (internal, rank 0) A <--> C (external, ranks 0+2)
 *       B <--> D (internal, rank 0) B <--> E (external, ranks 0+1)
 *       C <--> F (internal, rank 2) C <--> G (external, ranks 2+3)
 *
 *     * case b:
 *       process group rank:  0  1  2  3  4
 *             work group 2:     D  E  F  G
 *             work group 1:     B     C
 *             work group 0:  A
 *
 *       Communications:
 *       A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+3)
 *       B <--> D (internal, rank 1) B <--> E (external, ranks 1+2)
 *       C <--> F (internal, rank 3) C <--> G (external, ranks 3+4)
 *
 *     * case c:
 *       process group rank:  0  1  2  3  4  5  6
 *             work group 2:           D  E  F  G
 *             work group 1:     B  C
 *             work group 0:  A
 *
 *       Communications:
 *       A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+2)
 *       B <--> D (external, ranks 1+3) B <--> E (external, ranks 1+4)
 *
 *    Zunaechst mal zur Begrifflichkeit: In case a sind es tatsaechlich 3 Worker (A,B,D) auf Process 0 (und nicht etwa
 *    *ein* worker, der zu drei WorkGroups gehoert). Das heisst: ein Worker gehoert zu genau *einer* WorkGroup.
 *
 *    Die _finer RemoteWorker zu B in case a sind D und E, die zu zu C sind F und G. Der _coarser RemoteWorker zu
 *    B und C ist A. Man beachte, dass A, B und D auf dem selben Prozess "leben" (es handelt sich also eigentlich gar
 *    nicht um *Remote*Worker), waehrend B und E auf verschiedenen Prozessen leben. Dies muss noch irgendwie
 *    verwurschtelt werden.
 */

// BRAL: Nicht noetig! Siehe Kommentar zur class Coordinator
//  // rank of the coordinator process in this work group
//  // (MPI_PROC_NULL if the process is the coordinator)
//  int _rank_coordinator;
//    or
//  RemoteWorker _coordinator;

  public:
    /* **************
     * constructors *
     ****************/
    /**
    * \brief constructor requiring four parameters
    */
    Worker(
      const int rank_world,
      const int rank_master,
      const int rank_work_group,
      const MPI_Comm comm_work_group)
      : Process(rank_world, rank_master),
        _rank_work_group(rank_work_group),
        _comm_work_group(comm_work_group)
    {
    }

}; // class worker


// BRAL: In einer Coordinator Klasse sehe ich im Moment noch keinen Sinn. Es ist nichts anderes als der Worker mit rank
// 0 im entsprechenden Communicator. Seine einzige Extra-Aufgabe ist das Senden von Nachrichten an den Master. Man wird
// im Laufe der Implementierung sehen, ob es Sinn macht, daraus vielleicht doch 'ne Extra-Klasse zu machen.
// // class defining a coordinator process
//class Coordinator
//  : public Worker
//{
//
//};

#endif // guard KERNEL_PROCESS_HPP

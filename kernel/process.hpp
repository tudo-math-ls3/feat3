// BRAL: In dieser Datei werden erstmal einige Klassen gesammelt. Es macht später bestimmt Sinn, einige davon
// auszulagern.

#pragma once
#ifndef KERNEL_PROCESS_HPP
#define KERNEL_PROCESS_HPP 1

#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <kernel/base_header.hpp>

using namespace std;

// class defining an MPI process
// BRAL: grundlegende
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
     * \brief rank of the process within MPI_COMM_WORLD
     */
    int _rank_world;
    /**
     * \brief rank of the master process within MPI_COMM_WORLD
     *
     * Every process has to know the rank of the master process in order to trigger screen output (although this will
     * be mainly done by some coordinator processes).
     */
    //
    int _rank_master;

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
      int& rank_world,
      int& rank_master)
      : _rank_world(rank_world),
        _rank_master(rank_master)
    {
    }

    /**
     * \brief getter for the MPI_COMM_WORLD rank
     */
    int get_rank_world()
    {
      return _rank_world;
    }
};

// class defining the master process
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
    Master(int& rank_world)
      : Process(rank_world, rank_world)
    {
    }

    // dummy routine
    void wait()
    {
      while (true)
      {
        sleep(3.3);
        cout << "Master process with world rank " << _rank_world <<" is waiting..." << endl;
      }
    }
};

/* class defining a group process
 * BRAL: Was ist ein GroupProcess?
 * Am Anfang des Programms werden die verfügbaren Prozesse vom User in Gruppen eingeteilt (z.B. für Multiphysics-
 * Aufgaben). Zu jeder Gruppe gehört ein load balancer. Da erst der load balancer entscheidet, wie er "seine" Prozesse
 * einteilt, machen diese Prozesse also erstmal nix außer darauf zu warten, eingeteilt zu werden. Dieses Zwischenstadium
 * wird durch die Klasse GroupProcess beschrieben. Jeder Prozess, der also nicht load balancer oder master ist, wird
 * also zunächst mal als GroupProcess angesehen und in eine Warteschleife versetzt.
 * Wenn der load balancer dann entschieden hat, wie er seine Prozesse einteilen will, schickt er den GroupProcesses
 * entsprechende Nachrichten, sie sollen entsprechende Worker Prozesse erstellen.
 */
class GroupProcess
  : public Process
{
  /* *****************
   * private members *
   *******************/
  private:
    /* ******************
     * member variables *
     ********************/

    // communicator shared by the load balancer and all processes of the corresponding group
    MPI_Comm _comm_local;

    // rank of this process within the local communicator
    int _rank_local;

    // rank of the load balancer responsible for this GroupProcess (with respect to the local communicator)
    int _rank_load_bal;

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
    GroupProcess(
      int& rank_world,
      int& rank_master,
      int& rank_load_bal,
      MPI_Comm& comm_local,
      int& rank_local)
      : Process(rank_world, rank_master),
        _rank_load_bal(rank_load_bal),
        _comm_local(comm_local),
        _rank_local(rank_local)
    {
    }

    // dummy routine
    void wait()
    {
      while (true)
      {
        sleep(2);
        cout << "GroupProcess with world rank " << _rank_world <<" is waiting..." << endl;
      }
    }
};







// class defining a remote worker process
class RemoteWorker
  : public Process
{
  private:
    // rank of the remote worker within the corresponding communicator
    int _rank_local;
};



// class defining a worker process
class Worker
  : public Process
{
  private:
    // communicator used by the work group this worker belongs to
    MPI_Comm _comm_work_group;

    // rank of this process within the local communicator
    int _rank_work_group;

    // workers treating the neighbour subdomains
    RemoteWorker* _neighbours;

    // array of workers in the "next finer" work group this worker has to communicate with
    RemoteWorker* _finer; // besseren Namen finden!

    // worker in the "next coarser" work group this worker has to communicate with (this is always only one!)
    RemoteWorker _coarser; // besseren Namen finden!

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
 *    Zunächst mal zur Begrifflichkeit: In case a sind es tatsächlich 3 Worker (A,B,D) auf Process 0 (und nicht etwa
 *    *ein* worker, der zu drei WorkGroups gehört). Das heißt: ein Worker gehört zu genau *einer* WorkGroup.
 *
 *    Die _finer RemoteWorker zu B in case a sind D und E, die zu zu C sind F und G. Der _coarser RemoteWorker zu
 *    B und C ist A. Man beachte, dass A, B und D auf dem selben Prozess "leben" (es handelt sich also eigentlich gar
 *    nicht um *Remote*Worker), während B und E auf verschiedenen Prozessen leben. Dies muss noch irgendwie
 *    verwurschtelt werden.
 */

    // pointer to the GroupProcess, this worker corresponds to
    // (note, that more than one Worker can "belong" to one GroupProcess)
    // BRAL: über diesesn pointer lässt sich Verbindung zum load balancer aufnehmen, falls nötig
    GroupProcess* _group_process;

// BRAL: Nicht noetig! Siehe Kommentar zur class Coordinator
//  // rank of the coordinator process in this work group
//  // (MPI_PROC_NULL if the process is the coordinator)
//  int _rank_coordinator;
//    or
//  RemoteWorker _coordinator;
};


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

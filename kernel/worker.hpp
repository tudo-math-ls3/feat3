#pragma once
#ifndef KERNEL_WORKER_HPP
#define KERNEL_WORKER_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>

/**
* \brief class defining a remote worker
*
* A RemoteWorker is simply a Worker that lives on a remote process. So, a WorkGroup object managing n processes
* typically consists of one Worker object and n-1 RemoteWorker objects.
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
class RemoteWorker
{
private:
  /// pointer to the MPI communicator to the RemoteWorker
  MPI_Comm* const _comm;

  /// rank of the remote worker with respect to the communicator RemoteWorker::_comm
  const int _rank;

public:
  /* *************************
  * constructor & destructor *
  ***************************/
  /// constructor
  RemoteWorker(
    MPI_Comm* const comm,
    const int rank)
    : _comm(comm),
      _rank(rank)
  {
  }
}; // class RemoteWorker


// /// forward declaration of WorkGroup class
// class WorkGroup;

/**
* \brief class defining a worker
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
class Worker
{

private:

  /* *****************
  * member variables *
  *******************/
  /// communicator of the work group this worker belongs to
  MPI_Comm _comm_work_group;

  /// this worker's rank w.r.t. the communicator Worker::_comm_work_group
  int _rank_work_group;

  /// communicator of the process group this worker belongs to
  MPI_Comm _comm_process_group;

  /// this worker's rank w.r.t. the communicator Worker::_comm_process_group
  int _rank_process_group;

// COMMENT_HILMAR, 22.9.2010:
// Schick waeren natuerlich pointer auf die entsprechende WorkGroup bzw. die entsprechende ProcessGroup anstatt der
// obigen vier Variablen, aber damit faengt man sich eine zyklische Abhaengigkeit zwischen WorkGroup/ProcessGroup
// einerseits und Worker andererseits ein. "forward declaration" allein hilft nicht, man muesste die benoetigten
// Funktionen dann ausserhalb der Klasse definieren bzw. in eine Extra-Datei auslagern. Das ist mir im Moment zu
// umstaendlich. Mal sehen, ob das spaeter noetig bzw. sinnvoll sein wird...

  /**
  * \brief vector of work group ranks of workers treating neighbour subdomains
  *
  * Note that all these workers live on remote processes.
  */
  std::vector<int> _ranks_neighbours;

  /**
  * \brief vector of process group ranks of workers in the "next finer" work group this worker has to exchange
  *        data with
  *
  * Since these workers live in a different work group, this worker has to communicate with them via the process
  * group communicator. Note that one of these workers may have the same process group rank (i.e., live on the same
  * MPI_COMM_WORLD process) as this worker (hence, no MPI communication is necessary in this case).
  */
  std::vector<int> _ranks_finer; // COMMENT_HILMAR: find a better variable name!

  /**
  * \brief process group rank of the worker in the "next coarser" work group this worker has to exchange data with
  *        (this is always only one!)
  *
  * Since the worker lives in a different work group, this worker has to communicate with it via the process
  * group communicator. Note that one of the worker may have the same process group rank (i.e., live on the same
  * MPI_COMM_WORLD process) as this worker (hence, no MPI communication is necessary in this case).
  */
  int _rank_coarser; // COMMENT_HILMAR: find a better variable name!

// COMMENT_HILMAR, 22.9.2010:
// Die letzten drei Variablen koennte man auch ueber Pointer of Worker bzw. RemoteWorker-Objekte realisieren.
// Allerdings weiss man a priori nicht, ob es sich um einen Worker oder RemoteWorker handelt. Moechte man alle in
// einen Vektor packen, muesste man sie also aus einer "kuenstlichen" gemeinsamen Elternklasse ableiten, diese im
// Vektor abspeichern und dann beim Zugriff jeweils ueberpruefen, worum es sich nun handelt, um die richtige Klasse
// anzusprechen. Das erscheint mir im Moment etwas umstaendlich, deswegen erstmal die Realisierung ueber
// integer-Vektoren. Natuerlich muss man auch hier abfragen, ob der entsprechende Worker auf demselben Prozess lebt,
// denn MPI_send/recv mit source=destination ist unsafe. Dazu Zitat aus dem MPI2.2-Standard, Seite 31:
//   Source = destination is allowed, that is, a process can send a message to itself.
//   (However, it is unsafe to do so with the blocking send and receive operations described above,
//   since this may lead to deadlock. See Section 3.5.)
// Fuer den Fall, dass source und destination worker auf dem selben Prozess leben, muss kopiert anstatt kommuniziert
// werden. Dafuer muss natuerlich eine entsprechende Funktion vorgesehen werden, und die beiden muessen sich
// datenstruktur-technisch "kennen", also wird man wohl doch einen Pointer auf den entsprechenden Worker abspeichern
// muessen... na, mal sehen, wenn's an die Details geht.

/*
* COMMENT_HILMAR: zu _finer und _coarser: Dasselbe Beispiel wie in load_balancer.hpp:

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


public:

  /* *************************
  * constructor & destructor *
  ***************************/
  /// constructor
  Worker(
    MPI_Comm comm_work_group,
    int rank_work_group,
    MPI_Comm comm_process_group,
    int rank_process_group)
    : _comm_work_group(comm_work_group),
      _rank_work_group(rank_work_group),
      _comm_process_group(comm_process_group),
      _rank_process_group(rank_process_group)
  {
  }
}; // class worker

#endif // guard KERNEL_WORKER_HPP

/* GENERAL_REMARK_BY_HILMAR:
 * The description of this manager class especially sketches the general program flow.
 * See COMMENT_HILMAR throughout this file and issue_00037.txt.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_MANAGER_HPP
#define KERNEL_MANAGER_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/process_group.hpp>
#include <kernel/load_balancer.hpp>
#include <kernel/load_balancer_dedicated.hpp>
#include <kernel/manager_comp.hpp>
#include <kernel/manager_comp_coord.hpp>
#include <kernel/manager_comp_non_coord.hpp>

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

namespace FEAST
{
  /**
  * \brief class defining a process group manager
  *
  * The user creates one or more main process groups. These groups have nothing to do with each other. Every communication
  * between them has to be organised manually by the user (via the parent MPI communicator MPI_COMM_WORLD). More than
  * one process group is used in case the user wants to perform some multiphysics problem or similar where several
  * tasks can be performed (nearly) independent of each other. The user knows what each main process group
  * (distinguishable via its ID) should do, i.e., he calls something like
  \verbatim
  if(group_id == 0)
  {
    // solve Poisson equation
    ...
  }
  else
  {
    // brew coffee
    ...
  }
  \endverbatim
  * Each main process group is organised by a manager. It runs on all processes of the process group, which
  * eases work flow organisation significantly. This Manager class serves as 'interface' between user code and kernel
  * code. The routines under user control basically run on all processes of the main process group. In order to keep
  * user code as simple as possible the manager (and not the user) decides which process actually performs which task.
  * As an example, the computational mesh should only be read by one process (the coordinator process, see below). The
  * user only has to call
  \verbatim
     manager->read_mesh(...);
  \endverbatim
  * and then the manager performs something like
  \verbatim
     if(is_coordinator())
     {
       manager_coord->read_mesh(...);
     }
  \endverbatim
  * i.e., only on the coordinator process a corresponding read routine is called in some object manager_comp_coord
  * which lives on the coordinator process (see below). This is especially important when code is involved that calls
  * collective MPI routines (like MPI_Bcast(...)): The user just calls one general routine on all processes (without
  * any distinction who sends and who receives data), while the manager then guarantees that the correct routines
  * on sender and receiver side are called.
  *
  * All standard tasks within the simulation process should be organised in this way via this Manager class
  * (read mesh, define bilinear form, set boundary conditions, assembly matrices, start solver, perform visual output,
  * etc.). Of course, there may be special things the user wants to do and which can not be realised via the Manager
  * class (yet). This may require the user to do manual distinctions of the processes on his own. However, a 'standard'
  * user simulation code should do without ever explicitly asking 'are you the coordinator process?'.
  * (COMMENT_HILMAR: At least, that's the plan...)
  * (COMMENT_DOM: Maybe some sort of dummy callback routine is necessary. Also, all the infrastructure operations
  *  should be somehow steerable by the user, e.g., to decide if loadbalancing is performed every solver iteration,
  *  every solver call, or every timestep...)
  *
  * The processes of the group take different roles, and depending on that, different objects are created on the
  * processes.
  * - The coordinator process is the only one storing the complete computational mesh, reading config files etc. The
  *   coordinator has rank 0 within the main process group. The coordinator is not dedicated, but also serves
  *   as a compute process.
  * - The user can decide to use a dedicated load balancer process. It is the process with the largest rank in the
  *   main process group. COMMENT_DOM: I thought that's the master?
  * - All other processes are standard 'compute' processes.
  *
  * As 'compute process' we denote those processes which do the actual compute work (assembly, solving linear systems,
  * etc.).
  *
  * While the (optional) dedicated load balancer process is exclusively used for special load balancing tasks which
  * can be performed asynchroneously to the actual computations, the coordinator process also performs standard
  * compute tasks, i.e., it is also a compute process.
  *
  * Hence, all processes except the dedicated load balancer process build the group of compute processes. That is why
  * an extra process group is built that represents this group of compute processes. In case, there is no dedicated
  * load balancer, the main process group and the compute process group are identical.
  *
  * While this Manager class is responsible for all processes of the main process group (including the dedicated load
  * balancer), most (or even all) of the user instructions only concern the compute processes. Hence, there are
  * 'sub manager' classes, that are only active for these compute processes. These are:
  * - ManagerCompCoord: manager living on the coordinator of the compute process group. It is responsible for
  *   reading and distributing the mesh to the other compute processes, and for communication with the load balancer.
  *   The coordinator is part of at least one 'work group' (see below), i.e., it also performs compute tasks.
  * - ManagerCompNonCoord: manager living on all non-coordinator processes of the compute process group. It receives
  *   data from the coordinator (e.g., the relevant parts of the base mesh) and sends data back (e.g., statistics for
  *   the load balancer (see below).
  * - ManagerComp: Super class of the previous two which contains common data and functions.
  *
  * Usually, there are two matching routines in ManagerCompCoord and ManagerCompNonCoord, e.g., something like
  \verbatim
     send_data(...)
  \endverbatim
  * in ManagerCompCoord and
  \verbatim
     receive_data(...)
  \endverbatim
  * in ManagerCompNonCoord. The user calls a Manager function like
  \verbatim
     transfer_data()
  \endverbatim
  * (on all processes) which then does something like
  \verbatim
     if(is_coordinator())
     {
       manager_coord->send_data(...);
     }
     else
     {
       manager_non_coord->receive_data(...);
     }
  \endverbatim
  * This especially means, that the non-coordinator processes are always under user control. It's not like they are
  * in an infinite wait loop to be woken up by the coordinator for special tasks.
  *
  * On the coordinator process, a further object is created: a LoadBalancer object. It lives here because it needs
  * access to the global mesh/patch information. It is responsible for performing tasks that have to do with work
  * distribution / load balancing (process matrix patch (MP) statistics, compute partitioning, define work groups, ...).
  * It lives on the same process as the compute coordinator ManagerCompCoord since the two share some data, e.g.,
  * the base mesh (BM)/matrix patch mesh (MPM).
  *
  * Input for the load balancer is (probably among others) the connectivity of the MPM and MP statistics. Based on the
  * provided information, the load balancer distributes MPs to MPI processes, i.e., it creates process patches (PP).
  * When no balanced load distribution is possible, the load balancer might ask for a 'refined version' of a MP, i.e.,
  * a split of the MP into four MPs, for example, basing on the refinement of the corresponding base mesh cell. With
  * this, the load balancer may have better chances to equally distribute the work.
  *
  * In the case a multigrid solver is used, the coarse grid problem is defined as so called coarse grid matrix
  * patch (CGMP) and has to be treated in a special way by the load balancer. The standard case will be that the coarse
  * grid problem is solved on one process, i.e. one CGMP containing all base mesh cells. If this is not possible
  * (e.g., due to memory constraints), the load balancer asks for a refined version of CGMP. In contrast to the case
  * when the refined version of a standard MP is needed (which can be provided by the base mesh), the refined version
  * of the CGMP should be created by the load balancer itself. When, for example, it turns out that at least 4 processes
  * are necessary to solve the coarse grid problem (due to memory consumption), partitioning the whole MPM into four
  * sub-MPMs is in itself a load balancing task. The load balancer knows the MPM, hence this procedure should be
  * possible. (Note that basis for distributing the coarse grid problem is not the BM, but the MPM!) When the four
  * CGMPs are set up, they are distributed to four processes. Here, it does not make sense to schedule, e.g., two CGMPs
  * to one process (because we just created four since we need to distribute the original CGMP to four processes). So,
  * if you want to talk about 'coarse grid process patches' (CGPPs), then we always have CGMP = CGPP.
  *
  * In the case multigrid is used for solving, the load balancer defines so called WorkGroup objects, one for the
  * coarse grid problem and (at least) one for the finer grid problems. (In future, we eventually want to use different
  * number of processes (i.e., different work groups) on different levels, hence there may be more than one work group
  * for the finer grid problems.) These work groups are again process groups (with their own MPI communicator)
  * containing the processes of the compute process group responsible for the specific problem. (This is necessary
  * since collective communication is used for several computations during the computation.) However, there also has to
  * be an efficient way of communication between the coordinator of the (parent) compute process group and the work
  * groups. If this coordinator is part of the work group already, this is no problem: then the communicator of the work
  * group can be used. But when the coordinator is not part of the work group, we additionally build an 'extended work
  * group' (WorkGroupExt) which contains the compute processes of the work group and the coordinator of the parent
  * compute process group (see description of class WorkGroupExt for details). Note that one physical process can be
  * part of several work groups, i.e., work groups are not disjunct. (One process can solve a part of the fine grid
  * problem and also (a part of) the coarse grid problem.)
  *
  * While the load balancer defines (extended) work groups, the coordinator of the compute process group actually sets
  * them up and, for example, transfers to each worker process the relevant local part of the global communication
  * graph based on the process patch mesh (PPM). Each worker process then creates its local graph structure and calls
  * MPI_Dist_graph_create(...) to build the new MPI process topology in a distributed fashion.
  *
  * When work groups are set up, the coordinator of the compute process group also has to define for each worker within
  * a work group (using the information provided by the load balancer) which workers from other work groups it has to
  * communicate with. E.g., the fine mesh work group has to send the restricted defect vector to the
  * coarse mesh work group, while the coarse mesh work group has to send the coarse mesh correction to the fine mesh
  * work group. For this communication we again set up extra process groups, so called InterLevelGroup objects (see
  * class InterLevelGroup for details.)
  * Note, that two such communicating workers live either on the same process (internal communication = copy) or
  * on different processes (external communication = MPI send/recv). (See example below.)

  * The compute coordinator ManagerCompCoord can be seen as interface between the load balancer object and the compute
  * process group. Since load balancing and computation are performed after each other in an alternating fashion (e.g.,
  * after each time step, load balancing is initiated), it should not be a problem that load balancer and the compute
  * coordinater ManagerCompCoord live on the same process. Communication between the two is always triggered from the
  * manager object, i.e., the manager knows the load balancer, but the load balancer does not know the manager.
  *
  * COMMENT_HILMAR: It may turn out that it makes more sense, to let the LoadBalancer object live on all processes of
  *   the compute process group (maybe divided in LoadBalancerCoord and LoadBalancerNonCoord). Since load balancing
  *   details are not clear yet, I'm not sure what is the better way. Currently, I like the idea better that the
  *   LoadBalancer object exclusivly talks to the ManagerCompCoord, and only ManagerCompCoord talks to the
  *   non-coordinator processes.
  * COMMENT_DOM: A consequence of this alternative is that all process know the global mesh. Might be a problem to
  *   keep this info synchroneous, e.g., when things are split or merged. Also, statistics information must be
  *   allgathered (n:n communication).
  *
  * If a dedicated load balancer is used, there is an extra process which is not part of the compute process group
  * (actually the one process, which distinguishes the main process group from the compute process group). On this
  * process, an object of type LoadBalancerDedicated is created. This object performs 'special' load balancing tasks
  * which can be done asynchroneously to the actual computation (e.g., pre-calculating some possible scenarios based
  * on history data).
  *
  * COMMENT_HILMAR: It is not quite clear yet how communication to this object works. My current idea is that only the
  *   load balancer class is 'connected' to the dedicated load balancer class. I.e., data transfer etc. has to be
  *   organised by the load balancer class (which should have access to all necessary data since it lives on the
  *   coordinator process).
  *
  * COMMENT_HILMAR: In the current implementation the dedicated load balancer process is also under user control (in
  *   the sense that it is part of the main process group under control of the Manager object, and all these processes
  *   are forked to the user code). For the dedicated load balancer process, it might make sense to start an infinite
  *   wait loop, which can only be accessed by the LoadBalancer class...
  *   Since it is not clear yet, how all this is organised, I don't know what is the best way...
  *
  * In addition to the standard load balancer, there is also a 'numerical load balancer', which is responsible, e.g.,
  * for assessing base mesh cells w.r.t. their aspect ratio or for defining groups of BMCs that build one matrix patch.
  *
  * COMMENT_HILMAR: This numerical load balancer is not implemented at all yet. Probably, a corresponding object will
  *   also live on the coordinator process. As long as there is no numerical load balancer, one matrix patch will
  *   consist of exactly one base mesh cell (unless a temporary procedure for defining matrix patches is implemented).
  *
  * COMMENT_HILMAR: It might be more clever to also let the MPI implementation decide on which physical process the
  *   dedicated load balancer should reside (instead of pinning it to the last rank in the process group). To improve
  *   this is task of the ITMC.
  *
  *
  * FIRST EXAMPLE:
  * fine grid problems are solved by 4 processes, coarse grid problem by 2
  *
  * The matrix patch mesh (MPM) consists of 16 MPs. For the fine grid problems these are scheduled to 4 PPs (see
  * process patch mesh (PPM) in the second figure). The coarse grid problem is too large for 1 process, so the
  * original CGMP (which covered all 16 MPs) is split by the load balancer into 2 CGMPs (=CGPPs), each covering
  * 8 MPs.
  *
  * The letters A-G represent compute tasks to be scheduled to a process each. Note that processes A-G are not
  * necessarily disjunct, i.e., several of them can refer to the same physical process (see following example).
  *
  \verbatim
  ---------------         ---------------         ---------------
  |  |   |   |  |         |      |      |         |      |      |
  |------|------|         |  C   |  D   |         |      |      |
  |  |   |   |  |         |      |      |         |      |      |
  ---------------         ---------------         |  E   |  F   |
  |  |   |   |  |         |      |      |         |      |      |
  |------|------|         |  A   |  B   |         |      |      |
  |  |   |   |  |         |      |      |         |      |      |
  ---------------         ---------------         ---------------
        MPM               PPM fine grids        2 CGMPs coarse grid
  \endverbatim
  *
  * The user has inquired 9 processes and 2 main process groups, where he uses 2 processes for the second main process
  * group which computes something different. For the first process group he uses a dedicated load balancer. Of the
  * 5 remaining compute processes, one has to be used for fine and coarse grid problems. The following table shows the
  * distribution of the processes to different kind of process groups (pg), and which objects live on which
  * processes (obj). Numbers denote the local ranks within the corresponding process group, 'X' means, that an object
  * of this class is instantiated on this process. Letters in parentheses refer to the compute tasks in the figures,
  * an asterisk '*' means 'this is the root process of different collective MPI routines'.
  *
  \verbatim
  ---------------------------------------------------------------------------------------------------------
  | COMM_WORLD ranks                |  0      1      2      3      4      5    |    6      7    |    8    |
  |-------------------------------------------------------------------------------------------------------|
  | obj: Master                     |  -      -      -      -      -      -    |    -      -    |    X    |
  | pg:  main process group 0       |  0      1      2      3      4      5    |    -      -    |    -    |
  | pg:  main process group 1       |  -      -      -      -      -      -    |    0      1    |    -    |
  |-------------------------------------------------------------------------------------------------------|
  | pg:  compute process group 0    |  0*     1      2      3      4      -    |    -      -    |    -    |
  | pg:  compute process group 1    |  -      -      -      -      -      -    |    0*     1    |    -    |
  | obj: Manager                    |  X      X      X      X      X      X    |    X      X    |    -    |
  | obj: ManagerCompCoord           |  X      -      -      -      -      -    |    X      -    |    -    |
  | obj: ManagerCompNonCoord        |  -      X      X      X      X      -    |    -      X    |    -    |
  | obj: LoadBalancer               |  X      -      -      -      -      -    |    X      -    |    -    |
  | obj: LoadBalancerDedicated      |  -      -      -      -      -      X    |    -      -    |    -    |
  |-------------------------------------------------------------------------------------------------------|
  | pg:  ext. work group for FGP    |  0*(A)  1(B)   2(C)   3(D)   -      -    |    ?      ?    |    -    |
  | pg:  work group for FGP         |  0*(A)  1(B)   2(C)   3(D)   -      -    |    ?      ?    |    -    |
  | obj: WorkGroupExt (FGP)         |  X      X      X      X      -      -    |    ?      ?    |    -    |
  | obj: WorkGroup (FGP)            |  X      X      X      X      -      -    |    ?      ?    |    -    |
  |-------------------------------------------------------------------------------------------------------|
  | pg:  ext. work group for CGP    |  0*     -      -      1(E)   2(F)   -    |    ?      ?    |    -    |
  | pg:  work group for CGP         |  -      -      -      0*(E)  1(F)   -    |    ?      ?    |    -    |
  | obj: WorkGroupExt (CGP)         |  X      -      -      X      X      -    |    ?      ?    |    -    |
  | obj: WorkGroup (CGP)            |  -      -      -      X      X      -    |    ?      ?    |    -    |
  |-------------------------------------------------------------------------------------------------------|
  | pg:  inter level  group 'left'  |  0(A)   -      1(C)   2*(E)  -      -    |    ?      ?    |    -    |
  | obj: InterLevelGroup            |  X             X      X      -      -    |    ?      ?    |    -    |
  | pg:  inter level group 'right'  |  -      0(B)   -      1(D)   2*(F)  -    |    ?      ?    |    -    |
  | obj: InterLevelGroup            |  -      X      -      X      X      -    |    ?      ?    |    -    |
  ---------------------------------------------------------------------------------------------------------
  \endverbatim
  *
  * Usually, the coordinator of a process group (rank 0) is also the root process for collective MPI routines (e.g., for
  * logging or collecting statistics). For inter level groups, however, the process stemming from the coarse grid work
  * group is the root since here we always have a 1:n communication pattern (1 in the coarse grid work group, n in the
  * fine grid work group). So, coordinator and root don't have to be the same process in inter level groups.
  *
  * Now imagine that COMM_WORLD process 3 treats task F and process 4 treats task E and look at the work groups and
  * inter level groups:
  \verbatim
  ---------------------------------------------------------------------------------------------------------
  | COMM_WORLD ranks                |  0      1      2      3      4      5    |    6      7    |    8    |
  |-------------------------------------------------------------------------------------------------------|
  | pg:  ext. work group for FGP    |  0*(A)  1(B)   2(C)   3(D)   -      -    |    ?      ?    |    -    |
  | pg:  work group for FGP         |  0*(A)  1(B)   2(C)   3(D)   -      -    |    ?      ?    |    -    |
  | obj: WorkGroupExt (FGP)         |  X      X      X      X      -      -    |    ?      ?    |    -    |
  | obj: WorkGroup (FGP)            |  X      X      X      X      -      -    |    ?      ?    |    -    |
  |-------------------------------------------------------------------------------------------------------|
  | pg:  ext. work group for CGP    |  0*     -      -      1(F)   2(E)   -    |    ?      ?    |    -    |
  | pg:  work group for CGP         |  -      -      -      0*(F)  1(E)   -    |    ?      ?    |    -    |
  | obj: WorkGroupExt (CGP)         |  X      -      -      X      X      -    |    ?      ?    |    -    |
  | obj: WorkGroup (CGP)            |  -      -      -      X      X      -    |    ?      ?    |    -    |
  |-------------------------------------------------------------------------------------------------------|
  | pg:  inter level  group 'left'  |  0(A)   -      1(C)   -      2*(E)  -    |    ?      ?    |    -    |
  | obj: InterLevelGroup            |  X             X      X      -      -    |    ?      ?    |    -    |
  | pg:  inter level group 'right'  |  -      0(B)   -     1*(D+F) -      -    |    ?      ?    |    -    |
  | obj: InterLevelGroup            |  -      X      -      X      -      -    |    ?      ?    |    -    |
  ---------------------------------------------------------------------------------------------------------
  \endverbatim
  * Now, tasks D and F are treated by the same process, namely 3. Hence, the inter level group for the right half
  * of the domain consists of only two processes. When it comes to 'communicating' data 'between' D and F, one has to
  * take this into account.
  *
  *
  * SECOND EXAMPLE:
  * scenarios for three work groups and varying number of processes (briefer representation than for the first example):
  *
  \verbatim
  ---------------      ---------------      ---------------
  |             |      |      |      |      |      |      |
  |             |      |      |      |      |  D   |  G   |
  |             |      |      |      |      |      |      |
  |      A      |      |  B   |  C   |      ---------------
  |             |      |      |      |      |      |      |
  |             |      |      |      |      |  E   |  F   |
  |             |      |      |      |      |      |      |
  ---------------      ---------------      ---------------
    level 0               level 1              levels 2-L
  \endverbatim
  *
  * -# case a, 4 physical processes:
  \verbatim
  comp. proc. group rank:   0    1    2    3
  WorkGroup 2:              D    E    F    G         (4 WorkGroup processes for the problems on lv. 2-L)
  WorkGroup 1:              B    -    C    -         (2 WorkGroup processes for the problem on lv. 1)
  WorkGroup 0:              A    -    -    -         (1 WorkGroup process for the coarse mesh problem on lv. 0)
  InterLevelGroup 1-2 le:  B+D   E    -    -         (lv. 1 worker B talks to lv. 2 work. D+E; B+D on same proc.)
  InterLevelGroup 1-2 ri:   -    -   C+F   G         (lv. 1 worker C talks to lv. 2 work. F+G; C+F on same proc.)
  InterLevelGroup 0-1:     A+B   -    C    -         (lv. 0 worker A talks to lv. 1 work. B+C; A+B on same proc.)
  \endverbatim
  *
  * -# case b, 5 physical processes::
  \verbatim
  comp. proc. group rank:   0    1    2    3    4
  WorkGroup 2:              -    D    E    F    G
  WorkGroup 1:              -    B    -    C    -
  WorkGroup 0:              A    -    -    -    -
  InterLevelGroup 1-2 le:   -   B+D   E    -    -
  InterLevelGroup 1-2 ri:   -    -    -   C+F   G
  InterLevelGroup 0-1:      A    B    -    C    -
  \endverbatim
  *
  * -# case c, 7 physical processes:
  \verbatim
  comp. proc. group rank:   0    1    2    3    4    5    6
  WorkGroup 2:              -    -    -    D    E    F    G
  WorkGroup 1:              -    B    C    -    -    -    -
  WorkGroup 0:              A    -    -    -    -    -    -
  InterLevelGroup 1-2 le:   -    B    -    D    E    -    -
  InterLevelGroup 1-2 ri:   -    -    C    -    -    F    G
  InterLevelGroup 0-1:      A    B    C    -    -    -    -
  \endverbatim
  *
  * \tparam space_dim_
  * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
  *
  * \tparam world_dim_
  * world dimension (determines the number of coordinates)
  *
  * \note Since a manager manages one (global) process group, it could also have been called 'ProcessGroupManager'.
  * But since 99 % of the code 'happens' within process groups (only application programmers dealing with multiphysics
  * or similar stuff have to organise different process groups), I omitted the leading 'ProcessGroup'.
  * Other names which came to my mind were 'ProgramFlowManager' (too long) and 'FlowManager' (misleading in a FE tool
  * dealing with CFD problems). (Hilmar)
  *
  * \author Hilmar Wobker
  * \author Dominik Goeddeke
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class Manager
  {

  private:

    /* *****************
    * member variables *
    *******************/
    /// process group the manager manages (initially created by the universe)
    ProcessGroup* _process_group_main;

    /// flag whether the main process group uses a dedicated load balancer process
    bool const _group_has_dedicated_load_bal;

    /**
    * \brief rank of the dedicated load balancer process (if there is one)
    *
    * If a dedicated load balancer process is used, then it has the largest rank in the main process group.
    */
    int _rank_dedicated_load_balancer;

    /// rank of the load balancer process (equal to rank of the coordinator)
    int _rank_load_balancer;

    /**
    * \brief process group containing only the worker processes of the main process group #_process_group
    *
    * If there is no dedicated load balancer, then the worker process group and the main process group are identical
    */
    ProcessGroup* _process_group_comp;

    ///
    ManagerCompCoord<space_dim_, world_dim_>* _manager_comp_coord;

    ///
    ManagerCompNonCoord<space_dim_, world_dim_>* _manager_comp_non_coord;

    ///
    ManagerComp<space_dim_, world_dim_>* _manager_comp;

    ///
    LoadBalancer<space_dim_, world_dim_>* _load_balancer;

    ///
    LoadBalancerDedicated<space_dim_, world_dim_>* _load_balancer_dedicated;


  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /// CTOR
    Manager(
      ProcessGroup* process_group_main,
      bool group_has_dedicated_load_bal)
      : _process_group_main(process_group_main),
        _group_has_dedicated_load_bal(group_has_dedicated_load_bal),
        _rank_dedicated_load_balancer(MPI_PROC_NULL),
        _rank_load_balancer(MPI_PROC_NULL),
        _process_group_comp(nullptr),
        _manager_comp_coord(nullptr),
        _manager_comp_non_coord(nullptr),
        _manager_comp(nullptr),
        _load_balancer(nullptr),
        _load_balancer_dedicated(nullptr)
    {
      CONTEXT("Manager::Manager()");

      // if the group contains a dedicated load balancer, it is the process with the highest rank in the process group
      if(_group_has_dedicated_load_bal)
      {
        _rank_dedicated_load_balancer = _process_group_main->num_processes()-1;

        // Create ProcessGroup object representing the group of all processes excluding the dedicated load bal. process.
        // Note that *all* processes of the parent MPI group have to call the MPI_Comm_create() routine (otherwise the
        // forking will deadlock), so let the ded. load balancer call it as well. Since the dedicated load bal. does not
        // belong to the group, the new communicator is MPI_COMM_NULL on the dedicated load bal. process.
        MPI_Group gr_without_dedicated_load_bal;
        MPI_Comm gr_comm;
        int mpi_error_code = MPI_Group_excl(_process_group_main->group(), 1, &_rank_dedicated_load_balancer,
                                            &gr_without_dedicated_load_bal);
        validate_error_code_mpi(mpi_error_code, "MPI_Group_excl");
        mpi_error_code = MPI_Comm_create(_process_group_main->comm(), gr_without_dedicated_load_bal, &gr_comm);
        validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");
        if(!is_dedicated_load_balancer())
        {
          _process_group_comp = new ProcessGroup(gr_comm);
        }
      }
      else
      {
        // otherwise, the compute group equals the main group
        _process_group_comp = _process_group_main;
      }

      // Set rank of the load balancer process to that of the coordinator. This makes sense since some data is shared
      // by load balancer and coordinator (as, for example, the base mesh). If the load balancer lives on a different
      // process, such data would have to be duplicated.
      _rank_load_balancer = _process_group_main->rank_coord();

      if(!is_dedicated_load_balancer())
      {
        if(is_coordinator_comp())
        {
          // create manager coordinator object in the coordinator process of the process group
          _manager_comp_coord = new ManagerCompCoord<space_dim_, world_dim_>(_process_group_comp);
          _manager_comp = _manager_comp_coord;
        }
        else
        {
          // create manager non-coordinator object in each non-coordinator process of the process group
          _manager_comp_non_coord = new ManagerCompNonCoord<space_dim_, world_dim_>(_process_group_comp);
          _manager_comp = _manager_comp_non_coord;
        }

        // load balancer is created on coordinator process
        if(is_load_balancer())
        {
          _load_balancer = new LoadBalancer<space_dim_, world_dim_>(_process_group_main, _group_has_dedicated_load_bal,
                                                                    _rank_dedicated_load_balancer);
          ASSERT(_manager_comp_coord != nullptr, "Coordinator and load balancer should live on the same process!");
          _manager_comp_coord->set_load_balancer(_load_balancer);
        }
      }
      else
      {
        // if there is a dedicated load balancer process create corresponding object on this process
        _load_balancer_dedicated = new LoadBalancerDedicated<space_dim_, world_dim_>(process_group_main);
      }
    }


    /// DTOR
    virtual ~Manager()
    {
      CONTEXT("Manager::~Manager()");
      if (_manager_comp_coord != nullptr)
      {
        delete _manager_comp_coord;
      }
      if (_manager_comp_non_coord != nullptr)
      {
        delete _manager_comp_non_coord;
      }
      if(is_load_balancer())
      {
        delete _load_balancer;
      }
      if(_group_has_dedicated_load_bal)
      {
        if(!is_dedicated_load_balancer())
        {
          delete _process_group_comp;
        }
        else
        {
          delete _load_balancer_dedicated;
        }
      }
    }


    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the flag whether the process group has a dedicated load balancer
    *
    * \return true when the process group has a dedicated load balancer, otherwise false
    */
    inline bool group_has_dedicated_load_bal() const
    {
      CONTEXT("Manager::group_has_dedicated_load_bal()");
      return _group_has_dedicated_load_bal;
    }


    /**
    * \brief getter for the coordinator manager of the compute process group living on this process
    *
    * \return pointer to coordinator manager of the compute process group
    */
    inline ManagerCompCoord<space_dim_, world_dim_>* manager_coord() const
    {
      CONTEXT("Manager::manager_coord()");
      return _manager_comp_coord;
    }

    /**
    * \brief getter for the non-coordinator manager of the compute process group living on this process
    *
    * \return pointer to non-coordinator manager of the compute process group
    */
    inline ManagerCompNonCoord<space_dim_, world_dim_>* manager_non_coord() const
    {
      CONTEXT("Manager::manager_non_coord()");
      return _manager_comp_non_coord;
    }


    /**
    * \brief getter for the main process group this process belongs to
    *
    * \return pointer to main process group
    */
    inline ProcessGroup* process_group_main() const
    {
      CONTEXT("Manager::process_group_main()");
      return _process_group_main;
    }


    /**
    * \brief getter for the compute process group this process belongs to
    *
    * \return pointer to compute process group
    */
    inline ProcessGroup* process_group_comp() const
    {
      CONTEXT("Manager::process_group_comp()");
      return _process_group_comp;
    }


    /**
    * \brief getter for the base mesh
    *
    * \return BaseMesh::BM<space_dim_, world_dim_> pointer #_base_mesh
    */
    inline BaseMesh::BM<space_dim_, world_dim_>* base_mesh() const
    {
      CONTEXT("Manager::base_mesh()");
      if(is_coordinator_comp())
      {
        return _manager_comp_coord->base_mesh();
      }
      else
      {
        return nullptr;
      }
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief checks whether this process is the dedicated load balancer of the process group
    *
    * \return true if this process is the dedicated load balancer of the process group, otherwise false
    */
    inline bool is_dedicated_load_balancer() const
    {
      CONTEXT("Manager::is_dedicated_load_balancer()");
      return _process_group_main->rank() == _rank_dedicated_load_balancer;
    }


    /**
    * \brief checks whether this process is the load balancer of the process group
    *
    * \return true if this process is the load balancer of the process group, otherwise false
    */
    inline bool is_load_balancer() const
    {
      CONTEXT("Manager::is_load_balancer()");
      return _process_group_main->rank() == _rank_load_balancer;
    }


    /**
    * \brief checks whether this process is the coordinator of the compute process group
    *
    * \return true if this process is the coordinator of the compute process group
    */
    inline bool is_coordinator_comp() const
    {
      CONTEXT("Manager::is_coordinator_comp()");
      if(!is_dedicated_load_balancer())
      {
        return _process_group_comp->is_coordinator();
      }
      else
      {
        return false;
      }
    }


    /// read in a mesh file and set up base mesh
    void read_mesh(std::string const & mesh_file)
    {
      CONTEXT("Manager::read_mesh()");

      // the mesh is read by the process group coordinator
      if(is_coordinator_comp())
      {
        ASSERT(_manager_comp_non_coord == nullptr, "_manager_comp_non_coord must be null!");
        ASSERT(_load_balancer_dedicated == nullptr, "_load_balancer_dedicated must be null!");
        ASSERT(_manager_comp_coord != nullptr, "_manager_comp_coord must not be null!");
        _manager_comp_coord->read_mesh(mesh_file);
        ASSERT(_load_balancer != nullptr, "_load_balancer must not be null!");
        _load_balancer->set_base_mesh(_manager_comp_coord->base_mesh());
      }
    }


    /// calls corresponding function in load balancer object
    void define_work_groups_1()
    {
      if(is_load_balancer())
      {
        ASSERT(_load_balancer != nullptr, "_load_balancer is null!");
        _load_balancer->define_work_groups_1();
      }
    }


    /// calls corresponding function in load balancer object
    void define_work_groups_2()
    {
      if(is_load_balancer())
      {
        ASSERT(_load_balancer != nullptr, "_load_balancer is null!");
        _load_balancer->define_work_groups_2();
      }
    }


    /// calls corresponding function in compute manager (coordinator and non-coordinators)
    void create_work_groups()
    {
      if(_manager_comp_coord != nullptr)
      {
        // Now let the manager create the work groups. Deallocation of arrays (except the graph array) and destruction
        // of objects is done within the manager class. This function also has to be called on the non-coordinator
        // processes of the process group.
        _manager_comp_coord->_create_work_groups();

        // now let the manager send local parts of the global graphs to the worker processes
        _manager_comp_coord->_transfer_graphs_to_workers(_load_balancer->graphs());
      }
      else if(_manager_comp_non_coord != nullptr)
      {
        // The non-coordinator processes also have to call the function for creating work groups. They receive the
        // necessary data from the coordinator.
        _manager_comp_non_coord->_create_work_groups();

        // let the non-coordinator processes receive the local parts of the global graphs
        _manager_comp_non_coord->_receive_and_set_graphs();
      }
    }


    /**
    * \brief dummy function for testing work group and interlevel group communication
    *
    * Calls corresponding function in compute manager (coordinator and non-coordinators).
    *
    * \return flag whether test was succesful (0: succesful, >0: failed)
    */
    unsigned int test_communication()
    {
      unsigned int test_result;
      CONTEXT("Manager::test_communication()");
      if(_manager_comp_coord != nullptr)
      {
        test_result = _manager_comp_coord->_test_communication();
      }
      else if(_manager_comp_non_coord != nullptr)
      {
        test_result = _manager_comp_non_coord->_test_communication();
      }
      else
      {
        // for the dedicated load balancing process, there is no senseful test (yet), to we simply return
        // the flag whether the group has a dedicated balancer (which should be true)
        if(_group_has_dedicated_load_bal)
        {
          test_result = 0;
        }
        else
        {
          test_result = 1;
        }
      }

      // Determine the maximum over the test results of all processes. Rationale: When one process failed, the whole
      // test should be considered as failed. So, let all processes return the same values.
      unsigned int test_result_max;
      MPI_Allreduce(&test_result, &test_result_max, 1, MPI_UNSIGNED, MPI_MAX, _process_group_main->comm());
      return test_result_max;
    } // test_communication
  };
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_MANAGER_HPP

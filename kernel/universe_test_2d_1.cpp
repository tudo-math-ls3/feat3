// Hilmar's todo list before push to master:
// code cleanup:
// - add/remove 'const' and other modifiers where necessary
// - completion of doxygen comments (\param, array dimensions, \return, \throws etc.)

// functionality:
// - tests splitten und an Dirks Test-System anpassen
// - weitere tests schreiben
// - check: missing output when using MPICH

// Done:
// - template all classes "above" BaseMesh
// - FileParser --> templates!
// - BaseMesh1D, BaseMesh2D, BaseMesh3D --> templates!
// - subdiv_data an subdivide uebergeben
// - move mpi_init() from Universe::create() to the very beginning of the program
// - Universe CTOR und DTOR private machen
// - proper cleanup of resources (deletes, DTORs)
// - check all the MPI-related valgrind warnings (--> Dom)
// - replace std::cout/std::err/exit(1) by Logger::... (or remove, or change into mpi_aborts/exceptions)
// - use ErrorHandler in base mesh code
// - CONTEXT("...") einfuegen
// - missing DTORs, virtual DTORs?
// - test_1d.cpp etc. ersetzen

// includes, system
#include <iostream>
#include <cstdlib> // for exit()

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/comm.hpp>
#include <kernel/graph.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>

using namespace FEAST;

// set space and world dimension manually to 2
#define SDIM 2
#define WDIM 2


/**
* \brief dummy function in preparation of a function for defining work groups (which will reside in the load bal. class)
*
* This dummy function creates two work groups: one consisting of two workers responsible for the coarse grid
* problem and one consisting of all the other workers responsible for the fine grid problem. Currently, everything
* is hard-coded. Later, the user must be able to control the creation of work groups and even later the load
* balancer has to apply clever strategies to create these work groups automatically so that the user doesn't have
* to do anything.
*
* To optimise the communication between the coordinator of the main process group and the work groups, we add
* this coordinator to a work group if it is not a compute process of this work group anyway. Hence, for each work
* group, there are three different possibilities:
* 1) There is a dedicated load balancer process, which is automatically the coordinator of the main process group
*    and belongs to no work group:
*    --> work group adds the coordinator as extra process
* 2) There is no dedicated load balancer process, and the coordinator of the main process group ...
*   a) ... is not part of the work group:
*     --> work group adds the coordinator as extra process
*   b) ... is part of the work group:
*     --> work group does not have to add the coordinator
* Thus the 1-to-n or n-to-1 communication between coordinator and n work group processes can be performed via
* MPI_Scatter() and MPI_Gather() (which always have to be called by all members of an MPI process group). This
* is more efficient then using n calls of MPI_send() / MPI_recv() via the communicator of the main process group.
*
* The following test code is semi hard-wired in the sense that we schedule one BMC to each processor.
* Furthermore, two extra processes are required for testing the creation of a second work group and the use
* of a dedicated load balancing process.
* Later, this has to be all done in some auto-magically way.
*
* Two different tests can be performed (n is the number of base mesh cells):
* 1) with dedicated load balancer process
*    - work group for coarse grid: 2 processes: {0, 1}
*    - work group for fine grid: n processes: {1, ..., n}
*    - i.e. process 1 is in both work groups
*    - dedicated load balancer and coordinator process: n+1
* 2) without dedicated load balancer process
*    - work group for coarse grid: 2 processes: {0, 1}
*    - work group for fine grid: n processes: {2, ..., n+1}
*    - i.e. the two work groups are disjunct
*    - coordinator process: n+1
* Both tests need n+2 processes in total. To choose the test, change the first entry of the boolean array
* includes_dedicated_load_bal[] in the main() method.
*
* \param[in] load_balancer
* pointer to the load balancer object to get some further information
*
* \param[out] num_subgroups
* number of subgroups
*
* \param[out] num_proc_in_subgroup
* array of numbers of processes per subgroup
*
* \param[out] group_contains_extra_coord
* Array indicating whether the subgroups contain an extra process for the coordinator (which will then
* not be a compute process in the corresponding work group).
* Usually a boolean would do the trick, but we cannot send boolean arrays via MPI. (C++ bindings are deprecated,
* hence we cannot use MPI::BOOL. MPI_LOGICAL cannot be used, since this is not equivalent to C++'s bool.)
* So, we use unsigned char here and treat is as if it were boolean.
*
* \param[out] subgroup_ranks
* array for defining the rank partitioning
*
* \param[out] graphs
* array of Graph pointers representing the connectivity of work group processes
*/
void define_work_groups(
  LoadBalancer<SDIM, WDIM>* load_balancer,
  unsigned int& num_subgroups,
  unsigned int*& num_proc_in_subgroup,
  unsigned char*& group_contains_extra_coord,
  int**& subgroup_ranks,
  Graph**& graphs)
{
  CONTEXT("universe_test_2d_1: define_work_groups()");
  // set number of subgroups manually to 2
  num_subgroups = 2;
  // allocate arrays
  num_proc_in_subgroup = new unsigned int[num_subgroups];
  group_contains_extra_coord = new unsigned char[num_subgroups];
  subgroup_ranks = new int*[num_subgroups];
  graphs = new Graph*[num_subgroups];

  // shortcut to the number of processes in the load balancer's process group
  unsigned int num_processes = load_balancer->process_group()->num_processes();

  // shortcut to the number of cells in the base mesh
  unsigned int num_cells = load_balancer->base_mesh()->num_cells();

  // debug output
  Logger::log_master("num_cells: " + StringUtils::stringify(num_cells) + "\n");
  // assert that the number of processes is n+2
  assert(num_processes == num_cells + 2);

  // set up the two test cases

  if(load_balancer->group_has_dedicated_load_bal())
  {
    // test case 1 (n is the number of base mesh cells)
    // with dedicated load balancer process
    //  - work group for coarse grid: 2 processes: {0, 1}
    //  - work group for fine grid: n processes: {1, ..., n}
    //  - i.e. process 1 is in both work groups
    //  - dedicated load balancer and coordinator process: n+1

    // since there is a dedicated load balancer process, this has to be added to both work groups as extra
    // coordinator process
    group_contains_extra_coord[0] = true;
    group_contains_extra_coord[1] = true;

    // set number of processes per group
    num_proc_in_subgroup[0] = 2 + 1;
    num_proc_in_subgroup[1] = num_cells + 1;

    // partition the process group ranks into work groups
    // coarse grid work group
    subgroup_ranks[0] = new int[num_proc_in_subgroup[0]];
    subgroup_ranks[0][0] = 0;
    subgroup_ranks[0][1] = 1;
    subgroup_ranks[0][2] = num_cells + 1;

    // fine grid work group
    subgroup_ranks[1] = new int[num_proc_in_subgroup[1]];
    // set entries to {1, ..., n+1}
    for(unsigned int i(0) ; i < num_proc_in_subgroup[1] ; ++i)
    {
      subgroup_ranks[1][i] = i+1;
    }
  }
  else
  {
    // test case 2 (n is the number of base mesh cells)
    // without dedicated load balancer process
    //  - work group for coarse grid: 2 processes: {0, 1}
    //  - work group for fine grid: n processes: {2, ..., n+1}
    //  - i.e. the two work groups are disjunct
    //  - coordinator process: n+1

    // the coordinator is at the same time a compute process of the second work group, so only the first work group
    // has to add an extra process
    group_contains_extra_coord[0] = true;
    group_contains_extra_coord[1] = false;

    // set number of processes per group
    num_proc_in_subgroup[0] = 2 + 1;
    num_proc_in_subgroup[1] = num_cells;

    // partition the process group ranks into work groups
    subgroup_ranks[0] = new int[num_proc_in_subgroup[0]];
    subgroup_ranks[0][0] = 0;
    subgroup_ranks[0][1] = 1;
    subgroup_ranks[0][2] = num_cells + 1;
    subgroup_ranks[1] = new int[num_proc_in_subgroup[1]];
    // set entries to {2, ..., n+1}
    for(unsigned int i(0) ; i < num_proc_in_subgroup[1] ; ++i)
    {
      subgroup_ranks[1][i] = i+2;
    }
  }

  /* *********************************************************
  * create graph structures corresponding to the work groups *
  ***********************************************************/
  // build an artificial graph mimicing the distribution of the 16 base mesh cells to two processors
  // (e.g. BMCs 0-7 on proc 1 and BMCs 8-15 on proc 2) which start an imagined coarse grid solver; this graph will
  // be used for the coarse grid work group
  unsigned int* index = new unsigned int[3];
  unsigned int* neighbours = new unsigned int[2];
  index[0] = 0;
  index[1] = 1;
  index[2] = 2;
  neighbours[0] = 1;
  neighbours[1] = 0;
  // Artificially create a graph object here. Usually, this comes from somewhere else.
  graphs[0] = new Graph(2, index, neighbours);
  // arrays index and neighbours are copied within the Graph CTOR, hence they can be deallocated here
  delete [] index;
  delete [] neighbours;

  // get connectivity graph of the base mesh; this one will be used for the fine grid work group
  graphs[1] = load_balancer->base_mesh()->graph();

// COMMENT_HILMAR:
// We assume here that each process receives exactly one BMC and that the index of the cell in the graph structure
// equals the local rank within the work group. Later, there will be the matrix patch layer and the process patch
// layer, which both have their own connectivity structure. Here, we actually need the connectivity graph of the
// process patch layer.
}


/**
* \brief main routine - test driver for universe
*
* Call the routine via "mpirun -np n+5 <binary_name> <relative path to base mesh> n" where n is the number of base
* mesh cells in the mesh. The semi-hard-wired example then creates two 'main process groups' (for two completely
* independent problems). Each process group has its own load balancer. Only the first process group (consisting of
* n+2 processes does some useful things, the second one (consisting of 2 processes) is idle.
*
* \todo set up separate tests for universe, work group creation, logging/pretty printer, cell subdivision
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
int main(int argc, char* argv[])
{
  CONTEXT("universe_test_2d_1: main()");
  if (argc < 3)
  {
    std::cerr << "Call the program with \"mpirun -np n+5 " << argv[0] << " <relative_path_to_mesh_file> n\", "
              << "where n is the number of base mesh cells." << std::endl;
    exit(1);
  }

  // init MPI
  MPIUtils::init_MPI(argc, argv);

  // COMMENT_HILMAR: the following information will later be read from some parameter file

  // number of process groups (if not provided, then 1)
  unsigned int num_process_groups(2);

  // COMMENT_HILMAR:
  // For the semi-hard-wired example using a mesh with n base mesh cells we need n+5 processes in total
  // (n+2 for the first process group, 2 for the second process group and 1 for the master process).
  // (See description of the example in routine define_work_groups(...).)

  // COMMENT_HILMAR:
  // As an intermediate hack, the number of base mesh cells has to be provided as first argument to the program call.

  unsigned int num_cells = atoi(argv[2]);
  assert(num_cells > 0);

  // The number of processes for the first process group must equal num_cells + 2.
  unsigned int num_processes_in_first_group(num_cells + 2);

  // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
  unsigned int* num_processes_in_group = new unsigned int[num_process_groups];
  num_processes_in_group[0] = num_processes_in_first_group;
  num_processes_in_group[1] = 2;

  // array of flags whether a dedicated load balancer process is needed in process groups
  // (must be provided when num_process_groups > 1)
  // Set the first entry either to false or to true to test two different configurations. (You don't have to change
  // the number of processes for that.) (See description of the example in routine LoadBalancer::create_work_groups().)
  bool* includes_dedicated_load_bal = new bool[num_process_groups];
  includes_dedicated_load_bal[0] = false;
  includes_dedicated_load_bal[1] = false;

  // set shortcut to the one and only instance of Universe (since this is the first call of
  // Universe<SDIM, WDIM>::instance(), it also calls the constructor of the Universe singleton class)
  Universe<SDIM, WDIM>* universe = Universe<SDIM, WDIM>::instance();

  try
  {
    universe->create(num_process_groups, num_processes_in_group, includes_dedicated_load_bal);
  }
  catch (Exception& e)
  {
    // abort the program
    ErrorHandler::exception_occured(e);
  }

  // Get process objects. Note that on each process only one of the following two exists (the other one is the
  // null pointer).
  LoadBalancer<SDIM, WDIM>* load_balancer = universe->load_balancer();
  Master* master = universe->master();

  int rank_world = Process::rank;

  if(load_balancer != nullptr)
  {
    ProcessGroup* process_group = load_balancer->process_group();
    unsigned int group_id = process_group->group_id();
//    // debug output
//    int rank_process_group = process_group->rank();
//    std::string s("Process " + StringUtils::stringify(rank_world) + " is the ");
//    if(load_balancer->is_dedicated_load_bal())
//    {
//      s += "DEDICATED ";
//    }
//    s += "load balancer with local rank " + StringUtils::stringify(rank_process_group);
//    s += " in group " + StringUtils::stringify(group_id) + ".\n";
//    process_group->log_indiv_master(s);

    // perform actions depending on the group id
    if(group_id == 0)
    {
      // get name of the mesh file from the command line
      std::string mesh_file(argv[1]);
      // let the load balancer read the mesh file and create a base mesh (this is only done by the coordinator process)
      load_balancer->read_mesh(mesh_file);

      // necessary data to be set for creating work groups (see function define_work_groups(...) for details)
      unsigned int num_subgroups;
      unsigned int* num_proc_in_subgroup(nullptr);
      unsigned char* group_contains_extra_coord(nullptr);
      int** subgroup_ranks(nullptr);
      Graph** graphs;

      // define work groups
      if(process_group->is_coordinator())
      {
        // The coordinator is the only one knowing the base mesh, so only the coordinator decides over the number of
        // work groups and the process distribution to them. The passed arrays are allocated within this routine.
        define_work_groups(load_balancer, num_subgroups, num_proc_in_subgroup, group_contains_extra_coord,
                           subgroup_ranks, graphs);
      }
      // Now let the load balancer create the work groups (this function is called on all processes of the process
      // group). Deallocation of arrays (except the graph array) and destruction of objects is done within the load
      // balancer class.
      load_balancer->create_work_groups(num_subgroups, num_proc_in_subgroup, group_contains_extra_coord,
                                        subgroup_ranks, graphs);

      // let some process test the PrettyPrinter, the vector version of the function log_master_array() and the
      // standard file logging functions
      if (Process::rank % 7 == 0)
      {
        // test the PrettyPrinter
        std::string prefix(std::string("Proc" + StringUtils::stringify(Process::rank)));
        PrettyPrinter pp(40, '#', prefix + " ");
        pp.add_line_sep();
        pp.add_line_centered("Testing pp and logging!");
        pp.add_line_sep();
        pp.add_line("left bla blub");
        pp.add_line("too long too long too long too long too long too long");
        pp.add_line_no_right_delim("also too long too long too long too long too long too long");
        pp.add_line("left bla blub");
        pp.add_line_sep();
        // print it like this...
        Logger::log(pp.block());
        // ... or like this...
        pp.print(Logger::file);
        // send it to the master for file output
        Logger::log_master(pp.block(), Logger::FILE);

        // test vector version of log_master_array()
        std::vector<std::string> messages(3);
        messages[0] = prefix + ": Testing this...\n";
        messages[1] = prefix + ": ...vector logging...\n";
        messages[2] = prefix + ": ...feature!\n";
        Logger::log_master_array(messages, Logger::FILE);
        Logger::log(messages);

        // test standard log feature
        Logger::log("BRAL\n");
      }

      if(process_group->is_coordinator())
      {
        // manually destroy here the artificially created graph object (the other one is destroyed by the base mesh)
        delete graphs[0];
        // delete the graphs array
        delete [] graphs;
      }

      // Everything done, call universe destruction routine.
      universe->destroy();
    }
    else
    {
      // the second process group does something else, programmed by the application outside the kernel...
      // ...
      // Everything done, call universe destruction routine.
      universe->destroy();
   }
  }
  else if(master != nullptr)
  {
    // This branch is entered when the infinite service loop of the master has been finished.
    // This, however, usually happens only at program end.
    // Everything done, call universe destruction routine.
    universe->destroy();
  }
  else
  {
    MPIUtils::abort("Process with rank " + StringUtils::stringify(rank_world)
                    + " has no particular role, this should not happen.");
  }
}

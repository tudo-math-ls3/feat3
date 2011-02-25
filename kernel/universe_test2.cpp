// includes, system
#include <iostream>
#include <cstdlib> // for exit()

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/comm.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>

using namespace FEAST;

#define WDIM 2
#define SDIM 2

/**
* \brief dummy function in preparation of a function for defining work groups (which will reside in the load bal. class)
*
* This dummy function creates one work group consisting of all available processes.
* This test code is semi hard-wired in the sense that we schedule one BMC to each processor.
* We need n processes for the following (where n is the number of base mesh cells).
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
  // set number of subgroups manually to 1
  num_subgroups = 1;
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
  // assert that the number of processes is n
  assert(num_processes == num_cells);
  // assert that no dedicated load balancer is used
  assert(!load_balancer->group_has_dedicated_load_bal());

  // set up the test case
  // test case 1 (n is the number of base mesh cells)
  // with dedicated load balancer process
  //  - work group for coarse grid: 2 processes: {0, 1}
  //  - work group for fine grid: n processes: {1, ..., n}
  //  - i.e. process 1 is in both work groups
  //  - dedicated load balancer and coordinator process: n+1

  // since there is a dedicated load balancer process, this has to be added to both work groups as extra
  // coordinator process
  group_contains_extra_coord[0] = false;

  // set number of processes per group
  num_proc_in_subgroup[0] = num_cells;

  // fine grid work group
  subgroup_ranks[0] = new int[num_proc_in_subgroup[0]];
  // set entries to {0, ..., n-1}
  for(unsigned int i(0) ; i < num_proc_in_subgroup[0] ; ++i)
  {
    subgroup_ranks[0][i] = i;
  }

  /* *********************************************************
  * create graph structures corresponding to the work groups *
  ***********************************************************/
  // get connectivity graph of the base mesh; this one will be used for the one and only work group
  graphs[0] = load_balancer->base_mesh()->graph();

// COMMENT_HILMAR:
// We assume here that each process receives exactly one BMC and that the index of the cell in the graph structure
// equals the local rank within the work group. Later, there will be the matrix patch layer and the process patch
// layer, which both have their own connectivity structure. Here, we actually need the connectivity graph of the
// process patch layer.
}



/**
* \brief main routine - test driver for universe
*
* Call the routine via "mpirun -np n+3 <binary_name> <relative path to base mesh> n" where n is the number of base
* mesh cells in the mesh.
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
int main(int argc, char* argv[])
{

  if (argc < 3)
  {
    std::cerr << "Call the program with \"mpirun -np n+1 " << argv[0] << " <relative_path_to_mesh_file> n\", "
              << "where n is the number of base mesh cells." << std::endl;
    exit(1);
  }

  // init MPI
  MPIUtils::init_MPI(argc, argv);

  // COMMENT_HILMAR:
  // For the semi-hard-wired example using a mesh with n base mesh cells we need n+3 processes in total
  // (n+2 for the one for the one and only process group and 1 for the master process).
  // (See description of the example in routine LoadBalancer::create_subgroups().)

  // COMMENT_HILMAR:
  // As an intermediate hack, the number of base mesh cells has to be provided as first argument to the program call.

  unsigned int num_cells;
  num_cells = atoi(argv[2]);
  assert(num_cells > 0);

  // set shortcut to the one and only instance of Universe (since this is the first call of
  // Universe<SDIM, WDIM>::instance(), it also calls the constructor of the Universe singleton class)
  Universe<SDIM, WDIM>* universe = Universe<SDIM, WDIM>::instance();

  try
  {
    universe->create();
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
//    // debug output
//    int rank_process_group = process_group->rank();
//    std::string s("Process " + StringUtils::stringify(rank_world) + " is the ");
//    if(load_balancer->is_dedicated_load_bal())
//    {
//      s += "DEDICATED ";
//    }
//    s += "load balancer with local rank " + StringUtils::stringify(rank_process_group);
//    s += " in group " + StringUtils::stringify(group_id) + ".";
//    process_group->log_indiv_master(s);

    // get name of the mesh file from the command line
    std::string mesh_file(argv[1]);
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
     // delete the graphs array
      delete [] graphs;
    }

    // Everything done, call universe destruction routine.
    universe->destroy();
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

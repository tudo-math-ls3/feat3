//// includes, system
//#include <iostream>
//#include <cstdlib> // for exit()
//
// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/comm.hpp>
#include <kernel/graph.hpp>
#include <kernel/process.hpp>
#include <kernel/manager.hpp>
#include <kernel/universe.hpp>

using namespace FEAST;

#define WDIM 2
#define SDIM 2


/**
* \brief main routine - test driver for universe
*
* Call the routine via "mpirun -np n+5 <binary_name> <absolute path to base mesh> n" where n is the number of base
* mesh cells in the mesh. The semi-hard-wired example then creates two 'main process groups' (for two completely
* independent problems). Each process group has its own manager. Only the first process group (consisting of
* n+2 processes does some useful things, the second one (consisting of 2 processes) is idle.
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
int main(int argc, char* argv[])
{
  CONTEXT("pseudo_app: main()");
  if (argc < 3)
  {
    std::cerr << "Call the program with \"mpirun -np n+5 " << argv[0] << " <absolute_path_to_mesh_file> n\", "
              << "where n is the number of base mesh cells." << std::endl;
    exit(1);
  }

  // init MPI
  init_mpi(argc, argv);

  // COMMENT_HILMAR: the following information will later be read from some parameter file

  // number of process groups (if not provided, then 1)
  unsigned int num_process_groups(2);

  // COMMENT_HILMAR:
  // For the semi-hard-wired example using a mesh with n base mesh cells we need n+5 processes in total
  // (n+2 for the first process group, 2 for the second process group and 1 for the master process).
  // (See description of the example in routine define_work_groups().)

  // COMMENT_HILMAR:
  // As an intermediate hack, the number of base mesh cells has to be provided as first argument to the program call.

  unsigned int num_cells;
  num_cells = atoi(argv[2]);
  ASSERT(num_cells > 0, "Number of cells must not be zero.");

  // The number of processes for the first process group must equal num_cells + 2.
  unsigned int num_processes_in_first_group(num_cells + 2);

  // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
  unsigned int* num_processes_in_group = new unsigned int[num_process_groups];
  num_processes_in_group[0] = num_processes_in_first_group;
  num_processes_in_group[1] = 2;

  // array of flags whether a dedicated load balancer process is needed in process groups
  // (must be provided when num_process_groups > 1)
  // Set the first entry either to false or to true to test two different configurations. (You don't have to change
  // the number of processes for that.) (See description of the example in routine define_work_groups().)
  bool* includes_dedicated_load_bal = new bool[num_process_groups];
  includes_dedicated_load_bal[0] = true;
  includes_dedicated_load_bal[1] = false;

  // set shortcut to the one and only instance of Universe (since this is the first call of
  // Universe<SDIM, WDIM>::instance(), it also calls the constructor of the Universe singleton class)
  Universe<SDIM, WDIM>* universe = Universe<SDIM, WDIM>::instance();

  try
  {
    universe->create(num_process_groups, num_processes_in_group, includes_dedicated_load_bal, "pseudo_app");
  }
  catch (Exception& e)
  {
    // abort the program
    ErrorHandler::exception_occured(e);
  }

  // Get process objects. Note that on each process only one of the following two exists (the other one is
  // automatically the null pointer).
  Manager<SDIM, WDIM>* manager = universe->manager();
  Master* master = universe->master();

  if(manager != nullptr)
  {
    // get pointer to the process group this process belongs to
    ProcessGroup* process_group = manager->process_group_main();
    // get group id
    unsigned int group_id = process_group->group_id();

    // debug output
    int rank_process_group = process_group->rank();
    std::string s("Process " + stringify(Process::rank) + " is the manager with local rank "
                  + stringify(rank_process_group) + " in group " + stringify(group_id) + ".\n");
    Logger::log(s);

    // perform actions depending on the group id
    if(group_id == 0)
    {
      // get name of the mesh file from the command line
      std::string mesh_file(argv[1]);
      // let the manager read the mesh file and create a base mesh (this is only done by the coordinator process)
      manager->read_mesh(mesh_file);

      // let the load balancer define the work groups (currently hard-coded)
      manager->define_work_groups_2();

      // Now let the manager create the work groups. Deallocation of arrays (except the graph array) and destruction
      // of objects is done within the manager class. This function also has to be called on the non-coordinator
      // processes of the process group.
      manager->create_work_groups();

      // let some process test the PrettyPrinter, the vector version of the function log_master_array() and the
      // standard file logging functions
      if (Process::rank % 7 == 0)
      {
        // test the PrettyPrinter
        std::string prefix(std::string("Proc" + stringify(Process::rank)));
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

      // Everything done, call universe destruction routine.
      universe->destroy();
    }
    else // group_id != 0
    {
      // the second process group does something else, programmed by the application outside the kernel...
      Logger::log("Proc " + stringify(Process::rank) + " in group with ID " + stringify(group_id)
                  + " does something useful...\n");
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
    // This branch must not be entered. Throw InternalError which is caught by outer test system.
    throw InternalError("Process with rank " + stringify(Process::rank)
                        + " has no particular role, this should not happen.");
  }
  // let all processes call function ot finalise MPI
  finalise_mpi();
}

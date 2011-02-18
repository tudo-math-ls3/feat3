// Hilmar's todo list before push to master:
// code cleanup:
// - use ExceptionHandler in base mesh code
// - replace std::cout/std::err/exit(1) by Logger::... (or remove, or change into
//   mpi_aborts/exceptions)
// - singleton-Universe via Dirk's util class
// - move mpi_init() from Universe::create() to the very beginning of the program
// - completion of doxygen comments (\param, array dimensions, etc.)
// - add/remove 'const' and other modifiers where necessary
// - proper cleanup of resources (deletes, DTORs)

// Done:
// - FileParser --> templates!
// - BaseMesh1D, BaseMesh2D, BaseMesh3D --> templates!

// Test driver for universe

// includes, system
#include <iostream>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/comm.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>

using namespace FEAST;

/**
* \brief main routine
*
* Call the routine via "mpirun -np n+5 <binary_name> <relative path to base mesh> n" where n is the number of base
* mesh cells in the mesh. The semi-hard-wired example then creates two 'main process groups' (for two completely
* independent problems). Each process group has its own load balancer. Only the first process group (consisting of
* n+2 processes does some senseful things, the second one (consisting of 2 processes) is idle.
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
int main(int argc, char* argv[])
{

  if (argc < 3)
  {
    std::cerr << "Call the program with \"mpirun -np n+5 " << argv[0] << " <relative_path_to_mesh_file> n\", "
              << "where n is the number of base mesh cells." << std::endl;
    return 1;
  }

  // create universe with one process group
//  Universe* universe = Universe::create(argc, argv);

  // COMMENT_HILMAR: the following information will later be read from some parameter file

  // number of process groups (if not provided, then 1)
  unsigned int num_process_groups(2);

  // COMMENT_HILMAR:
  // For the semi-hard-wired example using a mesh with n base mesh cells we need n+5 processes in total
  // (n+2 for the first process group, 2 for the second process group and 1 for the master process).
  // (See description of the example in routine LoadBalancer::create_subgroups().)

  // COMMENT_HILMAR:
  // As an intermediate hack, the number of base mesh cells has to be provided as first argument to the program call.
  // If no argument is provided, then it is assumed that the mesh contains 16 base mesh cells.

  unsigned int num_cells;
  num_cells = atoi(argv[2]);
  assert(num_cells > 0);
//  std::cout << "Number of base mesh cells read from command line: " << num_cells << std::endl;

  // The number of processes for the first process group must equal num_cells + 2.
  unsigned int num_processes_in_first_group(num_cells + 2);

  // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
  unsigned int num_processes_in_group[] = {num_processes_in_first_group, 2};

  // array of flags whether a dedicated load balancer process is needed in process groups
  // (must be provided when num_process_groups > 1)
  // Set the first entry either to false or to true to test two different configurations. (You don't have to change
  // the number of processes for that.) (See description of the example in routine LoadBalancer::create_subgroups().)
  bool includes_dedicated_load_bal[] = {false, false};

// COMMENT_HILMAR: Is it problematic to define and fill the variables above already *before* MPI_Init() is called?
// (Until now, no problems occured!) If it turns out to be problematic then the Universe::create(...) routine has to be
// split in such a way that it first calls MPI_Init() *only*, and then (after setting the variables above) doing
// all the remaining stuff in a second step.

  // create universe with several process groups
  Universe* universe;
  try
  {
    universe = Universe::create(argc, argv, num_process_groups, num_processes_in_group, includes_dedicated_load_bal);
  }
  catch (Exception& e)
  {
    // Assume that all critical errors are already caught within the Universe class.
    // Since it is not clear whether the communication system has already been set up, do not forward the error message
    // to the master.
    ErrorHandler::exception_occured(e, ErrorHandler::NON_CRITICAL);
  }

  // create a second universe to test the error handler
  Universe* universe2;
  try
  {
    universe2 = Universe::create(argc, argv, num_process_groups, num_processes_in_group, includes_dedicated_load_bal);
  }
  catch (Exception& e)
  {
    // Assume that all critical errors are already caught within the Universe class.
    ErrorHandler::exception_occured(e, ErrorHandler::NON_CRITICAL);
  }

  // Get process objects. Note that on each process only one of the following two exists (the other one is the
  // null pointer).
  LoadBalancer* load_balancer = universe->load_balancer();
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
//    s += " in group " + StringUtils::stringify(group_id) + ".";
//    std::cout << s <<  std::endl;
    if(group_id == 0)
    {
      // get name of the mesh file from the command line
      std::string mesh_file(argv[1]);
      load_balancer->read_mesh(mesh_file);

      load_balancer->create_subgroups();

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
      // everything done, destroy the universe
      Universe::destroy();
    }
    else
    {
      // the second process group does something else, programmed by the application outside the kernel...
      // ...
      // everything done, destroy the universe
      Universe::destroy();
    }
  }
  else if(master != nullptr)
  {
    // This branch is entered when the infinite service loop of the master has been finished.
    // This, however, usually happens only at program end. Hence, destroy the universe.
    Universe::destroy();
  }
  else
  {
    MPIUtils::abort("Process with rank " + StringUtils::stringify(rank_world)
                    + " has no particular role, this should not happen.");
  }
}

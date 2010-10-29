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
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
int main(int argc, char* argv[])
{

  // create universe with one process group
//  Universe* universe = Universe::create(argc, argv);

  // the following information will be read from some dat file

  // number of process groups (if not provided, then 1)
  int num_process_groups(2);

  // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
  // For the hard-coded example mesh we need 21 processes. (See description of the example in routine
  // LoadBalancer::create_subgroups().)
  int num_processes_in_group[] = {18, 2};

  // array of flags whether a dedicated load balancer process is needed in process groups
  // (must be provided when num_process_groups > 1)
  // Set the first entry either to false or to true to test two different configurations. (You don't have to change
  // the number of processes for that.) (See description of the example in routine LoadBalancer::create_subgroups().)
  bool includes_dedicated_load_bal[] = {true, false};

  // create universe with several process groups
  Universe* universe = Universe::create(argc, argv, num_process_groups, num_processes_in_group,
                                        includes_dedicated_load_bal);

  // Get process objects. Note that on each process only one of the following two exists (the other one is the
  // null pointer).
  LoadBalancer* load_balancer = universe->load_balancer();
  Master* master = universe->master();

  int rank_world = Process::rank;

  if(load_balancer != nullptr)
  {
    ProcessGroup* process_group = load_balancer->process_group();
    int group_id = process_group->group_id();
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
      load_balancer->read_mesh();
      load_balancer->create_subgroups();

      // let some process test the PrettyPrinter, the vector version of the function log_master_array() and the
      // standard file logging functions
      if (Process::rank % 7 == 0)
      {
        // test the PrettyPrinter
        std::string prefix(std::string("Proc"+StringUtils::stringify(Process::rank)));
        PrettyPrinter pp(40, '#', prefix);
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

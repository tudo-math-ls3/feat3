// Test driver for universe

// includes, system
#include <iostream>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/universe.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>

// main routine
// @author Hilmar Wobker
// @author Dominik Goeddeke
int main(int argc, char* argv[])
{

  // create universe with one process group
//  Universe* universe = Universe::create(argc, argv);
//  Universe::destroy();

  // the following information will be read from some dat file

  // number of process groups (if not provided, then 1)
  int num_process_groups(2);
  // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
  int num_processes_in_group[] = {2, 3};
  // array of flags whether a dedicated load balancer process is needed in process groups
  // (must be provided when num_process_groups > 1)
  bool includes_dedicated_load_bal[] = {false, true};


  // create universe with several process groups
  Universe* universe = Universe::create(argc, argv, num_process_groups, num_processes_in_group,
                                        includes_dedicated_load_bal);

  // Get process objects. Note that on each process only one of the following three exists (the other two are
  // null pointers.).
  LoadBalancer* load_balancer = universe->load_balancer();
  Master* master = universe->master();

  if (load_balancer != nullptr)
  {
    ProcessGroup* process_group = load_balancer->process_group();
    int rank_world = load_balancer->rank_world();
    int group_id = process_group->group_id();
    int rank_process_group = process_group->my_rank();
    std::string s("Process " + StringUtils::stringify(rank_world) + " is the ");
    if (load_balancer->dedicated_load_bal_process())
    {
      s += "DEDICATED ";
    }
    s += "load balancer with local rank ";
    s += StringUtils::stringify(rank_process_group);
    s += " in group ";
    s += StringUtils::stringify(group_id) + ".";
    std::cout << s <<  std::endl;
    if (group_id == 0)
    {
      load_balancer->read_mesh();
// BRAL: temporarily deactivated
//      load_balancer->create_work_groups();
    }
    else
    {
      // the second process group does something else, programmed by the application outside the kernel...
    }
  }
  else if (master != nullptr)
  {
    // not sure yet if it makes sense to let the user control the master process
    int rank_world = master->rank_world();
    std::cout << "Process " << rank_world << " is the MASTER OF THE UNIVERSE!" << std::endl;
  }
  else
  {
    // rank of this process
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPIUtils::abort("Process with rank " + StringUtils::stringify(my_rank)
                    + " has no particular role, this should not happen.");
  }

  Universe::destroy();
}

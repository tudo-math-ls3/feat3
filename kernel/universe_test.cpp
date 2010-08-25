// Test driver for universe

// includes, system
#include <iostream>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/universe.hpp>

// main routine
// @author Hilmar Wobker
// @author Dominik Goeddeke
int main(int argc, char* argv[])
{
  // the following information will be read from some dat file

  // number of process groups (if not provided, then 1)
  const int num_process_groups;
  num_process_groups = 2;
  // number of processes in process groups (must be provided when num_process_groups > 1)
  int num_processes_in_group[num_process_groups];

  num_processes_in_group[0] = 2;
  num_processes_in_group[1] = 3;

  // create universe

  // simple constructor for creating one process group
  // Universe universe(argc, argv);

  // constructor for creating more than one process group
  // @Hilmar: So lebt das ganze Universum auf dem Stack und nicht auf dem Heap. Absicht?
  Universe universe(argc, argv, num_process_groups, num_processes_in_group);

  // Get process objects. Note that on each process only one of the following three exists (the other two are
  // null pointers.).
  LoadBalancer* load_balancer = universe.get_load_balancer();
  GroupProcess* group_process = universe.get_group_process();
  Master* master = universe.get_master();

  // rank of this process
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (load_balancer != nullptr)
  {
    int rw = load_balancer->get_rank_world();
    int gid = load_balancer->get_group_id();
    int rl = load_balancer->get_rank_local();
    std::cout << "process with world rank " << rw << " is the load balancer with id "
              << gid << " and local rank " << rl << std::endl;
    if (load_balancer->get_group_id() == 0)
    {
      load_balancer->read_mesh();
      load_balancer->create_work_groups();
    }
    else
    {
      // the second process group does something else, programmed by the application outside the kernel...
    }
  }
  else if (group_process != nullptr) // BUG: group_processes and masters are already in their infinite loop,
                                     // and thus this output never shows
  {
    // not sure yet if it makes sense to let the user control the group processes
    int rw = group_process->get_rank_world();
    std::cout << "Process with world rank "<< rw << " is a group process!" << std::endl;
  }
  else if (master != nullptr)
  {
    // not sure yet if it makes sense to let the user control the master process
    int rw = group_process->get_rank_world();
    std::cout << "Process with world rank "<< rw << " is the MASTER OF THE UNIVERSE!" << std::endl;
  }
  else
  {
    // dom-debug
    std::cout << "Process with rank " << my_rank << " has no particular role, this should not happen." << std::endl;
  }
}

#include <iostream>
#include <kernel/base_header.hpp>
#include <kernel/universe.hpp>

using namespace std;

int main(int argc, char* argv[])
{
  // the following information will be read from some dat file

  // number of process groups (if not provided, then 1)
  int num_process_groups;
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

  //@Hilmar: ab hier scheint nichts mehr ausgefuehrt zu werden. Ich vermute das liegt daran, dass MPI_Init() nicht der
  //         erste Befehl ist der ausgefuehrt wird.
  cout << "Das hier erscheint nicht auf dem Bildschirm" << endl;

  // Get process objects. Note that on each process only one of the following three exists (the other two are
  // null pointers).
  LoadBalancer* load_balancer = universe.get_load_balancer();
  GroupProcess* group_process = universe.get_group_process();
  Master* master = universe.get_master();

//  // rank of this process
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//
//  cout << "Hi! I'm rank " << my_rank << endl;

  if (load_balancer != nullptr)
  {
    cout << "process with world rank "<< load_balancer->get_rank_world() << " is the load balancer with id "
         << load_balancer->get_group_id() << " and local rank " << load_balancer->get_rank_local() << endl;
    if (load_balancer->get_group_id() == 0)
    {
      load_balancer->read_mesh();
      load_balancer->create_work_groups();
    }
    else
    {
      // the second process group does something else...
    }
  }
  else if (group_process != nullptr)
  {
    // not sure yet if it makes sense to let the user control the group processes
    cout << "Process with world rank "<< group_process->get_rank_world() << " is a group process!" << endl;
  }
  else if (master != nullptr)
  {
    // not sure yet if it makes sense to let the user control the master process
    cout << "Process with world rank "<< master->get_rank_world() << " is the MASTER OF THE UNIVERSE!" << endl;
  }
  else
  {
    // dom-debug
    cout << "Process with rank " << my_rank << " has no particular role, this should not happen." << endl;
  }
}

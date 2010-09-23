// Test driver for universe

// includes, system
#include <iostream>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>

// main routine
// \author Hilmar Wobker
// \author Dominik Goeddeke
int main(int argc, char* argv[])
{

  // create universe with one process group
//  Universe* universe = Universe::create(argc, argv);

  // the following information will be read from some dat file

  // number of process groups (if not provided, then 1)
  int num_process_groups(2);
  // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
  int num_processes_in_group[] = {6, 2};
  // array of flags whether a dedicated load balancer process is needed in process groups
  // (must be provided when num_process_groups > 1)
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
      load_balancer->create_work_groups();

      // let all process with even world rank test the vector version of the function log_master_array()
      if (Process::rank % 2 == 0)
      {
        std::vector<std::string> messages(3);
        messages[0] = StringUtils::stringify(Process::rank) + ". process is testing...";
        messages[1] = "... this vector ...";
        messages[2] = "... logging feature!";
        Logger::log_master_array(messages, Logger::SCREEN);
      }

// COMMENT_HILMAR: TEMPORARY HACK to stop the infinite service loop of the master
      sleep(5);
      if (process_group->rank() == 0)
      {
        // the first MPI process of the process group tells the master to stop its service loop
        // reset buffer
        MPIUtils::buffer_pos = 0;
        // write the message id to the buffer
        int id = ServiceIDs::MASTER_FINISH_SERVICE;
        int mpi_error_code = MPI_Pack(&id, 1, MPI_INTEGER, MPIUtils::buffer, MPIUtils::BUFFERSIZE_BYTES,
                                      &MPIUtils::buffer_pos, MPI_COMM_WORLD);
        MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Pack");
        // send message
        mpi_error_code = MPI_Send(MPIUtils::buffer, MPIUtils::buffer_pos, MPI_PACKED, Process::rank_master, 0,
                                  MPI_COMM_WORLD);
        MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Send");
      }
// COMMENT_HILMAR: end of TEMPORARY HACK

    }
    else
    {
      // the second process group does something else, programmed by the application outside the kernel...
    }
  }
  else if(master != nullptr)
  {
    // This branch should only be entered when the infinite service loop of the master has been finished.
    // This, however, usually happens only program end.
  }
  else
  {
    MPIUtils::abort("Process with rank " + StringUtils::stringify(rank_world)
                    + " has no particular role, this should not happen.");
  }

  Universe::destroy();
}

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/util/assertion.hpp>
#include <test_system/test_system.hpp>
#include <kernel/comm.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>

using namespace TestSystem;
using namespace FEAST;

/**
* \brief testing creation of work groups based on a 2D mesh
*
* \test
* This test creates a universe, reads a 2D mesh, builds a base mesh and creates work groups.
*
* \author Hilmar Wobker
*/
template <typename Tag_, typename DT_, unsigned char world_dim_, unsigned char space_dim_>
class WorkGroupTest2D1
  : public TaggedTest<Tag_, DT_>
{

private:

  /// name of the mesh file to be read
  std::string _mesh_file;

public:

  /**
  * \brief CTOR
  *
  * \param[in] mesh_file
  * name of the mesh file to be read
  */
  WorkGroupTest2D1(std::string const mesh_file)
    : TaggedTest<Tag_, DT_>("cell_2d_quad_test_subdivision"),
      _mesh_file(mesh_file)
  {
  }

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
    LoadBalancer<space_dim_, world_dim_>* load_balancer,
    unsigned int& num_subgroups,
    unsigned int*& num_proc_in_subgroup,
    unsigned char*& group_contains_extra_coord,
    int**& subgroup_ranks,
    Graph**& graphs) const
  {
    CONTEXT("WorkGroupTest2D1::define_work_groups()");
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
    Logger::log_master("num_processes: " + stringify(num_processes) + "\nnum_cells: " + stringify(num_cells) + "\n");
    // assert that the number of processes is n
    ASSERT(num_processes == num_cells, "");
    // assert that no dedicated load balancer is used
    ASSERT(!load_balancer->group_has_dedicated_load_bal(), "");

    // set up the test case

    // no extra coordinator process is needed
    group_contains_extra_coord[0] = false;
    // set number of processes in the work group
    num_proc_in_subgroup[0] = num_cells;
    // work group ranks
    subgroup_ranks[0] = new int[num_proc_in_subgroup[0]];
    // set entries to {0, ..., n-1}
    for(unsigned int i(0) ; i < num_proc_in_subgroup[0] ; ++i)
    {
      subgroup_ranks[0][i] = i;
    }

    // get connectivity graph of the base mesh; this one will be used for the one and only work group
    graphs[0] = load_balancer->base_mesh()->graph();

  // COMMENT_HILMAR:
  // We assume here that each process receives exactly one BMC and that the index of the cell in the graph structure
  // equals the local rank within the work group. Later, there will be the matrix patch layer and the process patch
  // layer, which both have their own connectivity structure. Here, we actually need the connectivity graph of the
  // process patch layer.
  }


  /**
  * \brief sets the number of processes to use to 17 (16 BMCs + 1 master)
  *
  * For the semi-hard-wired example using a mesh with n base mesh cells we need n+1 processes in total
  * (n for the one and only process group and 1 for the master process).
  *
  * \todo As an intermediate hack, read the number of BMCs from the mesh file name and return a corresponding value.
  */
  unsigned long mpi_proc_count() const
  {
    return 17;
  }


  /// main routine
  void run() const
  {
    CONTEXT("WorkGroupTest2D1::run()");
    // init MPI
    MPIUtils::init_MPI();

    // COMMENT_HILMAR:
    // For this hard-wired example using the fixed 1D mesh defined in the file parser with 3 base mesh cells we need 4
    // processes in total (3 for the one and only process group and 1 for the master process).

    // set shortcut to the one and only instance of Universe (since this is the first call of
    // Universe<space_dim_, world_dim_>::instance(), it also calls the constructor of the Universe singleton class)
    Universe<space_dim_, world_dim_>* universe = Universe<space_dim_, world_dim_>::instance();

    try
    {
      universe->create("work_group_test_1d");
    }
    catch (Exception& e)
    {
      // abort the program
      ErrorHandler::exception_occured(e);
    }

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    LoadBalancer<space_dim_, world_dim_>* load_balancer = universe->load_balancer();
    Master* master = universe->master();

    int rank_world = Process::rank;

    if(load_balancer != nullptr)
    {
      ProcessGroup* process_group = load_balancer->process_group();

      // debug output
      int rank_process_group = process_group->rank();
      std::string s("Process " + stringify(rank_world) + " is the load balancer with local rank "
                    + stringify(rank_process_group) + ".");
      Logger::log(s);

      // let the load balancer read the mesh file
      load_balancer->read_mesh(_mesh_file);

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

// add some TEST_CHECK(...)

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
      MPIUtils::abort("Process with rank " + stringify(rank_world)
                      + " has no particular role, this should not happen.");
    }
  } // run()

}; // WorkGroupTest1D

// create test instance, using space and world dimension 2
WorkGroupTest2D1<Nil, Nil, 2, 2> work_group_test_2d_1("data/meshes/test_16bmc.feast");

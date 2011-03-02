// includes, system
#include <iostream>
#include <cstdlib> // for exit()

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <test_system/test_system.hpp>
#include <kernel/comm.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>
#include <kernel/base_mesh/cell_subdivision.hpp>

using namespace TestSystem;
using namespace FEAST;

/**
* \brief testing creation of work groups based on a 1D mesh
*
* \test
* This test creates a universe, reads a 1D mesh (currently hard coded), builds a base mesh and creates work groups.
*
* \author Hilmar Wobker
*/
template <typename Tag_, typename DT_, unsigned char world_dim_, unsigned char space_dim_>
class WorkGroupTest1D
  : public TaggedTest<Tag_, DT_>
{

public:

  /// CTOR
  WorkGroupTest1D()
    : TaggedTest<Tag_, DT_>("work_group_test_1d")
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
    CONTEXT("universe_test_1d_1: define_work_groups()");
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
    Logger::log_master("num_processes: " + StringUtils::stringify(num_processes) + "\nnum_cells: "
                       + StringUtils::stringify(num_cells) + "\n");
    // assert that the number of processes is n
    assert(num_processes == num_cells);
    // assert that no dedicated load balancer is used
    assert(!load_balancer->group_has_dedicated_load_bal());

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


  /// sets the number of processes to use to 4 (3 BMCs + 1 master)
  unsigned long mpi_proc_count() const
  {
    return 4;
  }


  /// main routine
  void run() const
  {
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
      universe->create();
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

      // let the load balancer "read" the mesh, which currently means: create a hard-wired mesh consisting of 3 edges
      load_balancer->read_mesh("dummy");

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

//      // get pointer to the base mesh
//      BaseMesh::BM<space_dim_, world_dim_>* bm = load_balancer->base_mesh();
//
//      // let the coordinator subdivide a cell
//      if(process_group->is_coordinator())
//      {
//       // subdivide cell 1
//        Logger::log("******************\n");
//        Logger::log("Subdividing cell 1\n");
//        Logger::log("******************\n");
//        BaseMesh::SubdivisionData<1, space_dim_, world_dim_>* subdiv_data
//          = new BaseMesh::SubdivisionData<1, space_dim_, world_dim_>(BaseMesh::CONFORM_SAME_TYPE);
//        bm->cell(1)->subdivide(subdiv_data);
//
//        // add created cells and subcells to the corresponding base mesh vectors
//        bm->add_created_items(bm->cell(1)->subdiv_data());
//        // set cell numbers (now they differ from indices)
//        bm->set_cell_numbers();
//        // print base mesh
//        bm->print(Logger::file);
//       // validate base mesh
//        bm->validate(Logger::file);
//
//        // TODO: neighbourhood update
//
//        Logger::log("!!! Neighbourhood update after subdivision not implemented yet!!!\n");
//        Logger::log("!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!\n");
//      }

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
  } // run()

}; // WorkGroupTest1D

// create test instance, using space and world dimension 1
WorkGroupTest1D<Nil, Nil, 1, 1> work_group_test_1d;

#include <kernel/base_header.hpp>

#include <kernel/util/dist.hpp>
#include <kernel/archs.hpp>
#include <kernel/foundation/psynch.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/util/exception.hpp>


#include <iostream>
#include <limits>

using namespace FEAT;
using namespace Foundation;

template<typename DT1_, typename DT2_, typename DT3_>
struct TestResult
{
  TestResult(DT1_ l, DT2_ r, DT3_ eps) :
    left(l),
    right(r),
    epsilon(eps)
  {
    //passed = std::abs(l - r) < eps ? true : false;
    passed = (l < r ? r - l : l - r) <= eps ? true : false;
  }

  TestResult()
  {
  }

  DT1_ left;
  DT2_ right;
  DT3_ epsilon;
  bool passed;
};

template<typename DT1_, typename DT2_, typename DT3_>
TestResult<DT1_, DT2_, DT3_> test_check_equal_within_eps(DT1_ l, DT2_ r, DT3_ eps)
{
  return TestResult<DT1_, DT2_, DT3_>(l, r, eps);
}

void check_psynch_presult()
{
#ifdef FEAT_HAVE_PARMETIS
  Dist::Comm comm(Dist::Comm::world());
  typename PExecutorParmetis<ParmetisModePartKway>::PResult local(1, comm);
  local.get()[0] = comm.rank() == 0 ? 1 : 0;

  local.get_vtxdist()[0] = 0;
  local.get_vtxdist()[1] = 1;
  local.get_vtxdist()[2] = 2;

  auto global( PSynch<PExecutorParmetis<ParmetisModePartKway> >::exec(local, 2) );

  TestResult<int, int, int>* res = new TestResult<int, int, int>[2];
  res[0] = test_check_equal_within_eps(global.get()[0], 1, 0);
  res[1] = test_check_equal_within_eps(global.get()[1], 0, 0);

  bool passed(true);
  for(unsigned long i(0) ; i < 2 ; ++i)
    if(!res[i].passed)
    {
      std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "! (psynch) on rank " << comm.rank() << std::endl;
      passed = false;
      break;
    }

  if(passed)
    std::cout << "PASSED (rank " << comm.rank() <<"): foundation_psystem_test (psynch)" << std::endl;

  delete[] res;
#endif
}

void check_psynch_stringstream()
{
#ifdef FEAT_HAVE_PARMETIS
  Dist::Comm comm(Dist::Comm::world());
  std::stringstream synchstream;

  if(comm.rank() == 0)
    synchstream << "MASTER";

  comm.bcast_stringstream(synchstream);

  if(synchstream.str() == "MASTER")
    std::cout << "PASSED (rank " << comm.rank() <<"): foundation_psystem_test (psynch stringstream)" << std::endl;
  else
  {
      std::cout << "FAILED: (psynch stringstream) on rank " << comm.rank() << " is '" << synchstream.str() << "'" << std::endl;
  }
#endif
}

void check_psynch_meshstreamer()
{
#ifdef FEAT_HAVE_PARMETIS
  Dist::Comm comm(Dist::Comm::world());

  std::stringstream synchstream;
  std::string file_prefix(FEAT_SRC_DIR);
  std::string filename(file_prefix);
  filename += "/data/meshes/flowbench_c2d_00_quad_130.xml";

  if(comm.rank() == 0)
  {
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if(!ifs.is_open())
    {
      throw FileNotFound(filename);
    }
    synchstream << ifs.rdbuf();
  }

  comm.bcast_stringstream(synchstream);
  Geometry::MeshFileReader mesh_file_reader(synchstream);

  // create an empty atlas and a root mesh node
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> MeshType;
  Geometry::MeshAtlas<MeshType>* atlas = new Geometry::MeshAtlas<MeshType>();
  Geometry::RootMeshNode<MeshType>* node = new Geometry::RootMeshNode<MeshType>(nullptr, atlas);

  mesh_file_reader.parse(*node, *atlas);

  std::cout << "PASSED (rank " << comm.rank() <<"): foundation_psystem_test (psynch meshstreamer)" << std::endl;

  delete node;
  delete atlas;
#endif
}

void check_pexecutor_rank_at_elem()
{
#ifdef FEAT_HAVE_PARMETIS
  Dist::Comm comm(Dist::Comm::world());
  typename PExecutorParmetis<ParmetisModePartKway>::PResult local(1, comm);
  local.get()[0] = comm.rank() == 0 ? 1 : 0;

  local.get_vtxdist()[0] = 0;
  local.get_vtxdist()[1] = 1;
  local.get_vtxdist()[2] = 2;

  auto global( PSynch<PExecutorParmetis<ParmetisModePartKway> >::exec(local, 2) );

  TestResult<int, int, int>* res = new TestResult<int, int, int>[2];
  res[0] = test_check_equal_within_eps(global.get()[0], 1, 0);
  res[1] = test_check_equal_within_eps(global.get()[1], 0, 0);

  bool passed(true);
  for(unsigned long i(0) ; i < 2 ; ++i)
    if(!res[i].passed)
    {
      std::cout << "FAILED: " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "! (pexecutor, rank_at_elem) on rank " << comm.rank() << std::endl;
      passed = false;
      break;
    }

  Adjacency::Graph graph(global.rank_at_element());

  if(passed)
    std::cout << "PASSED (rank " << comm.rank() <<"): foundation_psystem_test (pexecutor, rank_at_elem)" << std::endl;

  delete[] res;
#endif
}

int main(int argc, char* argv[])
{
  int me(0);
#ifdef FEAT_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
  (void)argc;
  (void)argv;
  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;

#ifdef FEAT_HAVE_MPI
  check_psynch_presult();
  check_psynch_stringstream();
  check_psynch_meshstreamer();
  check_pexecutor_rank_at_elem();
#else
  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
#endif

#ifdef FEAT_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

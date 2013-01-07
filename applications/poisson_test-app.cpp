#ifndef SERIAL
#include <mpi.h>
#endif

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include<kernel/foundation/mesh_control.hpp>
#include<kernel/foundation/dense_data_wrapper.hpp>
#include<kernel/geometry/macro_factory.hpp>
#include<kernel/geometry/patch_factory.hpp>
#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/vector_mirror.hpp>
#include<kernel/foundation/halo_control.hpp>
#include<kernel/foundation/halo.hpp>
#include<kernel/geometry/cell_sub_set.hpp>
#include<kernel/archs.hpp>
#include<deque>
#include<algorithm>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/space/dof_mirror.hpp>
#include <kernel/assembly/standard_operators.hpp>
#include <kernel/assembly/standard_functionals.hpp>
#include <kernel/assembly/dirichlet_bc.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;
using namespace FEAST::Foundation;
using namespace FEAST::Geometry;

int main(int argc, char* argv[])
{
  int me(0);

#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif


#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}

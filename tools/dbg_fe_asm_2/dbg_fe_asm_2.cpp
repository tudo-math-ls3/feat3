#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/util/stop_watch.hpp>
#include <cstdio>

using namespace FEAST;
using namespace FEAST::Geometry;

typedef ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef MeshStreamerFactory<QuadMesh> QuadMeshFactory;
typedef StandardRefinery<QuadMesh> QuadRefinery;
typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceQ1T;


template<typename Space_>
void test_asm(const Space_& space, const String& cubature_name, int imat)
{
  std::cout << "Assembling matrix structure..." << std::endl;

  Adjacency::Graph dof_adjacency(Assembly::SymbolicGraphAssembler<>::assemble_graph(space));
  std::cout << "NEQ : " << dof_adjacency.get_num_nodes_domain() << std::endl;
  std::cout << "NNZE: " << dof_adjacency.get_num_indices() << std::endl;

  StopWatch stop_watch;

  LAFEM::SparseMatrixCSR<Mem::Main, double> matrix_d;
  Assembly::SymbolicMatrixAssemblerBase::assemble(matrix_d, dof_adjacency);

  Cubature::DynamicFactory cubature_factory(cubature_name);

  std::cout << "Assembling " << (imat == 0 ? "stiffness" : "mass") << " matrix data..." << std::endl;
  matrix_d.clear(0.0f);
  stop_watch.reset();
  stop_watch.start();
  switch(imat)
  {
  case 0:
    {
      Assembly::Common::LaplaceOperator laplace;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_d, laplace, space, cubature_factory);
    }
    break;

  case 1:
    {
      Assembly::Common::IdentityOperator identity;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_d, identity, space, cubature_factory);
    }
    break;
  }
  stop_watch.stop();
  std::cout << "Assembly time: " << stop_watch.elapsed() << " seconds" << std::endl;
}


int main(int argc, char* argv[])
{
  int nrefs = -1;
  int imat = 0;
  int space = 0;

  if(argc < 2)
  {
    std::cout << std::endl << "  USAGE: dbg-fe-asm-2 [-q1|-q1t] [-sm|-mm] [-r:<count>] <mesh-name>" << std::endl << std::endl;
    std::cout << "-r:<n>            Refine mesh <n> times" << std::endl;
    std::cout << "-q1               Use Q1 element" << std::endl;
    std::cout << "-q1t              Use Q1~ element" << std::endl;
    std::cout << "-sm               Assemble stiffness matrix" << std::endl;
    std::cout << "-mm               Assemble mass matrix" << std::endl;
    std::cout << std::endl;
    return 0;
  }

  // loop over all arguments
  for(int i(1); i < argc-1; ++i)
  {
    String arg(argv[i]);
    if(arg.compare("-q1") == 0)
      space = 0;
    else if(arg.compare("-q1t") == 0)
      space = 1;
    else if(arg.compare("-sm") == 0)
      imat = 0;
    else if(arg.compare("-mm") == 0)
      imat = 1;
    else if(arg.substr(0,3).compare("-r:") == 0)
      nrefs = atoi(arg.substr(3).c_str());
    else
    {
      std::cout << "Invalid argument: '" << argv[i] << "'" << std::endl;
      return 1;
    }
  }

  QuadMesh* mesh = nullptr;

  // read input mesh
  {
    String mesh_name = String("data/meshes/") + String(argv[argc-1]) + String(".txt");
    std::cout << "Reading mesh from '" << mesh_name << "'..." << std::endl;
    MeshStreamer mesh_reader;
    mesh_reader.parse_mesh_file(mesh_name);

    // create mesh
    QuadMeshFactory mesh_factory(mesh_reader);
    mesh = new QuadMesh(mesh_factory);
  }

  // refine if desired
  if(nrefs > 0)
  {
    std::cout << "Refining mesh up to level " << nrefs << "..." << std::endl;
    for(int i = 0; i < nrefs; ++i)
    {
      QuadMesh* mesh2 = mesh;
      {
        Geometry::StandardRefinery<QuadMesh> refinery(*mesh2);
        mesh = new QuadMesh(refinery);
      }
      delete mesh2;
    }
  }

  std::cout << "Creating trafo..." << std::endl;

  QuadTrafo* trafo = new QuadTrafo(*mesh);

  switch(space)
  {
  case 0:
    {
      std::cout << "Creating Q1 space..." << std::endl;
      QuadSpaceQ1 space_q1(*trafo);
      test_asm(space_q1, "gauss-legendre:2", imat);
    }
    break;

  case 1:
    {
      std::cout << "Creating Q1~ space..." << std::endl;
      QuadSpaceQ1T space_q1t(*trafo);
      test_asm(space_q1t, "gauss-legendre:3", imat);
    }
    break;
  }

  delete trafo;
  delete mesh;
}

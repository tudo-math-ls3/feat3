#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
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
typedef StandardRefinery<QuadMesh> QuadRefinery;
typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
typedef Space::Lagrange2::Element<QuadTrafo> QuadSpaceQ2;
typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceQ1T;

int nrefs = -1;
bool verify_d = false;
bool verify_s = false;
int imat = 0;


QuadMesh* read_quad_mesh(const char* filename)
{
  FILE* file = fopen(filename, "rb");
  ASSERT_(file != NULL);

  int magic, ne[3], idx[4];
  double vtx[2];

  // read magic
  size_t status = fread(&magic, 4, 1, file);
  if (status != 1)
    throw InternalError("fread failure!");
  ASSERT_(magic == 0x48534D42); // "BMSH"

  // read entity count
  status = fread(ne, 12, 1, file);
  if (status != 1)
    throw InternalError("fread failure!");
  Index num_entities[3] = {Index(ne[0]), Index(ne[1]), Index(ne[2])};

  // create mesh
  QuadMesh* mesh = new QuadMesh(num_entities);

  // read vertex coords
  for(Index i(0); i < Index(ne[0]); ++i)
  {
    status = fread(vtx, 16, 1, file);
    if (status != 1)
      throw InternalError("fread failure!");
    mesh->get_vertex_set()[i][0] = vtx[0];
    mesh->get_vertex_set()[i][1] = vtx[1];
  }

  // read vertex-at-edge
  for(Index i(0); i < Index(ne[1]); ++i)
  {
    status = fread(idx, 8, 1, file);
    if (status != 1)
      throw InternalError("fread failure!");
    mesh->get_index_set<1,0>()[i][0] = Index(idx[0]);
    mesh->get_index_set<1,0>()[i][1] = Index(idx[1]);
  }

  // read vertex-at-quad
  for(Index i(0); i < Index(ne[2]); ++i)
  {
    status = fread(idx, 16, 1, file);
    if (status != 1)
      throw InternalError("fread failure!");
    mesh->get_index_set<2,0>()[i][0] = Index(idx[0]);
    mesh->get_index_set<2,0>()[i][1] = Index(idx[1]);
    mesh->get_index_set<2,0>()[i][2] = Index(idx[3]);
    mesh->get_index_set<2,0>()[i][3] = Index(idx[2]);
  }

  // read edge-at-quad
  for(Index i(0); i < Index(ne[2]); ++i)
  {
    status = fread(idx, 16, 1, file);
    if (status != 1)
      throw InternalError("fread failure!");
    mesh->get_index_set<2,1>()[i][0] = Index(idx[0]);
    mesh->get_index_set<2,1>()[i][1] = Index(idx[2]);
    mesh->get_index_set<2,1>()[i][2] = Index(idx[3]);
    mesh->get_index_set<2,1>()[i][3] = Index(idx[1]);
  }

  // okay
  fclose(file);

  return mesh;
}

bool verify_struct(const char* filename, const Adjacency::Graph& graph)
{
  FILE* file = fopen(filename, "rb");
  ASSERT_(file != NULL);

  int magic;
  unsigned int ne[2],idx;
  bool okay = false;

  // fetch graph data
  Index neq = graph.get_num_nodes_domain();
  Index nnze = graph.get_num_indices();
  const Index* row_ptr = graph.get_domain_ptr();
  const Index* col_idx = graph.get_image_idx();

  // read magic
  size_t status = fread(&magic, 4, 1, file);
  if (status != 1)
    throw InternalError("fread failure!");
  ASSERT_(magic == 0x58544D42); // "BMTX"

  // read entity count
  status = fread(ne, 8, 1, file);
  if (status != 1)
    throw InternalError("fread failure!");

  // validate counts
  if(neq != ne[0])
  {
    std::cout << "ERROR: neq mismatch!" << std::endl;
    goto _end;
  }
  if(nnze != ne[1])
  {
    std::cout << "ERROR: nnze mismatch!" << std::endl;
    goto _end;
  }

  // skip assembly time
  fseek(file, 4, SEEK_CUR);

  // validate row pointer
  for(Index i(0); i <= neq; ++i)
  {
    status = fread(&idx, 4, 1, file);
    if (status != 1)
      throw InternalError("fread failure!");
    if(idx != row_ptr[i])
    {
      std::cout << "ERROR: row_ptr[" << i << "] mismatch!" << std::endl;
      goto _end;
    }
  }

  // validate column index
  for(Index i(0); i < nnze; ++i)
  {
    status = fread(&idx, 4, 1, file);
    if (status != 1)
      throw InternalError("fread failure!");
    if(idx != col_idx[i])
    {
      std::cout << "ERROR: col_idx[" << i << "] mismatch!" << std::endl;
      goto _end;
    }
  }

  okay = true;

_end:
  fclose(file);
  return okay;
}

template<typename DataType_>
double verify_data(const char* filename, const LAFEM::SparseMatrixCSR<Mem::Main, DataType_>& matrix, double& asm_time)
{
  FILE* file = fopen(filename, "rb");
  ASSERT_(file != NULL);

  int magic;
  unsigned int ne[2];
  float fasm_time = 0.0f;
  double dx,da,db;

  double dev_rel = 0.0;

  // fetch matrix data
  Index neq = matrix.rows();
  Index nnze = matrix.used_elements();
  const Index* row_ptr = matrix.row_ptr();
  const DataType_* data = matrix.val();

  // read magic
  size_t status = fread(&magic, 4, 1, file);
  if (status != 1)
    throw InternalError("fread failure!");
  ASSERT_(magic == 0x58544D42); // "BMTX"

  // read entity count
  status = fread(ne, 8, 1, file);
  if (status != 1)
    throw InternalError("fread failure!");

  // validate counts
  if(neq != ne[0])
  {
    std::cout << "ERROR: neq mismatch!" << std::endl;
    goto _end;
  }
  if(nnze != ne[1])
  {
    std::cout << "ERROR: nnze mismatch!" << std::endl;
    goto _end;
  }

  // read assembly time
  status = fread(&fasm_time, 4, 1, file);
  if (status != 1)
    throw InternalError("fread failure!");
  asm_time = fasm_time;

  // skip row pointer
  fseek(file, long(4*(ne[0]+1)), SEEK_CUR);

  // skip column indices
  fseek(file, long(4*ne[1]), SEEK_CUR);

  // validate data
  for(Index i(0); i < neq; ++i)
  {
    da = db = 0.0;
    for(Index j(row_ptr[i]); j < row_ptr[i+1]; ++j)
    {
      status = fread(&dx, 8, 1, file);
      if (status != 1)
        throw InternalError("fread failure!");
      da = std::max(da, std::abs(dx));                   // ||A_i.||
      db = std::max(db, std::abs(dx - double(data[j]))); // ||A_i. - B_i.||
    }
    dev_rel = std::max(dev_rel, db / da);
  }

_end:
  fclose(file);
  return dev_rel;
}

template<typename Space_>
void test_asm(const Space_& space, const String& cubature_name, const char* matx_name)
{
  std::cout << "Assembling matrix structure..." << std::endl;

  Adjacency::Graph dof_adjacency(Assembly::SymbolicGraphAssembler<>::assemble_graph(space));
  std::cout << "NEQ : " << dof_adjacency.get_num_nodes_domain() << std::endl;
  std::cout << "NNZE: " << dof_adjacency.get_num_indices() << std::endl;

  if(verify_s)
  {
    std::cout << "Verifying matrix structure..." << std::endl;
    if(!verify_struct(matx_name, dof_adjacency))
      return;
  }

  StopWatch stop_watch;

  LAFEM::SparseMatrixCSR<Mem::Main, double> matrix_d;
  Assembly::SymbolicMatrixAssemblerBase::assemble(matrix_d, dof_adjacency);

  Cubature::DynamicFactory cubature_factory(cubature_name);

  std::cout << "Assembling double precision matrix data..." << std::endl;
  matrix_d.format(0.0f);
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
  std::cout << "Assembly time for double precision: " << stop_watch.elapsed() << " seconds" << std::endl;

  if(verify_d)
  {
    double feat_time(0.0);
    std::cout << "Validating matrix against '" << matx_name << "'..." << std::endl;
    std::cout << "Deviation: " << verify_data(matx_name, matrix_d, feat_time) << std::endl;
    std::cout << "Feat2 Time: " << feat_time << " seconds" << std::endl;
  }
}


int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    std::cout << std::endl << "  USAGE: dbg-fe-asm-1 [-v[s|d]] [-q1[t]] [-sm|mm] [-r:<count>] <base-name>"
      << std::endl << std::endl;
    std::cout << "Valid Options:" << std::endl;
    std::cout << "-vd        verify matrix data" << std::endl;
    std::cout << "-vs        verify matrix structure" << std::endl;
    std::cout << "-v         verify both matrix structure and data" << std::endl;
    std::cout << "-q1        use Q1 element" << std::endl;
    std::cout << "-q2        use Q2 element" << std::endl;
    std::cout << "-q1t       use Q1T element" << std::endl;
    std::cout << "-sm        assemble stiffness matrix" << std::endl;
    std::cout << "-mm        assemble mass matrix" << std::endl;
    std::cout << "-r:<n>     refine mesh <n> times" << std::endl;
    std::cout << std::endl;
    return 0;
  }

  int space = 0;

  // loop over all arguments
  for(int i(1); i < argc-1; ++i)
  {
    String arg(argv[i]);
    if(arg.compare("-v") == 0)
      verify_d = verify_s = true;
    else if(arg.compare("-vd") == 0)
      verify_d = true;
    else if(arg.compare("-vs") == 0)
      verify_s = true;
    else if(arg.compare("-q1") == 0)
      space = 0;
    else if(arg.compare("-q2") == 0)
      space = 1;
    else if(arg.compare("-q1t") == 0)
      space = 2;
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

  if(nrefs >= 0)
    verify_d = verify_s = false;

  const char* base_name = argv[argc-1];
  char mesh_name[256], matx_name[256];
  sprintf(mesh_name, "data/%s_mesh.bin", base_name);
  switch(space)
  {
  case 0:
    // Q1 element
    switch(imat)
    {
    case 0:
      sprintf(matx_name, "data/%s_laplace_q1_g2.bin", base_name);
      break;
    case 1:
      sprintf(matx_name, "data/%s_mass_q1_g2.bin", base_name);
      break;
    }
    break;

  case 1:
    // Q2 element
    switch(imat)
    {
    case 0:
      sprintf(matx_name, "data/%s_laplace_q2_g3.bin", base_name);
      break;
    case 1:
      sprintf(matx_name, "data/%s_mass_q2_g3.bin", base_name);
      break;
    }
    break;

  case 2:
    // Q1~ element
    switch(imat)
    {
    case 0:
      sprintf(matx_name, "data/%s_laplace_q1t_g3.bin", base_name);
      break;
    case 1:
      sprintf(matx_name, "data/%s_mass_q1t_g3.bin", base_name);
      break;
    }
    break;
  }

  // read input mesh
  std::cout << "Reading mesh from '" << mesh_name << "'..." << std::endl;
  QuadMesh* mesh = read_quad_mesh(mesh_name);

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
      test_asm(space_q1, "gauss-legendre:2", matx_name);
    }
    break;

  case 1:
    {
      std::cout << "Creating Q2 space..." << std::endl;
      QuadSpaceQ2 space_q2(*trafo);
      test_asm(space_q2, "gauss-legendre:3", matx_name);
    }
    break;

  case 2:
    {
      std::cout << "Creating Q1~ space..." << std::endl;
      QuadSpaceQ1T space_q1t(*trafo);
      test_asm(space_q1t, "gauss-legendre:3", matx_name);
    }
    break;
  }

  delete trafo;
  delete mesh;
}

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/lafem/matrix_mirror.hpp>

using namespace FEAST;

typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef Geometry::CellSubSet<Shape::Quadrilateral> QuadCellSet;
typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceQ1T;

typedef LAFEM::DenseVector<Mem::Main, double> MyVectorType;
typedef LAFEM::SparseMatrixCSR<Mem::Main, double> MyMatrixType;
typedef LAFEM::SparseMatrixCSR<Mem::Main, double> BufMatrixType;

typedef LAFEM::VectorMirror<Mem::Main, double> VectorMirrorType;
typedef LAFEM::MatrixMirror<VectorMirrorType> MatrixMirrorType;

void fill_cell_set(QuadCellSet& cell, int face);
void fill_quad_mesh_2d(QuadMesh& mesh, Real x = 0.0, Real y = 0.0);

void sync_add(MyVectorType& a, MyVectorType& b)
{
  double * av(a.elements());
  double * bv(b.elements());
  Index n = a.size();
  for(Index i(0); i < n; ++i)
    av[i] = bv[i] = av[i] + bv[i];
}

void sync_add(BufMatrixType& a, BufMatrixType& b)
{
  double * av(a.val());
  double * bv(b.val());
  Index nnze = a.used_elements();
  for(Index i(0); i < nnze; ++i)
    av[i] = bv[i] = av[i] + bv[i];
}

template<typename Space_, typename CellSet_>
void test_mirror(
  const Space_& space_0,
  const Space_& space_1,
  const CellSet_& cell_0,
  const CellSet_& cell_1,
  const String cubature_name = "gauss-legendre:2")
{
  // assemble dof-adjacency graphs
  Adjacency::Graph dof_adj_0(Assembly::SymbolicGraphAssembler<>::assemble_graph(space_0));
  Adjacency::Graph dof_adj_1(Assembly::SymbolicGraphAssembler<>::assemble_graph(space_1));

  // assemble dof-mirror graphs
  Adjacency::Graph dof_mir_0(Assembly::MirrorAssembler::assemble_mirror_graph(space_0, cell_0));
  Adjacency::Graph dof_mir_1(Assembly::MirrorAssembler::assemble_mirror_graph(space_1, cell_1));

  // allocate local vectors
  MyVectorType loc_vec_0(space_0.get_num_dofs(), Real(0));
  MyVectorType loc_vec_1(space_1.get_num_dofs(), Real(0));

  // allocate local matrices
  MyMatrixType loc_mat_0;
  MyMatrixType loc_mat_1;
  Assembly::SymbolicMatrixAssemblerBase::assemble(loc_mat_0, dof_adj_0);
  Assembly::SymbolicMatrixAssemblerBase::assemble(loc_mat_1, dof_adj_1);

  Cubature::DynamicFactory cubature_factory(cubature_name);

  // assemble local vectors
  Assembly::Common::ConstantFunction rhs_func(1.0);
  Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> rhs_functional(rhs_func);
  Assembly::LinearFunctionalAssembler::assemble_vector(loc_vec_0, rhs_functional, space_0, cubature_factory, 16.0);
  Assembly::LinearFunctionalAssembler::assemble_vector(loc_vec_1, rhs_functional, space_1, cubature_factory, 16.0);

  // assemble local matrices
  loc_mat_0.format();
  loc_mat_1.format();
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(loc_mat_0, laplace, space_0, cubature_factory, 6.0);
  Assembly::BilinearOperatorAssembler::assemble_matrix1(loc_mat_1, laplace, space_1, cubature_factory, 6.0);

  // create vector mirrors
  VectorMirrorType vec_mir_0;
  VectorMirrorType vec_mir_1;
  Assembly::MirrorAssembler::assemble_mirror(vec_mir_0, dof_mir_0);
  Assembly::MirrorAssembler::assemble_mirror(vec_mir_1, dof_mir_1);

  // create matrix mirrors
  MatrixMirrorType mat_mir_0(vec_mir_0, vec_mir_0);
  MatrixMirrorType mat_mir_1(vec_mir_1, vec_mir_1);

  // allocate buffer vectors
  MyVectorType buf_vec_0;
  MyVectorType buf_vec_1;
  Assembly::MirrorAssembler::assemble_buffer_vector(buf_vec_0, vec_mir_0, loc_vec_0);
  Assembly::MirrorAssembler::assemble_buffer_vector(buf_vec_1, vec_mir_1, loc_vec_1);

  // allocate buffer matrices
  BufMatrixType buf_mat_0;
  BufMatrixType buf_mat_1;
  Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat_0, mat_mir_0, loc_mat_0);
  Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat_1, mat_mir_1, loc_mat_1);

  // gather vectors
  vec_mir_0.gather_dual<Algo::Generic>(buf_vec_0, loc_vec_0);
  vec_mir_1.gather_dual<Algo::Generic>(buf_vec_1, loc_vec_1);

  // gather matrices
  mat_mir_0.gather<Algo::Generic>(buf_mat_0, loc_mat_0);
  mat_mir_1.gather<Algo::Generic>(buf_mat_1, loc_mat_1);

  // synchronise vectors
  sync_add(buf_vec_0, buf_vec_1);

  // synchronise matrices
  sync_add(buf_mat_0, buf_mat_1);

  // clone local vectors
  MyVectorType glob_vec_0;
  glob_vec_0.clone(loc_vec_0);
  MyVectorType glob_vec_1;
  glob_vec_1.clone(loc_vec_1);

  // clone local matrices
  MyMatrixType glob_mat_0;
  glob_mat_0.clone(loc_mat_0);
  MyMatrixType glob_mat_1;
  glob_mat_1.clone(loc_mat_1);

  // scatter vectors
  vec_mir_0.scatter_dual<Algo::Generic>(glob_vec_0, buf_vec_0);
  vec_mir_1.scatter_dual<Algo::Generic>(glob_vec_1, buf_vec_1);

  // scatter matrices
  mat_mir_0.scatter<Algo::Generic>(glob_mat_0, buf_mat_0);
  mat_mir_1.scatter<Algo::Generic>(glob_mat_1, buf_mat_1);

  // print output
  std::cout << "Vector-Local :0" << std::endl << loc_vec_0 << std::endl;
  std::cout << "Vector-Local :1" << std::endl << loc_vec_1 << std::endl;
  std::cout << "Vector-Buffer:0" << std::endl << buf_vec_0 << std::endl;
  std::cout << "Vector-Buffer:1" << std::endl << buf_vec_1 << std::endl;
  std::cout << "Vector-Global:0" << std::endl << glob_vec_0 << std::endl;
  std::cout << "Vector-Global:1" << std::endl << glob_vec_1 << std::endl;
  std::cout << "Matrix-Local :0" << std::endl << loc_mat_0 << std::endl;
  std::cout << "Matrix-Local :1" << std::endl << loc_mat_1 << std::endl;
  std::cout << "Matrix-Buffer:0" << std::endl << buf_mat_0 << std::endl;
  std::cout << "Matrix-Buffer:1" << std::endl << buf_mat_1 << std::endl;
  std::cout << "Matrix-Global:0" << std::endl << glob_mat_0 << std::endl;
  std::cout << "Matrix-Global:1" << std::endl << glob_mat_1 << std::endl;
}

template<typename Space_>
void test_it(QuadMesh& mesh_0, QuadMesh& mesh_1, const QuadCellSet& cell_0, const QuadCellSet& cell_1)
{
  // create trafos
  QuadTrafo trafo_0(mesh_0);
  QuadTrafo trafo_1(mesh_1);

  // create spaces
  Space_ space_0(trafo_0);
  Space_ space_1(trafo_1);

  // test vector mirror
  test_mirror(space_0, space_1, cell_0, cell_1);
}

int main(int /*argc*/, char** /*argv*/)
{
  static const Index num_entities[] =
  {
    4, 4, 1
  };
  static const Index num_cellset_entities[] =
  {
    2, 1, 0
  };

  // create and fill meshes
  QuadMesh mesh_0(num_entities);
  QuadMesh mesh_1(num_entities);
  fill_quad_mesh_2d(mesh_0, -1.0);
  fill_quad_mesh_2d(mesh_1,  0.0);

  // create cell sets
  QuadCellSet cell_0(num_cellset_entities);
  QuadCellSet cell_1(num_cellset_entities);
  fill_cell_set(cell_0, 3);
  fill_cell_set(cell_1, 2);

  typedef QuadSpaceQ1 QuadSpace;

  // test coarse
  test_it<QuadSpace>(mesh_0, mesh_1, cell_0, cell_1);

  // refine meshes
  Geometry::StandardRefinery<QuadMesh> mesh_refinery0(mesh_0), mesh_refinery1(mesh_1);
  QuadMesh mesh_0f(mesh_refinery0);
  QuadMesh mesh_1f(mesh_refinery1);
  // refine cell sets
  Geometry::StandardRefinery<QuadCellSet, QuadMesh> cell_refinery0(cell_0, mesh_0), cell_refinery1(cell_1, mesh_1);
  QuadCellSet cell_0f(cell_refinery0);
  QuadCellSet cell_1f(cell_refinery1);

  // test refined
  test_it<QuadSpace>(mesh_0f, mesh_1f, cell_0f, cell_1f);
}


void fill_cell_set(QuadCellSet& cell, int face)
{
  switch(face)
  {
  case 0:
    cell.get_target_set<0>()[0] = 0;
    cell.get_target_set<0>()[1] = 1;
    cell.get_target_set<1>()[0] = 0;
    break;

  case 1:
    cell.get_target_set<0>()[0] = 2;
    cell.get_target_set<0>()[1] = 3;
    cell.get_target_set<1>()[0] = 1;
    break;

  case 2:
    cell.get_target_set<0>()[0] = 0;
    cell.get_target_set<0>()[1] = 2;
    cell.get_target_set<1>()[0] = 2;
    break;

  case 3:
    cell.get_target_set<0>()[0] = 1;
    cell.get_target_set<0>()[1] = 3;
    cell.get_target_set<1>()[0] = 3;
    break;
  }
}

void fill_quad_mesh_2d(QuadMesh& mesh, Real x, Real y)
{
  // set up vertex coordinates array
  const Real vtx0[4*2] =
  {
    0.0+x, 0.0+y,
    1.0+x, 0.0+y,
    0.0+x, 1.0+y,
    1.0+x, 1.0+y
  };

  // set up vertices-at-edge array
  static const Index v_e0[4*2] =
  {
    0, 1,
    2, 3,
    0, 2,
    1, 3
  };

  // set up vertices-at-quad array
  static const Index v_q0[1*4] =
  {
    0, 1, 2, 3
  };

  // set up edges-at-quad array
  static const Index e_q0[1*4] =
  {
    0,  1,  2,  3
  };

  Geometry::TestAux::copy_vtx(mesh.get_vertex_set(), vtx0);
  Geometry::TestAux::copy_idx(mesh.get_index_set<1,0>(), v_e0);
  Geometry::TestAux::copy_idx(mesh.get_index_set<2,0>(), v_q0);
  Geometry::TestAux::copy_idx(mesh.get_index_set<2,1>(), e_q0);

} // create_quad_mesh_2d

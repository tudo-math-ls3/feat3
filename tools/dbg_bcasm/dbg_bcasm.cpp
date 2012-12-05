#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/assembly/standard_operators.hpp>
#include <kernel/assembly/standard_functionals.hpp>
#include <kernel/assembly/dirichlet_bc.hpp>

using namespace FEAST;

typedef double DataType;

typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef Geometry::CellSubSet<Shape::Quadrilateral> QuadCellSet;
typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceQ1T;

typedef LAFEM::DenseVector<Mem::Main, DataType> VectorType;
typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType> MatrixType;

typedef LAFEM::UnitFilter<Mem::Main, DataType> UnitFilterType;

void fill_cell_set(QuadCellSet& cell, int face);
void fill_quad_mesh_2d(QuadMesh& mesh, Real x = 0.0, Real y = 0.0);

template<typename T_>
class RhsFunc
{
public:
  static T_ eval(T_ /*x*/, T_ /*y*/)
  {
    return T_(1);
  }
};

template<typename T_>
class BCFunc
{
public:
  static T_ eval(T_ /*x*/, T_ /*y*/)
  {
    return T_(17);
  }
};

template<typename Space_, typename CellSet_>
void test_bcasm(
  const Space_& space,
  const CellSet_& cell,
  const String cubature = "gauss-legendre:2")
{
  // assemble system matrix
  MatrixType mat_sys(Space::DofAdjacency<>::assemble(space));
  mat_sys.clear();
  Assembly::BilinearScalarLaplaceFunctor::assemble(mat_sys, space, cubature);

  // assemble rhs vector
  VectorType vec_rhs(space.get_num_dofs(), Real(0));
  Assembly::LinearScalarIntegralFunctor<RhsFunc>::assemble(vec_rhs, space, cubature);

  // allocate solution vector
  VectorType vec_sol(space.get_num_dofs(), Real(1));

  // assemble homogene Dirichlet BCs
  Assembly::DirichletBC<Space_> dirichlet(space);
  dirichlet.add_cell_set(cell);

  // assemble filter:
  // a) homogene Dirichlet BCs
  //UnitFilterType filter(dirichlet.assemble<Mem::Main, DataType>());
  // b) inhomogene Dirichlet BCs
  Analytic::StaticWrapperFunctor<BCFunc> bc_func;
  UnitFilterType filter(dirichlet.assemble<Mem::Main, DataType>(bc_func));

  // filter system
  filter.filter_mat(mat_sys);
  filter.filter_rhs(vec_rhs);
  filter.filter_sol(vec_sol);

  // print output
  std::cout << "System Matrix:" << std::endl << mat_sys << std::endl;
  std::cout << "Rhs Vector:" << std::endl << vec_rhs << std::endl;
  std::cout << "Sol Vector:" << std::endl << vec_sol << std::endl;
}

template<typename Space_>
void test_it(const QuadMesh& mesh_0, const QuadMesh& mesh_1, const QuadCellSet& cell_0, const QuadCellSet& cell_1)
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
  QuadMesh mesh(num_entities);
  fill_quad_mesh_2d(mesh,  0.0);

  // create cell sets
  QuadCellSet cell(num_cellset_entities);
  fill_cell_set(cell, 3);

  // create trafo
  QuadTrafo trafo(mesh);

  // create space
  typedef QuadSpaceQ1 QuadSpace;
  QuadSpace space(trafo);

  // test coarse
  test_bcasm(space, cell);

  if(true)
  {
    // refine meshes
    Geometry::StandardRefinery<QuadMesh> mesh_refinery(mesh);
    QuadMesh mesh_f(mesh_refinery);
    // refine cell sets
    Geometry::StandardRefinery<QuadCellSet, QuadMesh> cell_refinery(cell, mesh);
    QuadCellSet cell_f(cell_refinery);
    // create trafo
    QuadTrafo trafo_f(mesh_f);

    // create space
    QuadSpace space_f(trafo_f);

    // test refined
    test_bcasm(space_f, cell_f);
  }
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

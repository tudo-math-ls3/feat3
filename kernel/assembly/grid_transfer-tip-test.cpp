#include <test_system/test_system.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename ShapeType, template<typename> class Element_, Index level_coarse_>
class GridTransferTipTest :
  public TestSystem::BaseTest
{
  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;

  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Element_<TrafoType> SpaceType;

public:
  explicit GridTransferTipTest() :
    TestSystem::BaseTest("GridTransferTipTest<" + ShapeType::name() + "," + SpaceType::name() + ">")
  {
  }

  virtual ~GridTransferTipTest()
  {
  }

  virtual void run() const override
  {
    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<MeshType> coarse_factory(level_coarse_);
    MeshType mesh_c(coarse_factory);

    // refine the mesh
    Geometry::StandardRefinery<MeshType> refine_factory(mesh_c);
    MeshType mesh_f(refine_factory);

    // compute eps
    const DataType eps = Math::pow(Math::eps<DataType>(), DataType(0.8));

    // create trafos
    TrafoType trafo_f(mesh_f);
    TrafoType trafo_c(mesh_c);

    // create spaces
    SpaceType space_f(trafo_f);
    SpaceType space_c(trafo_c);

    // create a cubature factory of appropriate degree
    Cubature::DynamicFactory cubature_factory("auto-degree:" + stringify(Math::sqr(SpaceType::local_degree+1)+2));

    // assemble prolongation matrix "P"
    MatrixType prol_matrix;
    {
      Assembly::SymbolicAssembler::assemble_matrix_2lvl(prol_matrix, space_f, space_c);
      VectorType weight_vector(prol_matrix.create_vector_l());
      prol_matrix.format();
      weight_vector.format();
      Assembly::GridTransfer::assemble_prolongation(prol_matrix, weight_vector, space_f, space_c, cubature_factory);
      weight_vector.component_invert(weight_vector);
      prol_matrix.scale_rows(prol_matrix, weight_vector);
    }

    // assemble truncation matrix "T"
    MatrixType trunc_matrix(prol_matrix.transpose());
    {
      VectorType weight_vector(trunc_matrix.create_vector_l());
      trunc_matrix.format();
      weight_vector.format();
      Assembly::GridTransfer::assemble_truncation(trunc_matrix, weight_vector, space_f, space_c, cubature_factory);
      weight_vector.component_invert(weight_vector);
      trunc_matrix.scale_rows(trunc_matrix, weight_vector);
    }

    // assemble extended matrix structure on coarse mesh
    MatrixType rip_matrix;
    Assembly::SymbolicAssembler::assemble_matrix_ext1(rip_matrix, space_c);

    // initialise to coarse-mesh identity matrix "I_c"
    {
      const IndexType num_rows = IndexType(rip_matrix.rows());
      const IndexType* row_ptr = rip_matrix.row_ptr();
      const IndexType* col_idx = rip_matrix.col_ind();
      DataType* vals = rip_matrix.val();
      for(IndexType i(0); i < num_rows; ++i)
        for(IndexType j(row_ptr[i]); j < row_ptr[i+1]; ++j)
          vals[j] = DataType(col_idx[j] == i ? 1 : 0);
    }

    // compute test matrix: X := I_c - T * I_f * P
    VectorType vec_id(space_f.get_num_dofs(), DataType(1));
    rip_matrix.add_double_mat_mult(trunc_matrix, vec_id, prol_matrix, -DataType(1), true);

    // the resulting matrix should now be the null matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqr(rip_matrix.norm_frobenius()), DataType(0), eps);
  }
}; // GridTransferTipTest<...>

// Lagrange-1 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4> grid_transfer_truncate_test_hy1_lagrange1;
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2> grid_transfer_truncate_test_hy2_lagrange1;
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1> grid_transfer_truncate_test_hy3_lagrange1;
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2> grid_transfer_truncate_test_sx2_lagrange1;
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1> grid_transfer_truncate_test_sx3_lagrange1;

// Lagrange-2 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4> grid_transfer_truncate_test_hy1_lagrange2;
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2> grid_transfer_truncate_test_hy2_lagrange2;
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1> grid_transfer_truncate_test_hy3_lagrange2;
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2> grid_transfer_truncate_test_sx2_lagrange2;
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1> grid_transfer_truncate_test_sx3_lagrange2;

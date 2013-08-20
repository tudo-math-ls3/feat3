#include <test_system/test_system.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename DataType_>
class GridTransferTest :
  public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_> MatrixType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
  typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceQ1T;

public:
  GridTransferTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("GridTransferTest")
  {
  }

  virtual void run() const
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    // create coarse mesh
    Geometry::UnitCubeFactory<QuadMesh> unit_factory;
    QuadMesh mesh_coarse(unit_factory);

    // refine the mesh
    Geometry::StandardRefinery<QuadMesh> refine_factory(mesh_coarse);
    QuadMesh mesh_fine(refine_factory);

    // run tests
    test_unit_2d_q1(mesh_fine, mesh_coarse);
    test_unit_2d_q1t(mesh_fine, mesh_coarse);
  }

  void test_unit_2d_q1(const QuadMesh& mesh_f, const QuadMesh& mesh_c) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

    // create trafos
    QuadTrafo trafo_f(mesh_f);
    QuadTrafo trafo_c(mesh_c);

    // create spaces
    QuadSpaceQ1 space_f(trafo_f);
    QuadSpaceQ1 space_c(trafo_c);

    // assemble matrix structure
    typedef Space::DofAdjacency<Space::Stencil::StandardRefinement> DofAdjacency;
    MatrixType prol_matrix(DofAdjacency::assemble(space_f, space_c));

    // check matrix dimensions
    TEST_CHECK_EQUAL(prol_matrix.rows(), 9u);
    TEST_CHECK_EQUAL(prol_matrix.columns(), 4u);
    TEST_CHECK_EQUAL(prol_matrix.used_elements(), 36u);

    // fetch matrix arrays
    const Index* row_ptr = prol_matrix.row_ptr();
    const Index* col_idx = prol_matrix.col_ind();

    // check matrix structure; has to be a dense 9x4 matrix
    for(Index i(0); i < 9; ++i)
    {
      TEST_CHECK_EQUAL(row_ptr[i], 4*i);
      for(Index j(0); j < 4; ++j)
      {
        TEST_CHECK_EQUAL(col_idx[4*i+j], j);
      }
    }
    TEST_CHECK_EQUAL(row_ptr[9], 36u);

    // clear matrix
    prol_matrix.clear(DataType_(0));

    // assemble prolongation matrix
    Assembly::GridTransfer::
      assemble_prolongation(prol_matrix, Cubature::DynamicFactory("gauss-legendre:2"), space_f, space_c);

    // fetch matrix data
    const DataType_* data = prol_matrix.val();

    // reference data array
    static const DataType_ d0 = DataType_(0);
    static const DataType_ d1 = DataType_(1);
    static const DataType_ d2 = DataType_(0.5);
    static const DataType_ d4 = DataType_(0.25);
    static const DataType_ data_ref[] =
    {
      d1, d0, d0, d0,
      d0, d1, d0, d0,
      d0, d0, d1, d0,
      d0, d0, d0, d1,
      d2, d2, d0, d0,
      d0, d0, d2, d2,
      d2, d0, d2, d0,
      d0, d2, d0, d2,
      d4, d4, d4, d4
    };
    for(Index i(0); i < 36; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(data[i], data_ref[i], eps);
    }
  }

  void test_unit_2d_q1t(const QuadMesh& mesh_f, const QuadMesh& mesh_c) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.8));

    // create trafos
    QuadTrafo trafo_f(mesh_f);
    QuadTrafo trafo_c(mesh_c);

    // create spaces
    QuadSpaceQ1T space_f(trafo_f);
    QuadSpaceQ1T space_c(trafo_c);

    // assemble matrix structure
    typedef Space::DofAdjacency<Space::Stencil::StandardRefinement> DofAdjacency;
    MatrixType prol_matrix(DofAdjacency::assemble(space_f, space_c));

    // check matrix dimensions
    TEST_CHECK_EQUAL(prol_matrix.rows(), 12u);
    TEST_CHECK_EQUAL(prol_matrix.columns(), 4u);
    TEST_CHECK_EQUAL(prol_matrix.used_elements(), 48u);

    // fetch matrix arrays
    const Index* row_ptr = prol_matrix.row_ptr();
    const Index* col_idx = prol_matrix.col_ind();

    // check matrix structure; has to be a dense 9x4 matrix
    for(Index i(0); i < 12; ++i)
    {
      TEST_CHECK_EQUAL(row_ptr[i], 4*i);
      for(Index j(0); j < 4; ++j)
      {
        TEST_CHECK_EQUAL(col_idx[4*i+j], j);
      }
    }
    TEST_CHECK_EQUAL(row_ptr[12], 48u);

    // clear matrix
    prol_matrix.clear(DataType_(0));

    // assemble prolongation matrix
    Assembly::GridTransfer::
      assemble_prolongation(prol_matrix, Cubature::DynamicFactory("gauss-legendre:3"), space_f, space_c);

    // fetch matrix data
    const DataType_* data = prol_matrix.val();

    // reference data array
    static const DataType_ d0 = DataType_(0);
    static const DataType_ d1 = DataType_(1);
    static const DataType_ d4 = DataType_(0.25);   // = 1/4
    static const DataType_ d18 = DataType_(0.125); // = 1/8
    static const DataType_ d58 = DataType_(0.625); // = 5/8
    static const DataType_ data_ref[] =
    {
      d1, d0,  d4, -d4,
      d1, d0, -d4,  d4,
      d0, d1,  d4, -d4,
      d0, d1, -d4,  d4,
       d4, -d4, d1, d0,
      -d4,  d4, d1, d0,
       d4, -d4, d0, d1,
      -d4,  d4, d0, d1,
      d58, d18, d18, d18,
      d18, d58, d18, d18,
      d18, d18, d58, d18,
      d18, d18, d18, d58
    };
    for(Index i(0); i < 48; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(data[i], data_ref[i], eps);
    }
  }
};

GridTransferTest<float> grid_transfer_test_float;
GridTransferTest<double> grid_transfer_test_double;

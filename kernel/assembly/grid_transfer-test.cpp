#include <test_system/test_system.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename MatrixType_>
class GridTransferTest :
  public TestSystem::FullTaggedTest<Archs::None, typename MatrixType_::DataType, typename MatrixType_::IndexType>
{
  typedef typename MatrixType_::MemType MemType_;
  typedef typename MatrixType_::DataType DataType_;
  typedef typename MatrixType_::IndexType IndexType_;
  typedef LAFEM::DenseVector<MemType_, DataType_, IndexType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
  typedef Space::CroRavRanTur::Element<QuadTrafo> QuadSpaceQ1T;

public:
  GridTransferTest() :
    TestSystem::FullTaggedTest<Archs::None, DataType_, IndexType_>("GridTransferTest<" + MatrixType_::name() + ">")
  {
  }

  virtual ~GridTransferTest()
  {
  }

  virtual void run() const override
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

  void test_unit_2d_q1(QuadMesh& mesh_f, QuadMesh& mesh_c) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create trafos
    QuadTrafo trafo_f(mesh_f);
    QuadTrafo trafo_c(mesh_c);

    // create spaces
    QuadSpaceQ1 space_f(trafo_f);
    QuadSpaceQ1 space_c(trafo_c);

    // assemble matrix structure
    MatrixType_ prol_matrix;
    VectorType weight_vector(space_f.get_num_dofs());
    Assembly::SymbolicAssembler::assemble_matrix_2lvl(prol_matrix, space_f, space_c);

    // check matrix dimensions
    TEST_CHECK_EQUAL(prol_matrix.rows(), 9u);
    TEST_CHECK_EQUAL(prol_matrix.columns(), 4u);
    TEST_CHECK_EQUAL(prol_matrix.used_elements(), 36u);

    // fetch matrix arrays of matrix in CSR-format
    LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> tmp_prol_matrix;
    tmp_prol_matrix.convert(prol_matrix);
    const IndexType_ * row_ptr = tmp_prol_matrix.row_ptr();
    const IndexType_ * col_idx = tmp_prol_matrix.col_ind();

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
    prol_matrix.format();
    weight_vector.format();

    // assemble prolongation matrix
    Assembly::GridTransfer::
      assemble_prolongation(prol_matrix, weight_vector, space_f, space_c, Cubature::DynamicFactory("gauss-legendre:2"));

    // invert weight vector and scale matrix rows
    weight_vector.component_invert(weight_vector);
    prol_matrix.scale_rows(prol_matrix, weight_vector);

    // reference data array
    const DataType_ d0 = DataType_(0);
    const DataType_ d1 = DataType_(1);
    const DataType_ d2 = DataType_(0.5);
    const DataType_ d4 = DataType_(0.25);
    const DataType_ data_ref[] =
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

    for (Index i(0); i < prol_matrix.rows(); ++i)
    {
      for (Index j(0); j < prol_matrix.columns(); ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(prol_matrix(i,j), data_ref[i * prol_matrix.columns() + j], eps);
      }
    }
  }

  void test_unit_2d_q1t(QuadMesh& mesh_f, QuadMesh& mesh_c) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create trafos
    QuadTrafo trafo_f(mesh_f);
    QuadTrafo trafo_c(mesh_c);

    // create spaces
    QuadSpaceQ1T space_f(trafo_f);
    QuadSpaceQ1T space_c(trafo_c);

    // assemble matrix structure
    MatrixType_ prol_matrix;
    VectorType weight_vector(space_f.get_num_dofs());
    Assembly::SymbolicAssembler::assemble_matrix_2lvl(prol_matrix, space_f, space_c);

    // check matrix dimensions
    TEST_CHECK_EQUAL(prol_matrix.rows(), 12u);
    TEST_CHECK_EQUAL(prol_matrix.columns(), 4u);
    TEST_CHECK_EQUAL(prol_matrix.used_elements(), 48u);

    // fetch matrix arrays of matrix in CSR-format
    LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> tmp_prol_matrix;
    tmp_prol_matrix.convert(prol_matrix);
    const IndexType_ * row_ptr = tmp_prol_matrix.row_ptr();
    const IndexType_ * col_idx = tmp_prol_matrix.col_ind();

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
    prol_matrix.format();
    weight_vector.format();

    // assemble prolongation matrix
    Assembly::GridTransfer::
      assemble_prolongation(prol_matrix, weight_vector, space_f, space_c, Cubature::DynamicFactory("gauss-legendre:3"));

    // invert weight vector and scale matrix rows
    weight_vector.component_invert(weight_vector);
    prol_matrix.scale_rows(prol_matrix, weight_vector);

    // reference data array
    const DataType_ d0 = DataType_(0);
    const DataType_ d1 = DataType_(1);
    const DataType_ d4 = DataType_(0.25);   // = 1/4
    const DataType_ d18 = DataType_(0.125); // = 1/8
    const DataType_ d58 = DataType_(0.625); // = 5/8
    const DataType_ data_ref[] =
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
    for (Index i(0); i < prol_matrix.rows(); ++i)
    {
      for (Index j(0); j < prol_matrix.columns(); ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(prol_matrix(i,j), data_ref[i * prol_matrix.columns() + j], eps);
      }
    }
  }
};

GridTransferTest<LAFEM::SparseMatrixCSR<Mem::Main, float, unsigned int> > grid_transfer_test_csr_float_uint;
GridTransferTest<LAFEM::SparseMatrixCSR<Mem::Main, float, unsigned long> > grid_transfer_test_csr_float_ulong;
GridTransferTest<LAFEM::SparseMatrixCSR<Mem::Main, double, unsigned int> > grid_transfer_test_csr_double_uint;
GridTransferTest<LAFEM::SparseMatrixCSR<Mem::Main, double, unsigned long> > grid_transfer_test_csr_double_ulong;

GridTransferTest<LAFEM::SparseMatrixCOO<Mem::Main, float, unsigned int> > grid_transfer_test_coo_float_uint;
GridTransferTest<LAFEM::SparseMatrixCOO<Mem::Main, float, unsigned long> > grid_transfer_test_coo_float_ulong;
GridTransferTest<LAFEM::SparseMatrixCOO<Mem::Main, double, unsigned int> > grid_transfer_test_coo_double_uint;
GridTransferTest<LAFEM::SparseMatrixCOO<Mem::Main, double, unsigned long> > grid_transfer_test_coo_double_ulong;

GridTransferTest<LAFEM::SparseMatrixELL<Mem::Main, float, unsigned int> > grid_transfer_test_ell_float_uint;
GridTransferTest<LAFEM::SparseMatrixELL<Mem::Main, float, unsigned long> > grid_transfer_test_ell_float_ulong;
GridTransferTest<LAFEM::SparseMatrixELL<Mem::Main, double, unsigned int> > grid_transfer_test_ell_double_uint;
GridTransferTest<LAFEM::SparseMatrixELL<Mem::Main, double, unsigned long> > grid_transfer_test_ell_double_ulong;

#ifdef FEAT_HAVE_QUADMATH
GridTransferTest<LAFEM::SparseMatrixCSR<Mem::Main, __float128, unsigned int> > grid_transfer_test_csr_float128_uint;
GridTransferTest<LAFEM::SparseMatrixCSR<Mem::Main, __float128, unsigned long> > grid_transfer_test_csr_float128_ulong;

GridTransferTest<LAFEM::SparseMatrixCOO<Mem::Main, __float128, unsigned int> > grid_transfer_test_coo_float128_uint;
GridTransferTest<LAFEM::SparseMatrixCOO<Mem::Main, __float128, unsigned long> > grid_transfer_test_coo_float128_ulong;

GridTransferTest<LAFEM::SparseMatrixELL<Mem::Main, __float128, unsigned int> > grid_transfer_test_ell_float128_uint;
GridTransferTest<LAFEM::SparseMatrixELL<Mem::Main, __float128, unsigned long> > grid_transfer_test_ell_float128_ulong;
#endif

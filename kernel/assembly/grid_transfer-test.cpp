// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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

template<typename DT_, typename IT_>
class GridTransferTest :
  public UnitTest
{
  typedef DT_ DataType_;
  typedef IT_ IndexType_;
  typedef LAFEM::SparseMatrixCSR<DataType_, IndexType_> MatrixType_;
  typedef LAFEM::DenseVector<DataType_, IndexType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
  typedef Space::CroRavRanTur::Element<QuadTrafo> QuadSpaceQ1T;

public:
  GridTransferTest(PreferredBackend backend) :
    UnitTest("GridTransferTest<" + MatrixType_::name() + ">", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
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
    LAFEM::SparseMatrixCSR<DataType_, IndexType_> tmp_prol_matrix;
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
    LAFEM::SparseMatrixCSR<DataType_, IndexType_> tmp_prol_matrix;
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

GridTransferTest <float, std::uint32_t> grid_transfer_test_csr_float_uint32(PreferredBackend::generic);
GridTransferTest <float, std::uint64_t> grid_transfer_test_csr_float_uint64(PreferredBackend::generic);
GridTransferTest <double, std::uint32_t> grid_transfer_test_csr_double_uint32(PreferredBackend::generic);
GridTransferTest <double, std::uint64_t> grid_transfer_test_csr_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
GridTransferTest <float, std::uint64_t> mkl_grid_transfer_test_csr_float_uint64(PreferredBackend::mkl);
GridTransferTest <double, std::uint64_t> mkl_grid_transfer_test_csr_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
GridTransferTest <__float128, std::uint32_t> grid_transfer_test_csr_float128_uint32(PreferredBackend::generic);
GridTransferTest <__float128, std::uint64_t> grid_transfer_test_csr_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
GridTransferTest <Half, std::uint32_t> grid_transfer_test_csr_half_uint32(PreferredBackend::generic);
GridTransferTest <Half, std::uint64_t> grid_transfer_test_csr_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
GridTransferTest <float, std::uint32_t> cuda_grid_transfer_test_csr_float_uint32(PreferredBackend::cuda);
GridTransferTest <double, std::uint32_t> cuda_grid_transfer_test_csr_double_uint32(PreferredBackend::cuda);
GridTransferTest <float, std::uint64_t> cuda_grid_transfer_test_csr_float_uint64(PreferredBackend::cuda);
GridTransferTest <double, std::uint64_t> cuda_grid_transfer_test_csr_double_uint64(PreferredBackend::cuda);
#endif

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/space_transfer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename DT_, typename IT_>
class SpaceTransferTest :
  public UnitTest
{
  typedef DT_ DataType_;
  typedef IT_ IndexType_;
  typedef LAFEM::SparseMatrixCSR<DataType_, IndexType_> MatrixType;
  typedef LAFEM::DenseVector<DataType_, IndexType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Lagrange2::Element<QuadTrafo> QuadSpaceQ2;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  SpaceTransferTest(PreferredBackend backend) :
    UnitTest("SpaceTransferTest<" + MatrixType::name() + ">", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~SpaceTransferTest()
  {
  }

  virtual void run() const override
  {
    test_unit_2d_q2q1();
  }

  void test_unit_2d_q2q1() const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.75));

    String cubature = "gauss-legendre:3";

    // create mesh, trafo and spaces
    Geometry::RefinedUnitCubeFactory<QuadMesh> mesh_factory(2);
    QuadMesh mesh(mesh_factory);
    QuadTrafo trafo(mesh);
    QuadSpaceQ2 space_q2(trafo);
    QuadSpaceQ1 space_q1(trafo);

    // assemble matrix structure
    MatrixType matrix_q2_to_q1;
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_q2_to_q1, space_q1, space_q2);
    MatrixType matrix_q1_to_q2 = matrix_q2_to_q1.transpose();

    // assemble transfer matrices
    matrix_q1_to_q2.format();
    matrix_q2_to_q1.format();
    Assembly::SpaceTransfer::assemble_transfer_direct(matrix_q1_to_q2, space_q2, space_q1, cubature);
    Assembly::SpaceTransfer::assemble_transfer_direct(matrix_q2_to_q1, space_q1, space_q2, cubature);

    // interpolate a Q1 function
    VectorType vec_q2, vec_q1;
    auto func_q1 = Analytic::create_lambda_function_scalar_2d(
      [] (DataType_ x, DataType_ y) {return x*y + DataType_(2)*x - DataType_(3)*y;});
    Assembly::Interpolator::project(vec_q1, func_q1, space_q1);
    Assembly::Interpolator::project(vec_q2, func_q1, space_q2);

    // project Q1 vector into Q2 space with matrix
    VectorType vec_q2t_1 = matrix_q1_to_q2.create_vector_l();
    matrix_q1_to_q2.apply(vec_q2t_1, vec_q1);

    // project Q1 vector into Q2 space directly
    VectorType vec_q2t_2 = matrix_q1_to_q2.create_vector_l();
    vec_q2t_2.format();
    Assembly::SpaceTransfer::transfer_vector_direct(vec_q2t_2, vec_q1, space_q2, space_q1, cubature);

    // compare vectors
    vec_q2t_1.axpy(vec_q2, -DataType_(1));
    TEST_CHECK(vec_q2t_1.norm2() < eps);

    vec_q2t_2.axpy(vec_q2, -DataType_(1));
    TEST_CHECK(vec_q2t_2.norm2() < eps);

    // project Q2 vector into Q1 space with matrix
    VectorType vec_q1t_1 = matrix_q2_to_q1.create_vector_l();
    matrix_q2_to_q1.apply(vec_q1t_1, vec_q2);

    // project Q2 vector into Q1 space directly
    VectorType vec_q1t_2 = matrix_q2_to_q1.create_vector_l();
    vec_q1t_2.format();
    Assembly::SpaceTransfer::transfer_vector_direct(vec_q1t_2, vec_q2, space_q1, space_q2, cubature);

    // compare vectors
    vec_q1t_1.axpy(vec_q1, -DataType_(1));
    TEST_CHECK(vec_q1t_1.norm2() < eps);

    vec_q1t_2.axpy(vec_q1, -DataType_(1));
    TEST_CHECK(vec_q1t_2.norm2() < eps);
  }
};

SpaceTransferTest <double, std::uint64_t> grid_transfer_test_csr_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SpaceTransferTest <double, std::uint64_t> mkl_grid_transfer_test_csr_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SpaceTransferTest <__float128, std::uint64_t> grid_transfer_test_csr_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SpaceTransferTest <Half, std::uint64_t> grid_transfer_test_csr_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SpaceTransferTest <float, std::uint32_t> cuda_grid_transfer_test_csr_float_uint32(PreferredBackend::cuda);
SpaceTransferTest <double, std::uint64_t> cuda_grid_transfer_test_csr_double_uint64(PreferredBackend::cuda);
#endif

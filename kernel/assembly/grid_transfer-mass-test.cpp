// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename ShapeType, template<typename> class Element_, Index level_coarse_, typename DT_, typename IT_>
class GridTransferMassTest :
  public UnitTest
{
  typedef DT_ DataType;
  typedef IT_ IndexType;

  typedef LAFEM::DenseVector<DataType, IndexType> VectorType;
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;

  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Element_<TrafoType> SpaceType;

public:
  explicit GridTransferMassTest(PreferredBackend backend) :
    UnitTest("GridTransferMassTest<" + ShapeType::name() + "," + SpaceType::name() + ">", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~GridTransferMassTest()
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

    // assemble fine/coarse mesh mass matrix
    MatrixType mass_f, mass_c;
    Assembly::SymbolicAssembler::assemble_matrix_std1(mass_f, space_f);
    Assembly::SymbolicAssembler::assemble_matrix_std1(mass_c, space_c);
    mass_f.format();
    mass_c.format();
    Assembly::Common::IdentityOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(mass_f, operat, space_f, cubature_factory);
    Assembly::BilinearOperatorAssembler::assemble_matrix1(mass_c, operat, space_c, cubature_factory);

    // assemble prolongation matrix
    MatrixType prol_matrix;
    VectorType weight_vector(space_f.get_num_dofs());
    Assembly::SymbolicAssembler::assemble_matrix_2lvl(prol_matrix, space_f, space_c);
    prol_matrix.format();
    weight_vector.format();
    Assembly::GridTransfer::assemble_prolongation(prol_matrix, weight_vector, space_f, space_c, cubature_factory);
    weight_vector.component_invert(weight_vector);
    prol_matrix.scale_rows(prol_matrix, weight_vector);

    // transpose to obtain restriction matrix
    MatrixType rest_matrix(prol_matrix.transpose());

    // finally, restrict the fine mesh mass matrix onto the coarse mesh and
    // subtract it from the coarse mesh mass matrix, i.e.
    // M_c <- M_c - R * M_f * P
    mass_c.add_double_mat_product(rest_matrix, mass_f, prol_matrix, -DataType(1), true);

    // the resulting matrix should now be the null matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqr(mass_c.norm_frobenius()), DataType(0), eps);
  }
};

// Lagrange-1 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned int> grid_transfer_mass_test_hy1_lagrange1_double_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned int> grid_transfer_mass_test_hy2_lagrange1_float_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned int> grid_transfer_mass_test_hy3_lagrange1_double_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned int> grid_transfer_mass_test_sx2_lagrange1_double_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned int> grid_transfer_mass_test_sx3_lagrange1_float_uint(PreferredBackend::generic);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned long> grid_transfer_mass_test_hy1_lagrange1_double_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned long> grid_transfer_mass_test_hy2_lagrange1_float_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned long> grid_transfer_mass_test_hy3_lagrange1_double_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned long> grid_transfer_mass_test_sx2_lagrange1_double_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned long> grid_transfer_mass_test_sx3_lagrange1_float_ulong(PreferredBackend::generic);

// Lagrange-2 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned int> grid_transfer_mass_test_hy1_lagrange2_double_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned int> grid_transfer_mass_test_hy2_lagrange2_double_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned int> grid_transfer_mass_test_hy3_lagrange2_float_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned int> grid_transfer_mass_test_sx2_lagrange2_float_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned int> grid_transfer_mass_test_sx3_lagrange2_double_uint(PreferredBackend::generic);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned long> grid_transfer_mass_test_hy1_lagrange2_double_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned long> grid_transfer_mass_test_hy2_lagrange2_double_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned long> grid_transfer_mass_test_hy3_lagrange2_float_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned long> grid_transfer_mass_test_sx2_lagrange2_float_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned long> grid_transfer_mass_test_sx3_lagrange2_double_ulong(PreferredBackend::generic);

#ifdef FEAT_HAVE_MKL
// Lagrange-1 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned long> mkl_grid_transfer_mass_test_hy1_lagrange1_double_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned long> mkl_grid_transfer_mass_test_hy2_lagrange1_float_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned long> mkl_grid_transfer_mass_test_hy3_lagrange1_double_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned long> mkl_grid_transfer_mass_test_sx2_lagrange1_double_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned long> mkl_grid_transfer_mass_test_sx3_lagrange1_float_ulong(PreferredBackend::mkl);
// Lagrange-2 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned long> mkl_grid_transfer_mass_test_hy1_lagrange2_double_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned long> mkl_grid_transfer_mass_test_hy2_lagrange2_double_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned long> mkl_grid_transfer_mass_test_hy3_lagrange2_float_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned long> mkl_grid_transfer_mass_test_sx2_lagrange2_float_ulong(PreferredBackend::mkl);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned long> mkl_grid_transfer_mass_test_sx3_lagrange2_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
// Lagrange-1 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, __float128, unsigned int> grid_transfer_mass_test_hy1_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, __float128, unsigned int> grid_transfer_mass_test_hy2_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, __float128, unsigned int> grid_transfer_mass_test_hy3_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, __float128, unsigned int> grid_transfer_mass_test_sx2_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, __float128, unsigned int> grid_transfer_mass_test_sx3_lagrange1_float128_uint(PreferredBackend::generic);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, __float128, unsigned long> grid_transfer_mass_test_hy1_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, __float128, unsigned long> grid_transfer_mass_test_hy2_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, __float128, unsigned long> grid_transfer_mass_test_hy3_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, __float128, unsigned long> grid_transfer_mass_test_sx2_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, __float128, unsigned long> grid_transfer_mass_test_sx3_lagrange1_float128_ulong(PreferredBackend::generic);
// Lagrange-2 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, __float128, unsigned int> grid_transfer_mass_test_hy1_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, __float128, unsigned int> grid_transfer_mass_test_hy2_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, __float128, unsigned int> grid_transfer_mass_test_hy3_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, __float128, unsigned int> grid_transfer_mass_test_sx2_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, __float128, unsigned int> grid_transfer_mass_test_sx3_lagrange2_float128_uint(PreferredBackend::generic);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, __float128, unsigned long> grid_transfer_mass_test_hy1_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, __float128, unsigned long> grid_transfer_mass_test_hy2_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, __float128, unsigned long> grid_transfer_mass_test_hy3_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, __float128, unsigned long> grid_transfer_mass_test_sx2_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, __float128, unsigned long> grid_transfer_mass_test_sx3_lagrange2_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
// Lagrange-1 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, Half, unsigned int> grid_transfer_mass_test_hy1_lagrange1_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, Half, unsigned int> grid_transfer_mass_test_hy2_lagrange1_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, Half, unsigned int> grid_transfer_mass_test_hy3_lagrange1_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, Half, unsigned int> grid_transfer_mass_test_sx2_lagrange1_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, Half, unsigned int> grid_transfer_mass_test_sx3_lagrange1_half_uint(PreferredBackend::generic);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, Half, unsigned long> grid_transfer_mass_test_hy1_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, Half, unsigned long> grid_transfer_mass_test_hy2_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, Half, unsigned long> grid_transfer_mass_test_hy3_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, Half, unsigned long> grid_transfer_mass_test_sx2_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, Half, unsigned long> grid_transfer_mass_test_sx3_lagrange1_half_ulong(PreferredBackend::generic);
// Lagrange-2 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, Half, unsigned int> grid_transfer_mass_test_hy1_lagrange2_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, Half, unsigned int> grid_transfer_mass_test_hy2_lagrange2_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, Half, unsigned int> grid_transfer_mass_test_hy3_lagrange2_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, Half, unsigned int> grid_transfer_mass_test_sx2_lagrange2_half_uint(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, Half, unsigned int> grid_transfer_mass_test_sx3_lagrange2_half_uint(PreferredBackend::generic);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, Half, unsigned long> grid_transfer_mass_test_hy1_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, Half, unsigned long> grid_transfer_mass_test_hy2_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, Half, unsigned long> grid_transfer_mass_test_hy3_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, Half, unsigned long> grid_transfer_mass_test_sx2_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, Half, unsigned long> grid_transfer_mass_test_sx3_lagrange2_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
// Lagrange-1 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned int> cuda_grid_transfer_mass_test_hy1_lagrange1_double_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned int> cuda_grid_transfer_mass_test_hy2_lagrange1_float_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned int> cuda_grid_transfer_mass_test_hy3_lagrange1_double_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned int> cuda_grid_transfer_mass_test_sx2_lagrange1_double_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned int> cuda_grid_transfer_mass_test_sx3_lagrange1_float_uint(PreferredBackend::cuda);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned long> cuda_grid_transfer_mass_test_hy1_lagrange1_double_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned long> cuda_grid_transfer_mass_test_hy2_lagrange1_float_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned long> cuda_grid_transfer_mass_test_hy3_lagrange1_double_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned long> cuda_grid_transfer_mass_test_sx2_lagrange1_double_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned long> cuda_grid_transfer_mass_test_sx3_lagrange1_float_ulong(PreferredBackend::cuda);
// Lagrange-2 element
GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned int> cuda_grid_transfer_mass_test_hy1_lagrange2_double_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned int> cuda_grid_transfer_mass_test_hy2_lagrange2_double_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned int> cuda_grid_transfer_mass_test_hy3_lagrange2_float_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned int> cuda_grid_transfer_mass_test_sx2_lagrange2_float_uint(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned int> cuda_grid_transfer_mass_test_sx3_lagrange2_double_uint(PreferredBackend::cuda);

GridTransferMassTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned long> cuda_grid_transfer_mass_test_hy1_lagrange2_double_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned long> cuda_grid_transfer_mass_test_hy2_lagrange2_double_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned long> cuda_grid_transfer_mass_test_hy3_lagrange2_float_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned long> cuda_grid_transfer_mass_test_sx2_lagrange2_float_ulong(PreferredBackend::cuda);
GridTransferMassTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned long> cuda_grid_transfer_mass_test_sx3_lagrange2_double_ulong(PreferredBackend::cuda);
#endif

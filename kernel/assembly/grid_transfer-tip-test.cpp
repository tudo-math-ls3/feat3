// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
class GridTransferTipTest :
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
  explicit GridTransferTipTest(PreferredBackend backend) :
    UnitTest("GridTransferTipTest<" + ShapeType::name() + "," + SpaceType::name() + ">", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
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
    Assembly::SymbolicAssembler::assemble_matrix_ext_facet1(rip_matrix, space_c);

    // initialize to coarse-mesh identity matrix "I_c"
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
    rip_matrix.add_double_mat_product(trunc_matrix, vec_id, prol_matrix, -DataType(1), true);

    // the resulting matrix should now be the null matrix
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqr(rip_matrix.norm_frobenius()), DataType(0), eps);
  }
}; // GridTransferTipTest<...>

// Lagrange-1 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned int> grid_transfer_truncate_test_hy1_lagrange1_double_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned int> grid_transfer_truncate_test_hy2_lagrange1_float_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned int> grid_transfer_truncate_test_hy3_lagrange1_double_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned int> grid_transfer_truncate_test_sx2_lagrange1_double_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned int> grid_transfer_truncate_test_sx3_lagrange1_float_uint(PreferredBackend::generic);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned long> grid_transfer_truncate_test_hy1_lagrange1_double_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned long> grid_transfer_truncate_test_hy2_lagrange1_float_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned long> grid_transfer_truncate_test_hy3_lagrange1_double_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned long> grid_transfer_truncate_test_sx2_lagrange1_double_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned long> grid_transfer_truncate_test_sx3_lagrange1_float_ulong(PreferredBackend::generic);

// Lagrange-2 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned int> grid_transfer_truncate_test_hy1_lagrange2_double_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned int> grid_transfer_truncate_test_hy2_lagrange2_double_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned int> grid_transfer_truncate_test_hy3_lagrange2_float_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned int> grid_transfer_truncate_test_sx2_lagrange2_float_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned int> grid_transfer_truncate_test_sx3_lagrange2_double_uint(PreferredBackend::generic);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned long> grid_transfer_truncate_test_hy1_lagrange2_double_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned long> grid_transfer_truncate_test_hy2_lagrange2_double_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned long> grid_transfer_truncate_test_hy3_lagrange2_float_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned long> grid_transfer_truncate_test_sx2_lagrange2_float_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned long> grid_transfer_truncate_test_sx3_lagrange2_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
// Lagrange-1 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned long> mkl_grid_transfer_truncate_test_hy1_lagrange1_double_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned long> mkl_grid_transfer_truncate_test_hy2_lagrange1_float_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned long> mkl_grid_transfer_truncate_test_hy3_lagrange1_double_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned long> mkl_grid_transfer_truncate_test_sx2_lagrange1_double_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned long> mkl_grid_transfer_truncate_test_sx3_lagrange1_float_ulong(PreferredBackend::mkl);
// Lagrange-2 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned long> mkl_grid_transfer_truncate_test_hy1_lagrange2_double_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned long> mkl_grid_transfer_truncate_test_hy2_lagrange2_double_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned long> mkl_grid_transfer_truncate_test_hy3_lagrange2_float_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned long> mkl_grid_transfer_truncate_test_sx2_lagrange2_float_ulong(PreferredBackend::mkl);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned long> mkl_grid_transfer_truncate_test_sx3_lagrange2_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
// Lagrange-1 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, __float128, unsigned int> grid_transfer_truncate_test_hy1_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, __float128, unsigned int> grid_transfer_truncate_test_hy2_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, __float128, unsigned int> grid_transfer_truncate_test_hy3_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, __float128, unsigned int> grid_transfer_truncate_test_sx2_lagrange1_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, __float128, unsigned int> grid_transfer_truncate_test_sx3_lagrange1_float128_uint(PreferredBackend::generic);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, __float128, unsigned long> grid_transfer_truncate_test_hy1_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, __float128, unsigned long> grid_transfer_truncate_test_hy2_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, __float128, unsigned long> grid_transfer_truncate_test_hy3_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, __float128, unsigned long> grid_transfer_truncate_test_sx2_lagrange1_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, __float128, unsigned long> grid_transfer_truncate_test_sx3_lagrange1_float128_ulong(PreferredBackend::generic);
// Lagrange-2 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, __float128, unsigned int> grid_transfer_truncate_test_hy1_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, __float128, unsigned int> grid_transfer_truncate_test_hy2_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, __float128, unsigned int> grid_transfer_truncate_test_hy3_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, __float128, unsigned int> grid_transfer_truncate_test_sx2_lagrange2_float128_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, __float128, unsigned int> grid_transfer_truncate_test_sx3_lagrange2_float128_uint(PreferredBackend::generic);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, __float128, unsigned long> grid_transfer_truncate_test_hy1_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, __float128, unsigned long> grid_transfer_truncate_test_hy2_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, __float128, unsigned long> grid_transfer_truncate_test_hy3_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, __float128, unsigned long> grid_transfer_truncate_test_sx2_lagrange2_float128_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, __float128, unsigned long> grid_transfer_truncate_test_sx3_lagrange2_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
// Lagrange-1 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, Half, unsigned int> grid_transfer_truncate_test_hy1_lagrange1_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, Half, unsigned int> grid_transfer_truncate_test_hy2_lagrange1_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, Half, unsigned int> grid_transfer_truncate_test_hy3_lagrange1_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, Half, unsigned int> grid_transfer_truncate_test_sx2_lagrange1_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, Half, unsigned int> grid_transfer_truncate_test_sx3_lagrange1_half_uint(PreferredBackend::generic);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, Half, unsigned long> grid_transfer_truncate_test_hy1_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, Half, unsigned long> grid_transfer_truncate_test_hy2_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, Half, unsigned long> grid_transfer_truncate_test_hy3_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, Half, unsigned long> grid_transfer_truncate_test_sx2_lagrange1_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, Half, unsigned long> grid_transfer_truncate_test_sx3_lagrange1_half_ulong(PreferredBackend::generic);
// Lagrange-2 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, Half, unsigned int> grid_transfer_truncate_test_hy1_lagrange2_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, Half, unsigned int> grid_transfer_truncate_test_hy2_lagrange2_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, Half, unsigned int> grid_transfer_truncate_test_hy3_lagrange2_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, Half, unsigned int> grid_transfer_truncate_test_sx2_lagrange2_half_uint(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, Half, unsigned int> grid_transfer_truncate_test_sx3_lagrange2_half_uint(PreferredBackend::generic);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, Half, unsigned long> grid_transfer_truncate_test_hy1_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, Half, unsigned long> grid_transfer_truncate_test_hy2_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, Half, unsigned long> grid_transfer_truncate_test_hy3_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, Half, unsigned long> grid_transfer_truncate_test_sx2_lagrange2_half_ulong(PreferredBackend::generic);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, Half, unsigned long> grid_transfer_truncate_test_sx3_lagrange2_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
// Lagrange-1 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned int> cuda_grid_transfer_truncate_test_hy1_lagrange1_double_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned int> cuda_grid_transfer_truncate_test_hy2_lagrange1_float_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned int> cuda_grid_transfer_truncate_test_hy3_lagrange1_double_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned int> cuda_grid_transfer_truncate_test_sx2_lagrange1_double_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned int> cuda_grid_transfer_truncate_test_sx3_lagrange1_float_uint(PreferredBackend::cuda);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange1::Element, 4, double, unsigned long> cuda_grid_transfer_truncate_test_hy1_lagrange1_double_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange1::Element, 2, float, unsigned long> cuda_grid_transfer_truncate_test_hy2_lagrange1_float_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange1::Element, 1, double, unsigned long> cuda_grid_transfer_truncate_test_hy3_lagrange1_double_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange1::Element, 2, double, unsigned long> cuda_grid_transfer_truncate_test_sx2_lagrange1_double_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange1::Element, 1, float, unsigned long> cuda_grid_transfer_truncate_test_sx3_lagrange1_float_ulong(PreferredBackend::cuda);
// Lagrange-2 element
GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned int> cuda_grid_transfer_truncate_test_hy1_lagrange2_double_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned int> cuda_grid_transfer_truncate_test_hy2_lagrange2_double_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned int> cuda_grid_transfer_truncate_test_hy3_lagrange2_float_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned int> cuda_grid_transfer_truncate_test_sx2_lagrange2_float_uint(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned int> cuda_grid_transfer_truncate_test_sx3_lagrange2_double_uint(PreferredBackend::cuda);

GridTransferTipTest<Shape::Hypercube<1>, Space::Lagrange2::Element, 4, double, unsigned long> cuda_grid_transfer_truncate_test_hy1_lagrange2_double_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<2>, Space::Lagrange2::Element, 2, double, unsigned long> cuda_grid_transfer_truncate_test_hy2_lagrange2_double_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Hypercube<3>, Space::Lagrange2::Element, 1, float, unsigned long> cuda_grid_transfer_truncate_test_hy3_lagrange2_float_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<2>, Space::Lagrange2::Element, 2, float, unsigned long> cuda_grid_transfer_truncate_test_sx2_lagrange2_float_ulong(PreferredBackend::cuda);
GridTransferTipTest<Shape::Simplex<3>, Space::Lagrange2::Element, 1, double, unsigned long> cuda_grid_transfer_truncate_test_sx3_lagrange2_double_ulong(PreferredBackend::cuda);
#endif

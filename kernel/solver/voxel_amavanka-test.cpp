// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/unit_cube_patch_generator.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/solver/amavanka_base.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/voxel_amavanka.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/gpdv_assembler.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/util/time_stamp.hpp>


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

template<typename T_>
struct VeloFuncX
{
  static T_ eval (T_, T_ y) {return y * (T_(1) - y);}
};

template<typename Trafo_>
struct RTQ0
{
  typedef Space::CroRavRanTur::Element<Trafo_> V;
  typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<0>> P;
  static const char* name() {return "RT/Q0";}
};

template<typename Trafo_>
struct Q2P1
{
  typedef Space::Lagrange2::Element<Trafo_> V;
  typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<1>> P;
  static const char* name() {return "Q2/P1";}
};

template<typename DT_, typename IT_>
class MatrixWrapperTest :
  public UnitTest
{
public:
  explicit MatrixWrapperTest(PreferredBackend backend) :
    UnitTest("BaseAmaVankaMatrixWrapperTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~MatrixWrapperTest()
  {
  }

  template<int dim_>
  void test() const
  {
    //init big meta matrix
    typedef LAFEM::NullMatrix<DT_, IT_> BaseZeroScalar;
    // typedef LAFEM::NullMatrix<DT_, IT_, dim_, 1> BaseZeroRect;
    // typedef LAFEM::NullMatrix<DT_, IT_,  1, dim_> BaseZeroRectT;
    typedef LAFEM::SparseMatrixCSR<DT_, IT_> BaseCSR;
    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, dim_> BaseBCSRSquare;
    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, 1> BaseBCSRRect;
    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, 1, dim_> BaseBCSRRectT;
    typedef SaddlePointMatrix<BaseBCSRSquare, BaseBCSRRect, BaseBCSRRectT> SaddlePointType;
    //base test
    {
      BaseZeroScalar test_mat(30, 59);
      auto mat_wrapper = Solver::Intern::get_meta_matrix_wrapper(test_mat, Memory::Location::main);
      TEST_CHECK_EQUAL(mat_wrapper.n, 1);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[0], nullptr);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[0], nullptr);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[0], nullptr);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[0], IT_(0));
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[0], test_mat.num_rows());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[1], test_mat.num_cols());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[0], 1);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[1], 1);
    }
    {
      DenseVector<DT_, IT_> dv1(2*dim_*dim_);
      DenseVector<IT_, IT_> dv2(2);
      DenseVector<IT_, IT_> dv3(3);
      {
        Memory::TypedView<DT_> v1(dv1.elements_view_w());
        Memory::TypedView<IT_> v2(dv2.elements_view_w());
        Memory::TypedView<IT_> v3(dv3.elements_view_w());
        for (Index i(0) ; i < dv1.size() ; ++i)
          v1[i] = DT_(i+1);
        v2[0] = IT_(0);
        v2[1] = IT_(1);
        v3[0] = IT_(0);
        v3[1] = IT_(1);
        v3[2] = IT_(2);
      }
      SparseMatrixBCSR<DT_, IT_, dim_, dim_> test_mat(2, 2, dv3, dv2, dv1);
      auto mat_wrapper = Solver::Intern::get_meta_matrix_wrapper(test_mat, Memory::Location::main);
      TEST_CHECK_EQUAL(mat_wrapper.n, 1);
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[0], test_mat.val_view_raw_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[0], test_mat.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[0], test_mat.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[0], test_mat.num_nzes());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[0], test_mat.num_rows());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[1], test_mat.num_cols());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[0], dim_);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[1], dim_);
    }
    {
      DenseVector<DT_, IT_> dv1(2*dim_*dim_);
      DenseVector<DT_, IT_> dv4(2*dim_);
      DenseVector<DT_, IT_> dv5(2);
      DenseVector<IT_, IT_> dv2(2);
      DenseVector<IT_, IT_> dv3(3);
      {
        Memory::TypedView<DT_> v1(dv1.elements_view_w());
        Memory::TypedView<IT_> v2(dv2.elements_view_w());
        Memory::TypedView<IT_> v3(dv3.elements_view_w());
        Memory::TypedView<DT_> v4(dv4.elements_view_w());
        Memory::TypedView<DT_> v5(dv5.elements_view_w());
        for (Index i(0) ; i < dv1.size() ; ++i)
          v1[i] = DT_(i+1);
        for (Index i(0) ; i < dv4.size() ; ++i)
          v4[i] = DT_(i+1);
        for (Index i(0) ; i < dv5.size() ; ++i)
          v5[i] = DT_(i+1);
        v2[0] = IT_(0);
        v2[1] = IT_(1);
        v3[0] = IT_(0);
        v3[1] = IT_(1);
        v3[2] = IT_(2);
      }

      BaseBCSRSquare test_mat(2, 2, dv3, dv2, dv1);
      BaseBCSRRect test_mat2(2, 2, dv3, dv2, dv4);
      BaseBCSRRectT test_mat3(2, 2, dv3, dv2, dv4);
      BaseCSR test_mat4(2, 2, dv3, dv2, dv5);
      LAFEM::TupleMatrix<LAFEM::TupleMatrixRow<BaseBCSRSquare, BaseBCSRRect>, LAFEM::TupleMatrixRow<BaseBCSRRectT, BaseCSR>> meta_test_mat;
      meta_test_mat.template at<0,0>().clone(test_mat, LAFEM::CloneMode::Weak);
      meta_test_mat.template at<0,1>().clone(test_mat2, LAFEM::CloneMode::Weak);
      meta_test_mat.template at<1,0>().clone(test_mat3, LAFEM::CloneMode::Weak);
      meta_test_mat.template at<1,1>().clone(test_mat4, LAFEM::CloneMode::Weak);
      auto mat_wrapper = Solver::Intern::get_meta_matrix_wrapper(meta_test_mat, Memory::Location::main);
      TEST_CHECK_EQUAL(mat_wrapper.n, 2);
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[0], test_mat.val_view_raw_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[0], test_mat.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[0], test_mat.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[0], test_mat.num_nzes());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[0], test_mat.num_rows());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[1], test_mat.num_cols());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[0], dim_);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[1], dim_);
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[1], test_mat2.val_view_raw_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[1], test_mat2.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[1], test_mat2.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[1], test_mat2.num_nzes());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[2], test_mat2.num_cols());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[2], 1);
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[2], test_mat3.val_view_raw_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[2], test_mat3.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[2], test_mat3.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[2], test_mat3.num_nzes());
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[3], test_mat4.val_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[3], test_mat4.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[3], test_mat4.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[3], test_mat4.num_nzes());
    }
    {
      DenseVector<DT_, IT_> dv1(2*dim_*dim_);
      DenseVector<DT_, IT_> dv4(2*dim_);
      DenseVector<IT_, IT_> dv2(2);
      DenseVector<IT_, IT_> dv3(3);
      {
        Memory::TypedView<DT_> v1(dv1.elements_view_w());
        Memory::TypedView<IT_> v2(dv2.elements_view_w());
        Memory::TypedView<IT_> v3(dv3.elements_view_w());
        Memory::TypedView<DT_> v4(dv4.elements_view_w());
        for (Index i(0) ; i < dv1.size() ; ++i)
          v1[i] = DT_(i+1);
        for (Index i(0) ; i < dv4.size() ; ++i)
          v4[i] = DT_(i+1);
        v2[0] = IT_(0);
        v2[1] = IT_(1);
        v3[0] = IT_(0);
        v3[1] = IT_(1);
        v3[2] = IT_(2);
      }

      BaseBCSRSquare test_mat(2, 2, dv3, dv2, dv1);
      BaseBCSRRect test_mat2(2, 2, dv3, dv2, dv4);
      BaseBCSRRectT test_mat3(2, 2, dv3, dv2, dv4);
      SaddlePointType meta_test_mat;
      meta_test_mat.template at<0,0>().clone(test_mat, LAFEM::CloneMode::Weak);
      meta_test_mat.template at<0,1>().clone(test_mat2, LAFEM::CloneMode::Weak);
      meta_test_mat.template at<1,0>().clone(test_mat3, LAFEM::CloneMode::Weak);
      auto mat_wrapper = Solver::Intern::get_meta_matrix_wrapper(meta_test_mat, Memory::Location::main);
      TEST_CHECK_EQUAL(mat_wrapper.n, 2);
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[0], test_mat.val_view_raw_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[0], test_mat.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[0], test_mat.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[0], test_mat.num_nzes());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[0], test_mat.num_rows());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[1], test_mat.num_cols());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[0], dim_);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[1], dim_);
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[1], test_mat2.val_view_raw_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[1], test_mat2.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[1], test_mat2.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[1], test_mat2.num_nzes());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.tensor_counts[2], test_mat2.num_cols());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.blocksizes[2], 1);
      //TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[2], test_mat3.val_view_raw_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[2], test_mat3.col_idx_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[2], test_mat3.row_ptr_view_r().get_r());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[2], test_mat3.num_nzes());
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.data_arrays[3], nullptr);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.col_arrays[3], nullptr);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.row_arrays[3], nullptr);
      TEST_CHECK_EQUAL(mat_wrapper.raw_data.used_elements[3], IT_(0));
    }
  }

  virtual void run() const override
  {
    test<2>();
  }
}; // MatrixWrapperTest

SPAWN_UNIT_TEST_2T_P(MatrixWrapperTest, double, std::uint64_t, PreferredBackend::generic);

template<typename DT_, typename IT_, int dim_, FEAT::Intern::VankaAssemblyPolicy asm_pol_>
class VoxelAmaVankaTest :
  public UnitTest
{
public:
  explicit VoxelAmaVankaTest(PreferredBackend backend) :
    UnitTest("VoxelAmaVankaTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~VoxelAmaVankaTest()
  {
  }

  static constexpr int dim = dim_;
  typedef Shape::Hypercube<dim> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNode;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  /// tests Vanka for Deformation tensor formulation
  template<template<typename> class Space_, bool uniform_>
  void test_defo(MeshNode& mesh_node, const DT_ omega) const
  {
    // get our mesh
    MeshType& mesh = *mesh_node.get_mesh();

    // create our trafo
    TrafoType trafo(mesh);

    typename Space_<TrafoType>::V space_v(trafo);
    typename Space_<TrafoType>::P space_p(trafo);
    const String name = Space_<TrafoType>::name();

    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, dim, dim> MatrixTypeA;
    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, dim, 1> MatrixTypeB;
    typedef LAFEM::SparseMatrixBCSR<DT_, IT_, 1, dim> MatrixTypeD;
    typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;

    typedef typename MatrixType::VectorTypeL VectorType;

    typedef LAFEM::UnitFilterBlocked<DT_, IT_, dim> VeloFilterType;
    typedef LAFEM::NoneFilter<DT_, IT_> PresFilterType;
    typedef LAFEM::TupleFilter<VeloFilterType, PresFilterType> FilterType;

    MatrixType matrix;
    MatrixTypeA& mat_a = matrix.block_a();
    MatrixTypeB& mat_b = matrix.block_b();
    MatrixTypeD& mat_d = matrix.block_d();
    FilterType filter;

    Cubature::DynamicFactory cubature("auto-degree:5");

    // assemble A, B and D
    Assembly::SymbolicAssembler::assemble_matrix_std1(mat_a, space_v);
    mat_a.format();

    Assembly::Common::DuDvOperatorBlocked<dim> dudv_op;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a, dudv_op, space_v, cubature);
    Assembly::GradPresDivVeloAssembler::assemble(mat_b, mat_d, space_v, space_p, cubature);

    // assemble filter
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:0"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:1"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:2"));

    // finally, assemble the filters
    if constexpr(dim ==2)
    {
      Analytic::Common::ParProfileVector<DT_> inflow_func(0.0, 0.0, 0.0, 1.0, 1.0);
      unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);
    }
    else
    {
      Analytic::Common::PoiseuillePipeFlow<DT_, dim> inflow_func(Tiny::Vector<DT_, dim>{DT_(0.0), DT_(0.5), DT_(0.5)}, Tiny::Vector<DT_, dim>{DT_(1.0), DT_(0.), DT_(0.)}, DT_(0.5), DT_(1.3));
      unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);
    }

    // filter matrices
    filter.template at<0>().filter_mat(mat_a);
    filter.template at<0>().filter_offdiag_row_mat(mat_b);

    // create RHS and SOL
    VectorType vec_rhs = matrix.create_vector_l();
    VectorType vec_sol = matrix.create_vector_l();
    VectorType vec_def = matrix.create_vector_l();

    vec_rhs.format();
    vec_sol.format();

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    //create coloring from our mesh
    const auto& verts_at_elem = mesh.template get_index_set<dim, 0>();
    Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);
    Adjacency::Graph elems_at_elem(Adjacency::RenderType::injectify, verts_at_elem, elems_at_vert);

    // create coloring
    Adjacency::Coloring col(elems_at_elem);

    // create vanka
    auto vanka = Solver::new_amavanka(matrix, filter, omega, 10);
    constexpr FEAT::Intern::VankaMacroPolicy macro_pol = uniform_ ? FEAT::Intern::VankaMacroPolicy::uniformMacros : FEAT::Intern::VankaMacroPolicy::anisotropicMacros;
    auto voxel_vanka = Solver::new_voxel_amavanka<MatrixType, FilterType, decltype(col), asm_pol_, macro_pol>(matrix, filter, col, omega, 10);

    vanka->init();
    voxel_vanka->init();
    // // voxel_vanka->init_symbolic();

    TEST_CHECK(voxel_vanka->get_vanka_matrix().same_layout(vanka->get_vanka_matrix()));

    // // voxel_vanka->done_symbolic();
    voxel_vanka->done();
    vanka->done();


    vanka->set_skip_singular(true);
    voxel_vanka->set_skip_singular(true);

    vanka->init();
    voxel_vanka->init();

    TEST_CHECK(voxel_vanka->get_vanka_matrix().same_layout(vanka->get_vanka_matrix()));

    voxel_vanka->done();
    vanka->done();



  }

  virtual void run() const override
  {
    const int level = 4-dim;

    // create mesh node
    std::vector<int> ranks;
    std::unique_ptr<MeshNode> mesh_node;
    Geometry::UnitCubePatchGenerator<MeshType>::create_unique(0, 1, mesh_node, ranks);

    // refine a few times
    for(int i = 0; i < level; ++i)
    {
      mesh_node = mesh_node->refine_unique();
    }

    // test RT/Q0
    test_defo<RTQ0, true>(*mesh_node, DT_(1.0));
    if constexpr(FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock == asm_pol_)
      test_defo<RTQ0, false>(*mesh_node, DT_(1.0));

    // test Q2/P1dc
    test_defo<Q2P1, true>(*mesh_node, DT_(1.0));
    if constexpr(FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock == asm_pol_)
      test_defo<Q2P1, false>(*mesh_node, DT_(1.0));
  }
};

SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint32_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint64_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint32_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint64_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint32_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint64_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint32_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint64_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint32_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint64_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint32_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint64_t, 2, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint32_t, 2, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint64_t, 2, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint32_t, 2, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint64_t, 2, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint32_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint64_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint32_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint64_t, 3, FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint32_t, 3, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, double, std::uint64_t, 3, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint32_t, 3, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
SPAWN_UNIT_TEST_4T_P(VoxelAmaVankaTest, float, std::uint64_t, 3, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, PreferredBackend::cuda);
#endif

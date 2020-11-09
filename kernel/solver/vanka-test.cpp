// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
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
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/power_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/solver/vanka.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/util/time_stamp.hpp>

#include <iomanip>

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

template<typename Trafo_>
struct Q2Q1
{
  typedef Space::Lagrange2::Element<Trafo_> V;
  typedef Space::Lagrange1::Element<Trafo_> P;
  static const char* name() {return "Q2/Q1";}
};

template<typename DT_, typename IT_>
class VankaTest :
  public FullTaggedTest<Mem::Main, DT_, IT_>
{
public:
  VankaTest() :
    FullTaggedTest<Mem::Main, DT_, IT_>("VankaTest")
  {
  }

  virtual ~VankaTest()
  {
  }

  static constexpr int dim = 2;
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNode;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  /// tests Vanka for Gradient tensor formulation
  template<template<typename> class Space_>
  void test_grad(MeshNode& mesh_node, Solver::VankaType vtype, const DT_ omega) const
  {
    // get our mesh
    MeshType& mesh = *mesh_node.get_mesh();

    // create our trafo
    TrafoType trafo(mesh);

    typename Space_<TrafoType>::V space_v(trafo);
    typename Space_<TrafoType>::P space_p(trafo);
    const String name = Space_<TrafoType>::name();

    typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> ScalarMatrixType;
    typedef LAFEM::PowerDiagMatrix<ScalarMatrixType, dim> MatrixTypeA;
    typedef LAFEM::PowerColMatrix<ScalarMatrixType, dim> MatrixTypeB;
    typedef LAFEM::PowerRowMatrix<ScalarMatrixType, dim> MatrixTypeD;
    typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;
    typedef typename MatrixType::VectorTypeL VectorType;
    typedef LAFEM::NoneFilter<Mem::Main, double, Index> NoneFilterType;
    typedef LAFEM::UnitFilter<Mem::Main, double, Index> UnitFilterType;
    typedef LAFEM::PowerFilter<UnitFilterType, dim> VeloFilterType;
    typedef LAFEM::TupleFilter<VeloFilterType, NoneFilterType> FilterType;

    MatrixType matrix;
    MatrixTypeA& mat_a = matrix.block_a();
    MatrixTypeB& mat_b = matrix.block_b();
    MatrixTypeD& mat_d = matrix.block_d();
    FilterType filter;

    Cubature::DynamicFactory cubature("auto-degree:5");

    // assemble A, B and D
    {
      ScalarMatrixType mat;
      Assembly::SymbolicAssembler::assemble_matrix_std1(mat, space_v);
      mat_a = MatrixTypeA(mat.layout());
      mat_a.format();
      Assembly::Common::LaplaceOperator oper;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.template at<0,0>(), oper, space_v, cubature);
      mat_a.template at<1,1>().copy(mat_a.template at<0,0>());
    }
    {
      ScalarMatrixType mat;
      Assembly::SymbolicAssembler::assemble_matrix_std2(mat, space_v, space_p);
      mat_b = MatrixTypeB(mat.layout());
      mat_b.format();
      Assembly::Common::TestDerivativeOperator oper_x(0);
      Assembly::Common::TestDerivativeOperator oper_y(1);
      Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_b.template at<0,0>(), oper_x, space_v, space_p, cubature, -1.0);
      Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_b.template at<1,0>(), oper_y, space_v, space_p, cubature, -1.0);
    }
    {
      mat_d.template at<0,0>() = mat_b.template at<0,0>().transpose();
      mat_d.template at<0,1>() = mat_b.template at<1,0>().transpose();
    }

    // assemble filter
    {
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      // add mesh parts 0,1,2
      unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:0"));
      unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:1"));
      unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:2"));

      // finally, assemble the filters
      Analytic::StaticWrapperFunction<2, VeloFuncX> inflow_func;
      unit_asm.assemble(filter.template at<0>().template at<0>(), space_v, inflow_func);
      unit_asm.assemble(filter.template at<0>().template at<1>(), space_v);
    }

    // filter matrices
    filter.template at<0>().template at<0>().filter_mat(mat_a.template at<0,0>());
    filter.template at<0>().template at<1>().filter_mat(mat_a.template at<1,1>());

    // create RHS and SOL
    VectorType vec_rhs = matrix.create_vector_l();
    VectorType vec_sol = matrix.create_vector_l();
    VectorType vec_def = matrix.create_vector_l();

    vec_rhs.format();
    vec_sol.format();

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    // compute initial defect
    matrix.apply(vec_def, vec_sol, vec_rhs, -DT_(1));
    const DT_ def0 = vec_def.norm2();

    TimeStamp stamp1, stamp2;
    {
      // create vanka
      auto vanka = Solver::new_vanka(matrix, filter, vtype, omega, 10);

      // initialize and solve
      vanka->init();
      stamp1.stamp();
      Solver::solve(*vanka, vec_sol, vec_rhs, matrix, filter);
      stamp2.stamp();
    }

    // compute final defect
    matrix.apply(vec_def, vec_sol, vec_rhs, -DT_(1));
    const DT_ def1 = vec_def.norm2();

    // print some information
    std::cout << "Grad: " << name << ": ";
    switch(vtype)
    {
    case Solver::VankaType::nodal_diag_mult: std::cout << "NDM: "; break;
    case Solver::VankaType::nodal_full_mult: std::cout << "NFM: "; break;
    case Solver::VankaType::block_diag_mult: std::cout << "BDM: "; break;
    case Solver::VankaType::block_full_mult: std::cout << "BFM: "; break;
    case Solver::VankaType::nodal_diag_add:  std::cout << "NDA: "; break;
    case Solver::VankaType::nodal_full_add:  std::cout << "NFA: "; break;
    case Solver::VankaType::block_diag_add:  std::cout << "BDA: "; break;
    case Solver::VankaType::block_full_add:  std::cout << "BFA: "; break;
    }
    std::cout << stringify_fp_sci(def0) << " > " << stringify_fp_sci(def1);
    std::cout << " : " << std::fixed << std::setprecision(5) << (def1/def0);
    std::cout << " | " << stamp2.elapsed_string(stamp1, TimeFormat::s_m);
    std::cout << std::endl;

    // ensure that the defect decreased
    TEST_CHECK_IN_RANGE(def1/def0, DT_(0), DT_(0.7));
  }

  /// tests Vanka for Deformation tensor formulation
  template<template<typename> class Space_>
  void test_defo(MeshNode& mesh_node, Solver::VankaType vtype, const DT_ omega) const
  {
    // get our mesh
    MeshType& mesh = *mesh_node.get_mesh();

    // create our trafo
    TrafoType trafo(mesh);

    typename Space_<TrafoType>::V space_v(trafo);
    typename Space_<TrafoType>::P space_p(trafo);
    const String name = Space_<TrafoType>::name();

    typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> ScalarMatrixType;
    typedef LAFEM::PowerFullMatrix<ScalarMatrixType, dim, dim> MatrixTypeA;
    typedef LAFEM::PowerColMatrix<ScalarMatrixType, dim> MatrixTypeB;
    typedef LAFEM::PowerRowMatrix<ScalarMatrixType, dim> MatrixTypeD;
    typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;
    typedef typename MatrixType::VectorTypeL VectorType;
    typedef LAFEM::NoneFilter<Mem::Main, double, Index> NoneFilterType;
    typedef LAFEM::UnitFilter<Mem::Main, double, Index> UnitFilterType;
    typedef LAFEM::PowerFilter<UnitFilterType, dim> VeloFilterType;
    typedef LAFEM::TupleFilter<VeloFilterType, NoneFilterType> FilterType;

    MatrixType matrix;
    MatrixTypeA& mat_a = matrix.block_a();
    MatrixTypeB& mat_b = matrix.block_b();
    MatrixTypeD& mat_d = matrix.block_d();
    FilterType filter;

    Cubature::DynamicFactory cubature("auto-degree:5");

    // assemble A, B and D
    {
      ScalarMatrixType mat;
      Assembly::SymbolicAssembler::assemble_matrix_std1(mat, space_v);
      mat_a = MatrixTypeA(mat.layout());
      mat_a.format();
      Assembly::Common::DuDvOperator dudv00(0,0);
      Assembly::Common::DuDvOperator dudv01(0,1);
      Assembly::Common::DuDvOperator dudv10(1,0);
      Assembly::Common::DuDvOperator dudv11(1,1);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.template at<0,0>(), dudv00, space_v, cubature, 0.5);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.template at<0,1>(), dudv01, space_v, cubature, 0.5);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.template at<1,0>(), dudv10, space_v, cubature, 0.5);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.template at<1,1>(), dudv11, space_v, cubature, 0.5);
    }
    {
      ScalarMatrixType mat;
      Assembly::SymbolicAssembler::assemble_matrix_std2(mat, space_v, space_p);
      mat_b = MatrixTypeB(mat.layout());
      mat_b.format();
      Assembly::Common::TestDerivativeOperator oper_x(0);
      Assembly::Common::TestDerivativeOperator oper_y(1);
      Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_b.template at<0,0>(), oper_x, space_v, space_p, cubature, -1.0);
      Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_b.template at<1,0>(), oper_y, space_v, space_p, cubature, -1.0);
    }
    {
      mat_d.template at<0,0>() = mat_b.template at<0,0>().transpose();
      mat_d.template at<0,1>() = mat_b.template at<1,0>().transpose();
    }

    // assemble filter
    {
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      // add mesh parts 0,1,2
      unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:0"));
      unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:1"));
      unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:2"));

      // finally, assemble the filters
      Analytic::StaticWrapperFunction<2, VeloFuncX> inflow_func;
      unit_asm.assemble(filter.template at<0>().template at<0>(), space_v, inflow_func);
      unit_asm.assemble(filter.template at<0>().template at<1>(), space_v);
    }

    // filter matrices
    filter.template at<0>().template at<0>().filter_mat(mat_a.template at<0,0>());
    filter.template at<0>().template at<1>().filter_mat(mat_a.template at<1,1>());

    filter.template at<0>().template at<0>().filter_offdiag_row_mat(mat_a.template at<0,1>());
    filter.template at<0>().template at<1>().filter_offdiag_row_mat(mat_a.template at<1,0>());

    // create RHS and SOL
    VectorType vec_rhs = matrix.create_vector_l();
    VectorType vec_sol = matrix.create_vector_l();
    VectorType vec_def = matrix.create_vector_l();

    vec_rhs.format();
    vec_sol.format();

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    // compute initial defect
    matrix.apply(vec_def, vec_sol, vec_rhs, -DT_(1));
    const DT_ def0 = vec_def.norm2();

    TimeStamp stamp1, stamp2;
    {
      // create vanka
      auto vanka = Solver::new_vanka(matrix, filter, vtype, omega, 10);

      // initialize and solve
      vanka->init();
      stamp1.stamp();
      Solver::solve(*vanka, vec_sol, vec_rhs, matrix, filter);
      stamp2.stamp();
    }

    // compute final defect
    matrix.apply(vec_def, vec_sol, vec_rhs, -DT_(1));
    const DT_ def1 = vec_def.norm2();

    // print some information
    std::cout << "Defo: " << name << ": ";
    switch(vtype)
    {
    case Solver::VankaType::nodal_diag_mult: std::cout << "NDM: "; break;
    case Solver::VankaType::nodal_full_mult: std::cout << "NFM: "; break;
    case Solver::VankaType::block_diag_mult: std::cout << "BDM: "; break;
    case Solver::VankaType::block_full_mult: std::cout << "BFM: "; break;
    case Solver::VankaType::nodal_diag_add:  std::cout << "NDA: "; break;
    case Solver::VankaType::nodal_full_add:  std::cout << "NFA: "; break;
    case Solver::VankaType::block_diag_add:  std::cout << "BDA: "; break;
    case Solver::VankaType::block_full_add:  std::cout << "BFA: "; break;
    }
    std::cout << stringify_fp_sci(def0) << " > " << stringify_fp_sci(def1);
    std::cout << " : " << std::fixed << std::setprecision(5) << (def1/def0);
    std::cout << " | " << stamp2.elapsed_string(stamp1, TimeFormat::s_m);
    std::cout << std::endl;

    // ensure that the defect decreased
    TEST_CHECK_IN_RANGE(def1/def0, DT_(0), DT_(0.7));
  }

  virtual void run() const override
  {
    const int level = 2;

    // create mesh node
    std::shared_ptr<MeshNode> mesh_node;
    std::vector<int> ranks;
    Geometry::UnitCubePatchGenerator<MeshType>::create(0, 1, mesh_node, ranks);

    // refine a few times
    for(int i = 0; i < level; ++i)
    {
      mesh_node = std::shared_ptr<MeshNode>(mesh_node->refine());
    }

    // test RT/Q0
    test_grad<RTQ0>(*mesh_node, VankaType::nodal_diag_mult, DT_(1.0));
    test_grad<RTQ0>(*mesh_node, VankaType::block_full_mult, DT_(1.0));
    test_grad<RTQ0>(*mesh_node, VankaType::nodal_diag_add,  DT_(0.9)); // diverges for 1.0
    test_grad<RTQ0>(*mesh_node, VankaType::block_full_add,  DT_(1.0));

    // test Q2/P1dc
    test_grad<Q2P1>(*mesh_node, VankaType::nodal_diag_mult, DT_(1.0));
    test_grad<Q2P1>(*mesh_node, VankaType::nodal_full_mult, DT_(1.0));
    test_grad<Q2P1>(*mesh_node, VankaType::block_diag_mult, DT_(1.0));
    test_grad<Q2P1>(*mesh_node, VankaType::block_full_mult, DT_(1.0));
    test_grad<Q2P1>(*mesh_node, VankaType::nodal_diag_add,  DT_(1.0));
    test_grad<Q2P1>(*mesh_node, VankaType::nodal_full_add,  DT_(1.0));
    test_grad<Q2P1>(*mesh_node, VankaType::block_diag_add,  DT_(1.0));
    test_grad<Q2P1>(*mesh_node, VankaType::block_full_add,  DT_(1.0));

    // test Q2/Q1
    test_grad<Q2Q1>(*mesh_node, VankaType::nodal_diag_mult, DT_(0.9)); // diverges for 1.0
    test_grad<Q2Q1>(*mesh_node, VankaType::nodal_full_mult, DT_(1.0));
    test_grad<Q2Q1>(*mesh_node, VankaType::block_diag_mult, DT_(1.0));
    test_grad<Q2Q1>(*mesh_node, VankaType::block_full_mult, DT_(1.0));
    test_grad<Q2Q1>(*mesh_node, VankaType::nodal_diag_add,  DT_(1.0));
    test_grad<Q2Q1>(*mesh_node, VankaType::nodal_full_add,  DT_(1.0));
    test_grad<Q2Q1>(*mesh_node, VankaType::block_diag_add,  DT_(1.0));
    test_grad<Q2Q1>(*mesh_node, VankaType::block_full_add,  DT_(1.0));

    // test RT/Q0
    test_defo<RTQ0>(*mesh_node, VankaType::nodal_diag_mult, DT_(1.0));
    test_defo<RTQ0>(*mesh_node, VankaType::block_full_mult, DT_(1.0));
    test_defo<RTQ0>(*mesh_node, VankaType::nodal_diag_add,  DT_(0.9)); // diverges for 1.0
    test_defo<RTQ0>(*mesh_node, VankaType::block_full_add,  DT_(1.0));

    // test Q2/P1dc
    test_defo<Q2P1>(*mesh_node, VankaType::nodal_diag_mult, DT_(1.0));
    test_defo<Q2P1>(*mesh_node, VankaType::nodal_full_mult, DT_(1.0));
    test_defo<Q2P1>(*mesh_node, VankaType::block_diag_mult, DT_(1.0));
    test_defo<Q2P1>(*mesh_node, VankaType::block_full_mult, DT_(1.0));
    test_defo<Q2P1>(*mesh_node, VankaType::nodal_diag_add,  DT_(0.9)); // diverges for 1.0
    test_defo<Q2P1>(*mesh_node, VankaType::nodal_full_add,  DT_(1.0));
    test_defo<Q2P1>(*mesh_node, VankaType::block_diag_add,  DT_(1.0));
    test_defo<Q2P1>(*mesh_node, VankaType::block_full_add,  DT_(1.0));

    // test Q2/Q1
    test_defo<Q2Q1>(*mesh_node, VankaType::nodal_diag_mult, DT_(1.0));
    test_defo<Q2Q1>(*mesh_node, VankaType::nodal_full_mult, DT_(1.0));
    test_defo<Q2Q1>(*mesh_node, VankaType::block_diag_mult, DT_(1.0));
    test_defo<Q2Q1>(*mesh_node, VankaType::block_full_mult, DT_(1.0));
    test_defo<Q2Q1>(*mesh_node, VankaType::nodal_diag_add,  DT_(1.0));
    test_defo<Q2Q1>(*mesh_node, VankaType::nodal_full_add,  DT_(1.0));
    test_defo<Q2Q1>(*mesh_node, VankaType::block_diag_add,  DT_(1.0));
    test_defo<Q2Q1>(*mesh_node, VankaType::block_full_add,  DT_(1.0));
  }
};

VankaTest<double, Index> vanka_test_double_index;

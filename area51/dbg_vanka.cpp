// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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

using namespace FEAT;

template<typename T_>
struct VeloFuncX
{
  static T_ eval (T_, T_ y) {return y * (T_(1) - y);}
  static T_ der_x(T_, T_  ) {return T_(0);}
  static T_ der_y(T_, T_ y) {return T_(1) - T_(2)*y;}
};

template<typename T_>
struct VeloFuncY
{
  static T_ eval (T_, T_) {return T_(0);}
  static T_ der_x(T_, T_) {return T_(0);}
  static T_ der_y(T_, T_) {return T_(0);}
};

template<typename T_>
struct PresFunc
{
  static T_ eval (T_ x, T_) {return T_(2)*(T_(1) - x);}
};

template<typename Trafo_>
struct SpaceRTQ0
{
  typedef Space::CroRavRanTur::Element<Trafo_> V;
  typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<0>> P;
  static const char* name() {return "RT/Q0";}
};

template<typename Trafo_>
struct SpaceQ2P1
{
  typedef Space::Lagrange2::Element<Trafo_> V;
  typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<1>> P;
  static const char* name() {return "Q2/P1";}
};

template<typename Trafo_>
struct SpaceQ2Q1
{
  typedef Space::Lagrange2::Element<Trafo_> V;
  typedef Space::Lagrange1::Element<Trafo_> P;
  static const char* name() {return "Q2/Q1";}
};
//template<typename SpaceV_, typename SpaceP_>
template<template<typename> class Space_>
void test_poiseuille(int level, bool defo, Solver::VankaType vtype)
{
  static constexpr int dim = 2;

  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNode;

  // create mesh node
  MeshNode* mesh_node = nullptr;
  std::vector<Index> ranks, ctags;
  Geometry::UnitCubePatchGenerator<MeshType>::create(0, 1, mesh_node, ranks, ctags);

  // refine a few times
  for(int i = 0; i < level; ++i)
  {
    MeshNode* old_node = mesh_node;
    mesh_node = old_node->refine();
    delete old_node;
  }

  // get our mesh
  MeshType& mesh = *mesh_node->get_mesh();

  // create our trafo
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  TrafoType trafo(mesh);

  // create our spaces
  //typedef Space::RannacherTurek::Element<TrafoType> SpaceV;
  //typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<0>> SpaceP;

  //typedef Space::Lagrange2::Element<TrafoType> SpaceV;
  //typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpaceP;


  //SpaceV space_v(trafo);
  //SpaceP space_p(trafo);
  typename Space_<TrafoType>::V space_v(trafo);
  typename Space_<TrafoType>::P space_p(trafo);
  const String name = Space_<TrafoType>::name();

  typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> ScalarMatrixType;
  //typedef LAFEM::PowerDiagMatrix<ScalarMatrixType, dim> MatrixTypeA;
  typedef LAFEM::PowerFullMatrix<ScalarMatrixType, dim, dim> MatrixTypeA;
  typedef LAFEM::PowerColMatrix<ScalarMatrixType, dim> MatrixTypeB;
  typedef LAFEM::PowerRowMatrix<ScalarMatrixType, dim> MatrixTypeD;
  typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;
  typedef MatrixType::VectorTypeL VectorType;
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
    if(defo)
    {
      Assembly::Common::DuDvOperator dudv00(0,0);
      Assembly::Common::DuDvOperator dudv01(0,1);
      Assembly::Common::DuDvOperator dudv10(1,0);
      Assembly::Common::DuDvOperator dudv11(1,1);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.at<0,0>(), dudv00, space_v, cubature, 0.5);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.at<0,1>(), dudv01, space_v, cubature, 0.5);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.at<1,0>(), dudv10, space_v, cubature, 0.5);
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.at<1,1>(), dudv11, space_v, cubature, 0.5);
    }
    else
    {
      Assembly::Common::LaplaceOperator oper;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a.at<0,0>(), oper, space_v, cubature);
      mat_a.at<1,1>().copy(mat_a.at<0,0>());
      //mat_a.at<0,1>().clear();
      //mat_a.at<1,0>().clear();
    }
  }
  {
    ScalarMatrixType mat;
    Assembly::SymbolicAssembler::assemble_matrix_std2(mat, space_v, space_p);
    mat_b = MatrixTypeB(mat.layout());
    mat_b.format();
    Assembly::Common::TestDerivativeOperator oper_x(0);
    Assembly::Common::TestDerivativeOperator oper_y(1);
    Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_b.at<0,0>(), oper_x, space_v, space_p, cubature, -1.0);
    Assembly::BilinearOperatorAssembler::assemble_matrix2(mat_b.at<1,0>(), oper_y, space_v, space_p, cubature, -1.0);
  }
  {
    mat_d.at<0,0>() = mat_b.at<0,0>().transpose();
    mat_d.at<0,1>() = mat_b.at<1,0>().transpose();
  }

  {
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    // add mesh parts 0,1,2
    unit_asm.add_mesh_part(*mesh_node->find_mesh_part("bnd:0"));
    unit_asm.add_mesh_part(*mesh_node->find_mesh_part("bnd:1"));
    unit_asm.add_mesh_part(*mesh_node->find_mesh_part("bnd:2"));


    // finally, assemble the filters
    Analytic::StaticWrapperFunction<2, VeloFuncX> inflow_func;
    unit_asm.assemble(filter.at<0>().at<0>(), space_v, inflow_func);
    unit_asm.assemble(filter.at<0>().at<1>(), space_v);
  }

  // filter matrices
  filter.at<0>().at<0>().filter_mat(mat_a.at<0,0>());
  filter.at<0>().at<1>().filter_mat(mat_a.at<1,1>());

  filter.at<0>().at<0>().filter_offdiag_row_mat(mat_a.at<0,1>());
  filter.at<0>().at<1>().filter_offdiag_row_mat(mat_a.at<1,0>());

  filter.at<0>().at<0>().filter_offdiag_row_mat(mat_b.at<0,0>());
  filter.at<0>().at<1>().filter_offdiag_row_mat(mat_b.at<1,0>());

  // create RHS and SOL
  VectorType vec_rhs = matrix.create_vector_l();
  VectorType vec_sol = matrix.create_vector_l();

  vec_rhs.format();
  vec_sol.format();

  filter.filter_rhs(vec_rhs);
  filter.filter_sol(vec_sol);

  // select vanka type
  //Solver::VankaType vtype = Solver::VankaType::nodal_diag;
  //Solver::VankaType vtype = Solver::VankaType::nodal_full;
  //Solver::VankaType vtype = Solver::VankaType::block_diag;
  //Solver::VankaType vtype = Solver::VankaType::block_full;

  // create vanka
  //auto vanka = std::make_shared<Solver::Vanka<MatrixType, FilterType>>(matrix, filter, vtype, 1.0, 100);
  auto vanka = Solver::new_vanka(matrix, filter, vtype, 1.0, 10);

  // create richardson
  auto solver = Solver::new_richardson(matrix, filter, 1.0, vanka);
  //auto solver = Solver::new_fgmres(matrix, filter, 16, 0.0, vanka);

  solver->set_max_iter(1);
  //solver->set_plot(true);

  TimeStamp stamp1, stamp2;

  stamp1.stamp();
  solver->init_symbolic();
  stamp2.stamp();
  String symbo_time = stamp2.elapsed_string(stamp1, TimeFormat::s_m);

  stamp1.stamp();
  solver->init_numeric();
  stamp2.stamp();
  String numer_time = stamp2.elapsed_string(stamp1, TimeFormat::s_m);


  // solve
  stamp1.stamp();
  /*Solver::Status status =*/ Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);
  stamp2.stamp();
  String solve_time = stamp2.elapsed_string(stamp1, TimeFormat::s_m);

  /*std::cout << std::endl
            << "Solver Status: " << status << std::endl;
  std::cout << "Convergence..: " << solver->get_conv_rate() << std::endl;
  std::cout << "Elapsed Time.: " << stamp2.elapsed_string(stamp1, TimeFormat::s_m) << std::endl;*/
  std::cout << name << ": ";
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
  std::cout
    << stringify(solver->calc_convergence_rate()).pad_front(17) << " "
    << symbo_time << " "
    << numer_time << " "
    << solve_time << " "
    << stringify(vanka->data_size()).pad_front(10)
    << stringify(matrix.used_elements()).pad_front(10)
    << std::endl;
}

int main()
{
  int lvl = 2;
  //*
  test_poiseuille<SpaceRTQ0>(lvl, false, Solver::VankaType::nodal_diag_mult);
  test_poiseuille<SpaceRTQ0>(lvl, false, Solver::VankaType::block_full_mult);
  test_poiseuille<SpaceRTQ0>(lvl, false, Solver::VankaType::nodal_diag_add);
  test_poiseuille<SpaceRTQ0>(lvl, false, Solver::VankaType::block_full_add);
  test_poiseuille<SpaceRTQ0>(lvl, true , Solver::VankaType::nodal_diag_mult);
  test_poiseuille<SpaceRTQ0>(lvl, true , Solver::VankaType::block_full_mult);
  test_poiseuille<SpaceRTQ0>(lvl, true , Solver::VankaType::nodal_diag_add);
  test_poiseuille<SpaceRTQ0>(lvl, true , Solver::VankaType::block_full_add);
  //*/
  //*
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::nodal_diag_mult);
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::nodal_full_mult);
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::block_diag_mult);
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::block_full_mult);
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::nodal_diag_add);
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::nodal_full_add);
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::block_diag_add);
  test_poiseuille<SpaceQ2P1>(lvl, false, Solver::VankaType::block_full_add);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::nodal_diag_mult);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::nodal_full_mult);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::block_diag_mult);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::block_full_mult);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::nodal_diag_add);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::nodal_full_add);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::block_diag_add);
  test_poiseuille<SpaceQ2P1>(lvl, true , Solver::VankaType::block_full_add);
  //*/
  //*
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::nodal_diag_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::nodal_full_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::block_diag_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::block_full_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::nodal_diag_add);
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::nodal_full_add);
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::block_diag_add);
  test_poiseuille<SpaceQ2Q1>(lvl, false, Solver::VankaType::block_full_add);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::nodal_diag_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::nodal_full_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::block_diag_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::block_full_mult);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::nodal_diag_add);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::nodal_full_add);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::block_diag_add);
  test_poiseuille<SpaceQ2Q1>(lvl, true , Solver::VankaType::block_full_add);
  //*/
}

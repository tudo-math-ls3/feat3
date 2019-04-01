// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
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
#include <kernel/solver/amavanka.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/gpdv_assembler.hpp>
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

template<typename DT_, typename IT_>
class AmaVankaTest :
  public FullTaggedTest<Mem::Main, DT_, IT_>
{
public:
  AmaVankaTest() :
    FullTaggedTest<Mem::Main, DT_, IT_>("AmaVankaTest")
  {
  }

  virtual ~AmaVankaTest()
  {
  }

  static constexpr int dim = 2;
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNode;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  /// tests Vanka for Deformation tensor formulation
  template<template<typename> class Space_>
  void test_defo(MeshNode& mesh_node, const DT_ omega) const
  {
    // get our mesh
    MeshType& mesh = *mesh_node.get_mesh();

    // create our trafo
    TrafoType trafo(mesh);

    typename Space_<TrafoType>::V space_v(trafo);
    typename Space_<TrafoType>::P space_p(trafo);
    const String name = Space_<TrafoType>::name();

    typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim, dim> MatrixTypeA;
    typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim, 1> MatrixTypeB;
    typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, 1, dim> MatrixTypeD;
    typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;

    typedef typename MatrixType::VectorTypeL VectorType;

    typedef LAFEM::UnitFilterBlocked<Mem::Main, DT_, IT_, dim> VeloFilterType;
    typedef LAFEM::NoneFilter<Mem::Main, double, Index> PresFilterType;
    typedef LAFEM::TupleFilter<VeloFilterType, PresFilterType> FilterType;

    MatrixType matrix;
    MatrixTypeA& mat_a = matrix.block_a();
    MatrixTypeB& mat_b = matrix.block_b();
    MatrixTypeD& mat_d = matrix.block_d();
    FilterType filter;

    Cubature::DynamicFactory cubature("auto-degree:5");

    matrix.format();

    // assemble A, B and D
    Assembly::SymbolicAssembler::assemble_matrix_std1(mat_a, space_v);

    Assembly::Common::DuDvOperatorBlocked<2> dudv_op;
    Assembly::BilinearOperatorAssembler::assemble_block_matrix1(mat_a, dudv_op, space_v, cubature);
    Assembly::GradPresDivVeloAssembler::assemble(mat_b, mat_d, space_v, space_p, cubature);

    // assemble filter
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:0"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:1"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:2"));

    // finally, assemble the filters
    Analytic::Common::ParProfileVector inflow_func(0.0, 0.0, 0.0, 1.0, 1.0);
    unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);

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

    // compute initial defect
    matrix.apply(vec_def, vec_sol, vec_rhs, -DT_(1));
    const DT_ def0 = vec_def.norm2();

    TimeStamp stamp1, stamp2;
    {
      // create vanka
      auto vanka = Solver::new_amavanka(matrix, filter, omega, 10);

      // initialise and solve
      vanka->init();
      stamp1.stamp();
      Solver::solve(*vanka, vec_sol, vec_rhs, matrix, filter);
      stamp2.stamp();
    }

    // compute final defect
    matrix.apply(vec_def, vec_sol, vec_rhs, -DT_(1));
    const DT_ def1 = vec_def.norm2();

    // print some information
    std::cout << name << ": ";
    std::cout << stringify_fp_sci(def0) << " > " << stringify_fp_sci(def1);
    std::cout << " : " << std::fixed << std::setprecision(5) << (def1/def0);
    std::cout << " | " << stamp2.elapsed_string(stamp1, TimeFormat::s_m);
    std::cout << std::endl;

    // ensure that the defect decreased
    TEST_CHECK_IN_RANGE(def1/def0, DT_(0), DT_(0.5));
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
    test_defo<RTQ0>(*mesh_node, DT_(1.0));

    // test Q2/P1dc
    test_defo<Q2P1>(*mesh_node, DT_(1.0));
  }
};

AmaVankaTest<double, Index> amavanka_test_double_index;

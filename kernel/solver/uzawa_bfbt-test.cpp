// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>

#include <kernel/analytic/common.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/common_functionals.hpp> // for ForceFunctional
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/discrete_projector.hpp> // for DiscreteVertexProjector
#include <kernel/assembly/domain_assembler.hpp>   // for DomainAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp> // for Assembly::assemble_***
#include <kernel/assembly/gpdv_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>  // NEW: for GridTransfer
#include <kernel/assembly/mean_filter_assembler.hpp> // for MeanFilterAssembler
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/velocity_analyser.hpp> // NEW: for VelocityAnalyser
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/unit_cube_patch_generator.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/bfbt.hpp>
#include <kernel/solver/jacobi_precond.hpp> // NEW: for JacobiPrecond
#include <kernel/solver/richardson.hpp>     // NEW: for Richardson
#include <kernel/solver/uzawa_precond.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/time_stamp.hpp>

#include <iomanip>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

template <typename T_>
struct VeloFuncX
{
  static T_ eval(T_, T_ y) { return y * (T_(1) - y); }
};

template <typename Trafo_>
struct RTQ0
{
  typedef Space::CroRavRanTur::Element<Trafo_> V;
  typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<0>> P;
  static const char *name() { return "RT/Q0"; }
};

template <typename Trafo_>
struct Q2P1
{
  typedef Space::Lagrange2::Element<Trafo_> V;
  typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<1>> P;
  static const char *name() { return "Q2/P1"; }
};

template <typename DT_, typename IT_>
class UzawaBFBTTest :
  public UnitTest
{
public:
  UzawaBFBTTest(PreferredBackend backend) :
    UnitTest("UzawaBFBTTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~UzawaBFBTTest()
  {
  }

  static constexpr int dim = 2;
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNode;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  /// tests Uzawa-BFBT for Deformation tensor formulation
  template <template <typename> class Space_>
  void test_defo(MeshNode& mesh_node, UzawaType uzawa_type) const
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
    typedef LAFEM::DenseVectorBlocked<DT_, IT_, dim> VectorVeloType;
    typedef LAFEM::DenseVector<DT_, IT_> VectorPresType;
    MatrixType matrix;
    MatrixTypeA& mat_a = matrix.block_a();
    MatrixTypeB& mat_b = matrix.block_b();
    MatrixTypeD& mat_d = matrix.block_d();
    FilterType filter;

    Cubature::DynamicFactory cubature("auto-degree:5");
    String cubature_name = "auto-degree:5";

    // assemble A, B and D
    Assembly::SymbolicAssembler::assemble_matrix_std1(mat_a, space_v);
    mat_a.format();

    Assembly::Common::DuDvOperatorBlocked<2> dudv_op;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a, dudv_op, space_v, cubature);
    Assembly::GradPresDivVeloAssembler::assemble(mat_b, mat_d, space_v, space_p, cubature);

    // assemble filter
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:0"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:1"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:2"));

    // finally, assemble the filters
    Analytic::Common::ParProfileVector<DT_> inflow_func(0.0, 0.0, 0.0, 1.0, 1.0);
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

    // create and compile domain assembler
    Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
    domain_assembler.compile_all_elements();

    // assemble velocity mass matrix
    MatrixTypeA M_v;
    Assembly::SymbolicAssembler::assemble_matrix_std1(M_v, space_v);
    M_v.format();
    Assembly::Common::IdentityOperatorBlocked<dim> identity_op_blocked;
    Assembly::assemble_bilinear_operator_matrix_1(domain_assembler, M_v, identity_op_blocked, space_v, cubature_name);

    // assemble lumped velocity vector
    VectorVeloType lumped_velo_mass_vec;
    lumped_velo_mass_vec = VectorVeloType(M_v.rows());
    M_v.lump_rows(lumped_velo_mass_vec);

    // invert lumped velocity vector and filter it
    lumped_velo_mass_vec.component_invert(lumped_velo_mass_vec);
    filter.template at<0>().filter_def(lumped_velo_mass_vec);

    // assemble structure of L_p
    LAFEM::SparseMatrixCSR<DT_, IT_> L_p;

    FEAT::Adjacency::Graph b_d(Adjacency::RenderType::injectify_sorted, mat_d, mat_b);
    L_p.convert(b_d);
    L_p.format();

    // fill L_p
    L_p.add_double_mat_product(mat_d, lumped_velo_mass_vec, mat_b);


    // solver for matrix A
    auto jacobi_A = Solver::new_jacobi_precond(mat_a, filter.template at<0>(), 0.5);
    auto solver_A = Solver::new_richardson(mat_a, filter.template at<0>(), 1.0, jacobi_A);

    auto jacobi_L_p = Solver::new_jacobi_precond(L_p, filter.template at<1>(), 0.5);
    auto rich = Solver::new_richardson(L_p, filter.template at<1>(), 1.0, jacobi_L_p);
    rich->set_max_iter(4);
    std::shared_ptr<Solver::SolverBase<VectorPresType>> solver_L_p_left;
    std::shared_ptr<Solver::SolverBase<VectorPresType>> solver_L_p_right;

    solver_L_p_left   = rich;
    solver_L_p_right   = rich;

    std::cout << "Testing solver pair: " << "rich/rich" << "\n";

    TimeStamp stamp1, stamp2;
    {
      // create smoother and solver
      auto bfbt = Solver::new_bfbt(mat_a, mat_b, mat_d, filter.template at<0>(), filter.template at<1>(), solver_L_p_left, solver_L_p_right, lumped_velo_mass_vec);
      auto uzawa = Solver::new_uzawa_precond(mat_a, mat_b, mat_d, filter.template at<0>(), filter.template at<1>(), solver_A, bfbt, uzawa_type);
      auto richardson = Solver::new_richardson(matrix, filter, 0.5, uzawa);

      // initialize and solve
      richardson->init();
      stamp1.stamp();
      Solver::solve(*richardson, vec_sol, vec_rhs, matrix, filter);
      stamp2.stamp();
    }

    // compute final defect
    matrix.apply(vec_def, vec_sol, vec_rhs, -DT_(1));
    const DT_ def1 = vec_def.norm2();

    // print some information
    std::cout << name << " ";
    std::cout << stringify_fp_sci(def0) << " > " << stringify_fp_sci(def1);
    std::cout << " : " << stringify_fp_fix(def1 / def0, 5);
    std::cout << " | " << stamp2.elapsed_string(stamp1, TimeFormat::s_m);
    std::cout << '\n';

    // ensure that the defect decreased
    TEST_CHECK_IN_RANGE(def1 / def0, DT_(0), DT_(0.5));

  }

  virtual void run() const override
  {
    const int level = 2;

    // create mesh node
    std::vector<int> ranks;
    std::unique_ptr<MeshNode> mesh_node;
    Geometry::UnitCubePatchGenerator<MeshType>::create_unique(0, 1, mesh_node, ranks);

    // refine a few times
    for(int i = 0; i < level; ++i) {
      mesh_node = mesh_node->refine_unique();
    }

    // test RT/Q0 and Q2/P1d for UzawaType upper
    UzawaType type_upper = Solver::UzawaType::upper;
    test_defo<RTQ0>(*mesh_node, type_upper);
    test_defo<Q2P1>(*mesh_node, type_upper);
    std::cout << "Upper passed\n\n";

    // test RT/Q0 and Q2/P1d for UzawaType lower
    UzawaType type_lower = Solver::UzawaType::lower;
    test_defo<RTQ0>(*mesh_node, type_lower);
    test_defo<Q2P1>(*mesh_node, type_lower);
    std::cout << "Lower passed\n\n";

    // test RT/Q0 and Q2/P1d for UzawaType full
    UzawaType type_full = Solver::UzawaType::full;
    test_defo<RTQ0>(*mesh_node, type_full);
    test_defo<Q2P1>(*mesh_node, type_full);
    std::cout << "Full passed\n\n";

    // test RT/Q0 and Q2/P1d for UzawaType diagonal
    UzawaType type_diagonal = Solver::UzawaType::diagonal;
    test_defo<RTQ0>(*mesh_node, type_diagonal);
    test_defo<Q2P1>(*mesh_node, type_diagonal);
    std::cout << "Diagonal passed\n\n";
  }
};

// UzawaBFBTTest<float, unsigned int>
// bfbf_test_float_uint(PreferredBackend::generic); UzawaBFBTTest<double, unsigned
// int> bfbf_test_double_uint(PreferredBackend::generic); UzawaBFBTTest<float,
// unsigned long> bfbf_test_float_ulong(PreferredBackend::generic);
UzawaBFBTTest<double, std::uint64_t> bfbf_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
UzawaBFBTTest<float, std::uint64_t> mkl_bfbf_test_float_uint64(PreferredBackend::mkl);
UzawaBFBTTest<double, std::uint64_t> mkl_bfbf_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
UzawaBFBTTest<__float128, std::uint32_t> bfbf_test_float128_uint32(PreferredBackend::generic);
UzawaBFBTTest<__float128, std::uint64_t> bfbf_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
UzawaBFBTTest<Half, std::uint32_t> bfbf_test_half_uint32(PreferredBackend::generic);
UzawaBFBTTest<Half, std::uint64_t> bfbf_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
UzawaBFBTTest<float, std::uint32_t> cuda_bfbf_test_float_uint32(PreferredBackend::cuda);
UzawaBFBTTest<double, std::uint32_t> cuda_bfbf_test_double_uint32(PreferredBackend::cuda);
UzawaBFBTTest<float, std::uint64_t> cuda_bfbf_test_float_uint64(PreferredBackend::cuda);
UzawaBFBTTest<double, std::uint64_t> cuda_bfbf_test_double_uint64(PreferredBackend::cuda);
#endif

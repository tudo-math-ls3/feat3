// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/unit_cube_patch_generator.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/amavanka_richardson.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/gpdv_assembler.hpp>
#include <kernel/analytic/common.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

/**
 * \brief Validates the fused CUDA AmaVankaRichardson smoother against its generic (CPU) reference.
 *
 * The generic path of AmaVankaRichardson is literally the Richardson(Schwarz(AmaVanka)) loop it is
 * meant to replace; on a pure-Dirichlet system (velocity UnitFilterBlocked + pressure NoneFilter)
 * the phase-1 unit-only device filtering and the full filtering coincide, so both paths must agree.
 *
 * \author Maximilian Esser
 */
template<typename DT_, typename IT_>
class AmaVankaRichardsonTest :
  public UnitTest
{
  static constexpr int dim = 2;
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::RootMeshNode<MeshType> MeshNode;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

public:
  explicit AmaVankaRichardsonTest(PreferredBackend backend) :
    UnitTest("AmaVankaRichardsonTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    // build a small Q2/P1dc unit-cube system
    const int level = 2;
    std::vector<int> ranks;
    std::unique_ptr<MeshNode> mesh_node;
    Geometry::UnitCubePatchGenerator<MeshType>::create_unique(0, 1, mesh_node, ranks);
    for(int i = 0; i < level; ++i)
      mesh_node = mesh_node->refine_unique();

    MeshType& mesh = *mesh_node->get_mesh();
    TrafoType trafo(mesh);
    Space::Lagrange2::Element<TrafoType> space_v(trafo);
    Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> space_p(trafo);

    typedef SparseMatrixBCSR<DT_, IT_, dim, dim> MatrixTypeA;
    typedef SparseMatrixBCSR<DT_, IT_, dim, 1> MatrixTypeB;
    typedef SparseMatrixBCSR<DT_, IT_, 1, dim> MatrixTypeD;
    typedef SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> LocalMatrix;
    typedef typename LocalMatrix::VectorTypeL LocalVector;

    typedef UnitFilterBlocked<DT_, IT_, dim> VeloFilter;
    typedef NoneFilter<DT_, IT_> PresFilter;
    typedef TupleFilter<VeloFilter, PresFilter> LocalFilter;

    typedef VectorMirror<DT_, IT_> ScalarMirror;
    typedef TupleMirror<ScalarMirror, ScalarMirror> SysMirror;

    typedef Global::Gate<LocalVector, SysMirror> GateType;
    typedef Global::Matrix<LocalMatrix, SysMirror, SysMirror> GlobalMatrix;
    typedef Global::Filter<LocalFilter, SysMirror> GlobalFilter;
    typedef Global::Vector<LocalVector, SysMirror> GlobalVector;

    // assemble the local system
    LocalMatrix matrix;
    LocalFilter filter;
    Cubature::DynamicFactory cubature("auto-degree:5");
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix.block_a(), space_v);
    matrix.block_a().format();
    Assembly::Common::DuDvOperatorBlocked<dim> dudv_op;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix.block_a(), dudv_op, space_v, cubature);
    Assembly::GradPresDivVeloAssembler::assemble(matrix.block_b(), matrix.block_d(), space_v, space_p, cubature);

    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(*mesh_node->find_mesh_part("bnd:0"));
    unit_asm.add_mesh_part(*mesh_node->find_mesh_part("bnd:1"));
    unit_asm.add_mesh_part(*mesh_node->find_mesh_part("bnd:2"));
    unit_asm.add_mesh_part(*mesh_node->find_mesh_part("bnd:3"));
    Analytic::Common::ParProfileVector<DT_> inflow_func(0.0, 0.0, 0.0, 1.0, 1.0);
    unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);
    filter.template at<0>().filter_mat(matrix.block_a());
    filter.template at<0>().filter_offdiag_row_mat(matrix.block_b());

    // wrap into single-process global objects with a trivial (empty) gate
    Dist::Comm comm(Dist::Comm::world());
    GateType gate;
    gate.set_comm(&comm);
    {
      LocalVector freq = matrix.create_vector_l();
      freq.format(DT_(1));
      gate.compile(std::move(freq));
    }
    GlobalMatrix gmat(&gate, &gate, std::move(matrix));
    GlobalFilter gfilter(std::move(filter));

    // rhs / initial solution
    GlobalVector vec_rhs = gmat.create_vector_l();
    vec_rhs.format(DT_(1));
    gfilter.filter_rhs(vec_rhs);
    GlobalVector vec_ref = gmat.create_vector_l();
    vec_ref.format();
    gfilter.filter_sol(vec_ref);
    GlobalVector vec_cuda = vec_ref.clone(LAFEM::CloneMode::Deep);

    // build the smoother (fixed 3 Richardson steps, omega = 0.7)
    auto vanka = Solver::new_amavanka(gmat.local(), gfilter.local(), DT_(1), Index(1));
    vanka->set_skip_singular(true);
    auto rich = Solver::new_amavanka_richardson(gmat, gfilter, vanka, DT_(0.7));
    rich->set_min_iter(Index(3));
    rich->set_max_iter(Index(3));

    // assemble the vanka on the host so both paths start from identical numeric data
    Backend::set_preferred_backend(PreferredBackend::generic);
    rich->init();

    // reference: generic (CPU) Richardson(Schwarz(AmaVanka)) loop
    rich->set_use_cuda(false);
    rich->correct(vec_ref, vec_rhs);

    // candidate: fused CUDA driver
    #ifdef FEAT_HAVE_CUDA
    Backend::set_preferred_backend(PreferredBackend::cuda);
    rich->set_use_cuda(true);
    rich->correct(vec_cuda, vec_rhs);
    #endif

    rich->done();
    Backend::set_preferred_backend(PreferredBackend::generic);

    // compare
    GlobalVector diff = vec_ref.clone(LAFEM::CloneMode::Deep);
    diff.axpy(vec_cuda, -DT_(1));
    const DT_ dn = diff.norm2();
    const DT_ rn = vec_ref.norm2();
    const DT_ tol = (sizeof(DT_) >= 8u) ? DT_(1e-8) : DT_(1e-3);
    std::cout << "AmaVankaRichardson cuda-vs-ref: |diff|=" << stringify_fp_sci(dn)
              << " |ref|=" << stringify_fp_sci(rn) << "\n";
    TEST_CHECK(dn <= tol * (rn + DT_(1)));
  }
};

#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(AmaVankaRichardsonTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(AmaVankaRichardsonTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

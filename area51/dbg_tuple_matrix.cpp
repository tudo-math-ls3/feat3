// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>         // for RefinedUnitCubeFactor
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping
#include <kernel/space/lagrange3/element.hpp>              // the Lagrange-2 Element (aka "Q2")
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-2 Element (aka "Q2")
#include <kernel/space/discontinuous/element.hpp>          // the Discontinuous Element (aka "P1dc")
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/analytic/common.hpp>                      // for ParProfileFunction
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/mean_filter_assembler.hpp>       // for MeanFilterAssembler
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/gpdv_assembler.hpp>              // for GradPresDivVeloAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/velocity_analyser.hpp>           // NEW: for VelocityAnalyser
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/dense_vector_blocked.hpp>           // NEW: for DenseVectorBlocked
#include <kernel/lafem/tuple_vector.hpp>                   // NEW: for TupleVector
#include <kernel/lafem/sparse_matrix_bcsr.hpp>             // NEW: for SparseMatrixBCSR
#include <kernel/lafem/saddle_point_matrix.hpp>            // NEW: for SaddlePointMatrix
#include <kernel/lafem/null_matrix.hpp>
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/mean_filter.hpp>                    // NEW: for MeanFilter
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>                   // NEW: for TupleFilter
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/solver/bicgstab.hpp>

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT, so use the namespace here.
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial07
{
  // Once again, we use quadrilaterals.
  typedef Shape::Quadrilateral ShapeType;
  // Use the unstructured conformal mesh class
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  // Define the corresponding mesh-part type
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  // Use the standard transformation mapping
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // Use the Lagrange-2 element for velocity
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
  //typedef Space::Lagrange3::Element<TrafoType> SpaceStressType;
  typedef Space::Lagrange3::Element<TrafoType> SpaceStressType;

  // Use the P1dc element for pressure
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;

  // Finally, store the dimension of the shape as a constant here for convenience,
  // as we will require this several times in a moment:
  static constexpr int dim = ShapeType::dimension;

  class StokesHitTestFunction
  {
  public:
    template<typename Point_>
    bool operator()(const Point_& p) const
    {
      return (p[0] < 1E-5) || (p[1] < 1E-5) || (p[1] > 1-1E-5);
    }
  };

  class MyOperatorR :
    public Assembly::BilinearOperator
  {
  public:
    static constexpr int BlockHeight = 2;
    static constexpr int BlockWidth = 4;

    static constexpr TrafoTags trafo_config = TrafoTags::none;
    static constexpr SpaceTags test_config = SpaceTags::value;
    static constexpr SpaceTags trial_config = SpaceTags::grad;

    template<typename AsmTraits_>
    class Evaluator :
      public BilinearOperator::Evaluator<AsmTraits_>
    {
    public:
      /// the data type to be used
      typedef typename AsmTraits_::DataType DataType;
      /// the data type for the block system
      typedef typename AsmTraits_::OperatorValueType OperatorValueType;
      /// the assembler's trafo data type
      typedef typename AsmTraits_::TrafoData TrafoData;
      /// the assembler's test-function data type
      typedef typename AsmTraits_::TestBasisData TestBasisData;
      /// the assembler's trial-function data type
      typedef typename AsmTraits_::TrialBasisData TrialBasisData;

    public:
      explicit Evaluator(const MyOperatorR&)
      {
      }

      OperatorValueType operator()(const TrialBasisData& phi, const TestBasisData& psi)
      {
        OperatorValueType r(DataType(0));
        r(0,0) = r(1,2) = -phi.grad[0] * psi.value;
        r(0,1) = r(1,3) = -phi.grad[1] * psi.value;
        return r;
      }
    };
  };

  class MyOperatorK :
    public Assembly::BilinearOperator
  {
  public:
    static constexpr int BlockHeight = 4;
    static constexpr int BlockWidth = 2;

    static constexpr TrafoTags trafo_config = TrafoTags::none;
    static constexpr SpaceTags test_config = SpaceTags::value;
    static constexpr SpaceTags trial_config = SpaceTags::grad;

    template<typename AsmTraits_>
    class Evaluator :
      public BilinearOperator::Evaluator<AsmTraits_>
    {
    public:
      /// the data type to be used
      typedef typename AsmTraits_::DataType DataType;
      /// the data type for the block system
      typedef typename AsmTraits_::OperatorValueType OperatorValueType;
      /// the assembler's trafo data type
      typedef typename AsmTraits_::TrafoData TrafoData;
      /// the assembler's test-function data type
      typedef typename AsmTraits_::TestBasisData TestBasisData;
      /// the assembler's trial-function data type
      typedef typename AsmTraits_::TrialBasisData TrialBasisData;

    public:
      explicit Evaluator(const MyOperatorK&)
      {
      }

      OperatorValueType operator()(const TrialBasisData& phi, const TestBasisData& psi)
      {
        OperatorValueType r(DataType(0));
        r(0,0) = -phi.grad[0] * psi.value;
        r(3,1) = -phi.grad[1] * psi.value;
        r(1,0) = r(2,0) = -0.5 * phi.grad[1] * psi.value;
        r(1,1) = r(2,1) = -0.5 * phi.grad[0] * psi.value;
        return r;
      }
    };
  };

  template<typename DataType_, typename IndexType_, typename SpaceVelo_, typename SpaceStress_>
  void assemble_stress_matrix(
    LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, 4, 4>& matrix,
    const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, 2>& convect,
    const SpaceVelo_& space_velo,
    const SpaceStress_& space_stress,
    const Cubature::DynamicFactory& cubature_factory,
    const DataType_ alpha = DataType_(1),
    const DataType_ lambda = DataType_(0)
    )
  {
    // validate matrix and vector dimensions
    XASSERTM(matrix.rows() == space_stress.get_num_dofs(), "invalid matrix dimensions");
    XASSERTM(matrix.columns() == space_stress.get_num_dofs(), "invalid matrix dimensions");
    XASSERTM(convect.size() == space_velo.get_num_dofs(), "invalid vector size");

    typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, 2> VectorType;
    typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, 4, 4> MatrixType;

    // define our assembly traits
    typedef Assembly::AsmTraits2<
      DataType_,
      SpaceStress_,
      SpaceVelo_,
      TrafoTags::jac_det,
      SpaceTags::value|SpaceTags::grad,
      SpaceTags::value|SpaceTags::grad
    > AsmTraits;

    // fetch our trafo
    const typename AsmTraits::TrafoType& trafo = space_stress.get_trafo();

    // create a trafo evaluator
    typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

    // create a space evaluator and evaluation data
    typename AsmTraits::TestEvaluator space_eval_stress(space_stress);
    typename AsmTraits::TrialEvaluator space_eval_velo(space_velo);

    // create a dof-mapping
    typename AsmTraits::TestDofMapping dof_mapping_stress(space_stress);
    typename AsmTraits::TrialDofMapping dof_mapping_velo(space_velo);

    // create trafo evaluation data
    typename AsmTraits::TrafoEvalData trafo_data;

    // create space evaluation data
    typename AsmTraits::TestEvalData space_data_stress;
    typename AsmTraits::TrialEvalData space_data_velo;

    // create cubature rule
    typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

    // create matrix scatter-axpy
    typename MatrixType::ScatterAxpy scatter_matrix(matrix);

    // create convection gather-axpy
    typename VectorType::GatherAxpy gather_conv(convect);

    // get maximum number of local dofs
    static constexpr int max_local_stress_dofs = AsmTraits::max_local_test_dofs;
    static constexpr int max_local_velo_dofs = AsmTraits::max_local_trial_dofs;

    // create local matrix data
    typedef Tiny::Matrix<DataType_, 4, 4> MatrixValue;
    typedef Tiny::Matrix<MatrixValue, max_local_stress_dofs, max_local_stress_dofs> LocalMatrixType;
    LocalMatrixType local_matrix;

    // create local vector data
    typedef Tiny::Vector<DataType_, 2> VectorValue;
    typedef Tiny::Vector<VectorValue, max_local_velo_dofs> LocalVectorType;

    // local convection field dofs
    LocalVectorType local_conv_dofs;

    // our local velocity value
    Tiny::Vector<DataType_, 2> loc_v;

    // our local velocity gradient
    Tiny::Matrix<DataType_, 2, 2> loc_grad_v;

    loc_v.format();
    loc_grad_v.format();

    // loop over all cells of the mesh
    for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
    {
      // prepare trafo evaluator
      trafo_eval.prepare(cell);

      // prepare space evaluator
      space_eval_stress.prepare(trafo_eval);
      space_eval_velo.prepare(trafo_eval);

      // initialize dof-mapping
      dof_mapping_stress.prepare(cell);
      dof_mapping_velo.prepare(cell);

      // fetch number of local dofs
      const int num_loc_stress_dofs = space_eval_stress.get_num_local_dofs();
      const int num_loc_velo_dofs = space_eval_velo.get_num_local_dofs();

      // gather our local convection dofs
      local_conv_dofs.format();
      gather_conv(local_conv_dofs, dof_mapping_velo);

      // format our local matrix and vector
      local_matrix.format();

      // loop over all quadrature points and integrate
      for(int point(0); point < cubature_rule.get_num_points(); ++point)
      {
        // compute trafo data
        trafo_eval(trafo_data, cubature_rule.get_point(point));

        // compute basis function data
        space_eval_stress(space_data_stress, trafo_data);
        space_eval_velo(space_data_velo, trafo_data);

        // pre-compute cubature weight
        const DataType_ weight = trafo_data.jac_det * cubature_rule.get_weight(point);

        // evaluate convection function and its gradient (if required)
        loc_v.format();
        for(int i(0); i < num_loc_velo_dofs; ++i)
        {
          // update velocity value
          loc_v.axpy(space_data_velo.phi[i].value, local_conv_dofs[i]);
        }
        loc_grad_v.format();
        for(int i(0); i < num_loc_velo_dofs; ++i)
        {
          // update velocity gradient
          loc_grad_v.add_outer_product(local_conv_dofs[i], space_data_velo.phi[i].grad);
        }

        // test function loop
        for(int i(0); i < num_loc_stress_dofs; ++i)
        {
          // trial function loop
          for(int j(0); j < num_loc_stress_dofs; ++j)
          {
            const auto& psi = space_data_stress.phi[i];
            const auto& phi = space_data_stress.phi[j];

            // mass
            const DataType_ mass = weight * alpha * phi.value * psi.value;

            // for sigma_11
            local_matrix[i][j][0][0] += mass + weight * lambda * psi.value * (
              loc_v[0]*phi.grad[0] + loc_v[1] * phi.grad[1] // u1*dx + u2*dy
              + DataType_(2) * loc_grad_v(0,0) * phi.value // 2*dx(u1)
            );

            local_matrix[i][j][0][1] += weight * lambda * loc_grad_v(1,0) * phi.value * psi.value;
            local_matrix[i][j][0][2] += weight * lambda * loc_grad_v(1,0) * phi.value * psi.value;

            // for sigma_12
            local_matrix[i][j][1][1] += mass + weight * lambda * psi.value * (
              loc_v[0]*phi.grad[0] + loc_v[1] * phi.grad[1] // u1*dx + u2*dy
              + (loc_grad_v(0,0) + loc_grad_v(1,1)) * phi.value // dx(u1)+dy(u2)
            );

            local_matrix[i][j][1][0] += weight * lambda * loc_grad_v(0,1) * phi.value * psi.value;
            local_matrix[i][j][1][3] += weight * lambda * loc_grad_v(1,0) * phi.value * psi.value;

            // for sigma_12
            local_matrix[i][j][2][2] += mass + weight * lambda * psi.value * (
              loc_v[0]*phi.grad[0] + loc_v[1] * phi.grad[1] // u1*dx + u2*dy
              + (loc_grad_v(0,0) + loc_grad_v(1,1)) * phi.value // dx(u1)+dy(u2)
            );

            local_matrix[i][j][2][0] += weight * lambda * loc_grad_v(0,1) * phi.value * psi.value;
            local_matrix[i][j][2][3] += weight * lambda * loc_grad_v(1,0) * phi.value * psi.value;

            // for sigma_22
            local_matrix[i][j][3][3] += mass + weight * lambda * psi.value * (
              loc_v[0]*phi.grad[0] + loc_v[1] * phi.grad[1] // u1*dx + u2*dy
              + DataType_(2) * loc_grad_v(1,1) * phi.value // 2*dy(u2)
            );

            local_matrix[i][j][3][1] += weight * lambda * loc_grad_v(0,1) * phi.value * psi.value;
            local_matrix[i][j][3][2] += weight * lambda * loc_grad_v(0,1) * phi.value * psi.value;
          }
        }

        // continue with next cubature point
      }

      // scatter into matrix
      scatter_matrix(local_matrix, dof_mapping_stress, dof_mapping_stress, 1.0);

      // finish dof mapping
      dof_mapping_velo.finish();
      dof_mapping_stress.finish();

      // finish evaluators
      space_eval_stress.finish();
      space_eval_velo.finish();
      trafo_eval.finish();
    }
  }

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Linear System type definitions

  // Our LAFEM containers work in main memory.
  typedef Mem::Main MemType;
  // Our data arrays should be double precision.
  typedef double DataType;
  // Use the default index type for indexing.
  typedef Index IndexType;

  typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> VectorVeloType;
  typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, 4> VectorStressType;
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorPresType;

  typedef LAFEM::TupleVector<VectorVeloType, VectorPresType, VectorStressType> SystemVectorType;


  //typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixTypeS;
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, dim> MatrixTypeA;
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim,   1> MatrixTypeB;
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType,   1, dim> MatrixTypeD;

  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 4, 2> MatrixTypeK;
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 2, 4> MatrixTypeR;
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 4, 4> MatrixTypeM;

  typedef LAFEM::NullMatrix<MemType, DataType, IndexType, 1, 1> Matrix0_1x1;
  typedef LAFEM::NullMatrix<MemType, DataType, IndexType, 4, 1> Matrix0_4x1;
  typedef LAFEM::NullMatrix<MemType, DataType, IndexType, 1, 4> Matrix0_1x4;
  //typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 4, 1> MatrixType4x1;
  //typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 1, 4> MatrixType1x4;

  typedef LAFEM::TupleMatrix
  <
    LAFEM::TupleMatrixRow<MatrixTypeA, MatrixTypeB, MatrixTypeR>,
    LAFEM::TupleMatrixRow<MatrixTypeD, Matrix0_1x1, Matrix0_1x4>,
    LAFEM::TupleMatrixRow<MatrixTypeK, Matrix0_4x1, MatrixTypeM>
  > SystemMatrixType;

  typedef LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, dim> FilterVeloType;
  typedef LAFEM::NoneFilter<MemType, DataType, IndexType> FilterPresType;
  typedef LAFEM::NoneFilterBlocked<MemType, DataType, IndexType, 4> FilterStressType;

  typedef LAFEM::TupleFilter<FilterVeloType, FilterPresType, FilterStressType> SystemFilterType;

  // Here's our tutorial's main function
  void main(Index level)
  {
    // As usual, create a mesh and the corresponding boundary mesh-part first.

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;

    // Create the mesh
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    // Create the boundary mesh-part
    StokesHitTestFunction stokes_hit_func;
    Geometry::HitTestFactory<StokesHitTestFunction, MeshType> boundary_factory(stokes_hit_func, mesh);
    MeshPartType boundary(boundary_factory);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Trafo and Finite Element Space initialization

    std::cout << "Creating Trafo and Spaces..." << std::endl;

    // Create the trafo
    TrafoType trafo(mesh);

    // In the case of Stokes, there are two finite element spaces involved here:
    SpaceVeloType space_velo(trafo); // velocity space
    SpacePresType space_pres(trafo); // pressure space
    SpaceStressType space_stress(trafo);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly
    std::cout << "Allocating matrix and vectors..." << std::endl;

    // As usual, we create first an empty matrix.
    SystemMatrixType matrix;

    MatrixTypeA& matrix_a = matrix.at<0,0>();
    MatrixTypeB& matrix_b = matrix.at<0,1>();
    MatrixTypeD& matrix_d = matrix.at<1,0>();

    MatrixTypeK& matrix_k = matrix.at<2,0>();
    MatrixTypeR& matrix_r = matrix.at<0,2>();
    MatrixTypeM& matrix_m = matrix.at<2,2>();

    Matrix0_1x1& matrix_0_1x1 = matrix.at<1,1>();
    Matrix0_1x4& matrix_0_1x4 = matrix.at<1,2>();
    Matrix0_4x1& matrix_0_4x1 = matrix.at<2,1>();


    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_a, space_velo);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_0_1x1, space_pres);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_m, space_stress);

    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_b, space_velo, space_pres);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_d, space_pres, space_velo);

    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_r, space_velo, space_stress);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_k, space_stress, space_velo);

    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_0_1x4, space_pres, space_stress);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_0_4x1, space_stress, space_pres);

    // That's it for the symbolic assembly of the matrix sparsity patterns.
    // In analogy to all previous tutorials, we can now use our matrix to create compatible
    // vectors for us:
    SystemVectorType vec_sol = matrix.create_vector_r();
    SystemVectorType vec_rhs = matrix.create_vector_l();

    // Finally, format the matrix and the vectors for the upcoming numerical assembly.
    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    // First of all, let's create a cubature factory:
    Cubature::DynamicFactory cubature_factory("auto-degree:7");

    // Now we need to assemble the three sub-matrices individually.

    // Let us first assemble the sub-matrix A. The velocity diffusion operator of the Stokes
    // equations can be either the "full" deformation tensor (aka "Du:Dv") or its simplified
    // form, the gradient tensor (aka Laplace). We pick the gradient tensor here, but you
    // may as well use the deformation tensor here:
    Assembly::Common::LaplaceOperatorBlocked<dim> operator_a; // gradient tensor
    //Assembly::Common::DuDvOperatorBlocked<dim>  operator_a; // deformation tensor

    // In analogy to all previous (scalar) tutorials, we now use the BilinearOperatorAssembler
    // to do the dirty work for us. The only notable difference is that we have to call the
    // "assemble_block_matrix1" function instead of its scalar counterpart "assemble_matrix1":
    Assembly::BilinearOperatorAssembler::assemble_block_matrix1(
      matrix_a, operator_a, space_velo, cubature_factory);

    /*Assembly::Common::IdentityOperatorBlocked<4> operator_m;
    Assembly::BilinearOperatorAssembler::assemble_block_matrix1(
      matrix_m, operator_m, space_stress, cubature_factory);*/
    assemble_stress_matrix(matrix_m, vec_sol.at<0>(), space_velo, space_stress, cubature_factory, 1.0, 1.0);

    MyOperatorR operator_r;
    MyOperatorK operator_k;

    Assembly::BilinearOperatorAssembler::assemble_block_matrix2(
      matrix_r, operator_r, space_velo, space_stress, cubature_factory);
    Assembly::BilinearOperatorAssembler::assemble_block_matrix2(
      matrix_k, operator_k, space_stress, space_velo, cubature_factory);

    // Next, we need to assemble the sub-matrices B and D for the pressure gradient and
    // velocity divergence. The good news is that this operator pair is required so often
    // that someone (me, actually) decided to write a custom assembler for this task.
    // So instead of setting up some operators, we simply use the "GradPresDivVeloAssembler"
    // for that:
    Assembly::GradPresDivVeloAssembler::assemble(
      matrix_b,  // pressure gradient sub-matrix B
      matrix_d,  // velocity divergence sub-matrix D
      space_velo,        // velocity space
      space_pres,        // pressure space
      cubature_factory   // cubature factory
    );

    // Note that the right-rand-side of our Stokes equations is zero, so we don't have to assemble
    // the contents of the right-rand-side vector either.

    std::cout << "Matrix Frobenius Norms:" << std::endl;
    std::cout << "A..: " << stringify_fp_sci(matrix_a.norm_frobenius()) << std::endl;
    std::cout << "B..: " << stringify_fp_sci(matrix_b.norm_frobenius()) << std::endl;
    std::cout << "D..: " << stringify_fp_sci(matrix_d.norm_frobenius()) << std::endl;
    //std::cout << "S..: " << stringify_fp_sci(matrix_s.norm_frobenius()) << std::endl;
    std::cout << "K..: " << stringify_fp_sci(matrix_k.norm_frobenius()) << std::endl;
    std::cout << "R..: " << stringify_fp_sci(matrix_r.norm_frobenius()) << std::endl;
    std::cout << "M..: " << stringify_fp_sci(matrix_m.norm_frobenius()) << std::endl;
    //std::cout << "1x4: " << stringify_fp_sci(matrix_1x4.norm_frobenius()) << std::endl;
    //std::cout << "4x1: " << stringify_fp_sci(matrix_4x1.norm_frobenius()) << std::endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << std::endl;

    // Now we have to assemble the boundary condition filters, which is also performed by
    // assembling the individual parts of the system filter. So, first of all, let us create
    // a new empty filter:
    SystemFilterType filter;

    // We have Dirichlet boundary conditions for the velocity, so we first require a unit-filter
    // assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // As usual, add the whole boundary to the assembler:
    unit_asm.add_mesh_part(boundary);

    // Now we need to chose the boundary condition function. We want to have a parabolic profile
    // flux, so we can use the pre-defined "ParProfileVector" function, which gives us a parabolic
    // profile flux orthogonal to the line segment (x0,y0)--(x1,y1) with maximum magnitude 'v-max':
    Analytic::Common::ParProfileVector<DataType> profile_function(
      0.0, // x0
      0.0, // y0
      0.0, // x1
      1.0, // y1
      1.0  // v-max
    );

    // We can now assemble the unit-filter for the velocity component of the system filter.
    // As the system filter is a TupleFilter, we can access its first (i.e. velocity) component
    // by calling the 'at' member function template and passing the desired index as a template
    // parameter:
    unit_asm.assemble(
      filter.at<0>(), // velocity unit-filter component of system filter
      space_velo,              // velocity space
      profile_function         // parabolic profile function
    );

    // That's it for the velocity filter, so let us continue with the pressure filter.
    // As mentioned above, we require a "MeanFilter" for the pressure, which can be assembled
    // by the MeanFilterAssembler. Note that there is no special setup required for this assembler,
    // so we can call the "assemble" function directly:
    /*Assembly::MeanFilterAssembler::assemble(
      filter.at<1>(),  // pressure mean-filter component of system filter
      space_pres,               // pressure space
      cubature_factory          // cubature factory
    );*/

    // Okay, now we have a mean-filter for the pressure and therefore the assembly of the
    // system filter is complete.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << std::endl;

    // As usual, apply the filter onto the initial solution and right-hand-side vectors:
    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    // Moreover, we also want to filter the system matrix, just as we did in all previous
    // tutorials (except for tutorial 06). Unfortunately, this is not as easy as previously,
    // because the TupleFilter has no "filter_mat" function -- and this is for a reason.

    // The problem is that there is no general rule on whether and, if so, how filters have
    // to be applied on matrices. This gets even more complicated when you are dealing with
    // systems of PDEs with different (possibly chained) filters for the various components.

    // In this tutorial, we are dealing with a unit-filter/mean-filter pair for a Stokes system.
    // In this case, we have to apply the 'filter_mat' function of the velocity unit-filter
    // into the "A" sub-matrix (in analogy to the previous tutorials):
    filter.at<0>().filter_mat(matrix_a);

    // Moreover, we have also have to apply the unit-filter onto the sub-matrix "B", this
    // time calling the "filter_offdiag_row_mat" function to zero our all rows of B, which
    // correspond to test function DOFs affected by Dirichlet boundary conditions:
    filter.at<0>().filter_offdiag_row_mat(matrix_b);
    filter.at<0>().filter_offdiag_row_mat(matrix_r);

    // The mean-filter component of the pressure does not provide a means to filter matrices,
    // thus the "D" sub-matrix remains unfiltered -- this is, by the way, the reason why the
    // commonly (i.e. in literature) used identity "D = B^T" is not really correct from a
    // practical point of view.

    // The Stokes is system is fully assembled now, so we can continue with the linear solver.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    std::cout << "Solving linear system..." << std::endl;

#ifdef FEAT_HAVE_UMFPACK

    auto solver = Solver::new_generic_umfpack(matrix);

    solver->init_symbolic();

    //for(int i = 0; i < 3; ++i)
    {
      //matrix_m.format();
      //assemble_stress_matrix(matrix_m, vec_sol.at<0>(), space_velo, space_stress, cubature_factory, 1.0, 0.1);

      solver->init_numeric();
      Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);
      solver->done_numeric();
    }

    solver->done_symbolic();

#endif

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Analyse velocity field

    // At this point, we may compute the L2/H1-errors of the velocity and pressure against
    // reference solutions as in the previous tutorials. However, we want to do something else
    // for a change here, so let us analyse the velocity field, i.e. compute various quantities
    // from our field: divergence, vorticity and its H0/H1-norms.

    std::cout << std::endl << "Performing velocity field analysis..." << std::endl;

    // The class that performs this analysis is the "VelocityAnalyser" and it returns a
    // "VelicityInfo" object of the appropriate datatype and dimension:
    Assembly::VelocityInfo<DataType, dim> velo_info = Assembly::VelocityAnalyser::compute(
      vec_sol.at<0>(), // the velocity field
      space_velo,               // the velocity space
      cubature_factory          // a cubature factory
    );

    // Just as the "ScalarErrorInfo", that is returned by the "ScalarErrorComputer", we can now
    // simply push the velo_info object to cout:
    std::cout << velo_info << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // First of all, build the filename string
    String vtk_name(String("./tutorial-07-stokes-lvl") + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // First, let's add the velocity field to our exporter. As the velocity is a vector field
    // and not a scalar function, we have to call the "add_vertex_vector" function instead of
    // "add_vertex_scalar". Note that we can directly pass the Q2 solution vector to the function
    // without projecting it explicitly to the vertices:
    exporter.add_vertex_vector("velocity", vec_sol.at<0>());

    {
      // split stress into 4 components
      const Index nverts = mesh.get_num_vertices();
      double* sigma_xx = new double[nverts];
      double* sigma_xy = new double[nverts];
      double* sigma_yx = new double[nverts];
      double* sigma_yy = new double[nverts];

      const auto* sigma = vec_sol.at<2>().elements();
      for(Index i(0); i < nverts; ++i)
      {
        sigma_xx[i] = sigma[i][0];
        sigma_xy[i] = sigma[i][1];
        sigma_yx[i] = sigma[i][2];
        sigma_yy[i] = sigma[i][3];
      }
      exporter.add_vertex_scalar("sigma_xx", sigma_xx);
      exporter.add_vertex_scalar("sigma_xy", sigma_xy);
      exporter.add_vertex_scalar("sigma_yx", sigma_yx);
      exporter.add_vertex_scalar("sigma_yy", sigma_yy);
      delete [] sigma_yy;
      delete [] sigma_yx;
      delete [] sigma_xy;
      delete [] sigma_xx;
    }

    // Our pressure is a scalar function defined in the P1dc space. Unfortunately, we cannot simply
    // pass the pressure component of the solution vector to the 'add_scalar_cell' function, but we
    // need to project it onto the cells before. This works quite similar to the "vertex projection"
    // that has been demonstrated in the previous tutorials:
    VectorPresType cell_pres;
    Assembly::DiscreteCellProjector::project(
      cell_pres,                // the vector that receives the cell-projection of the pressure
      vec_sol.at<1>(), // the pressure-vector to be projection
      space_pres,               // the pressure space
      cubature_factory          // a cubature factory for the projection.
    );

    // Now we can add the cell-projected pressure to our VTK exporter:
    exporter.add_cell_scalar("pressure", cell_pres.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  } // int main(...)
} // namespace Tutorial07

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialize the FEAT runtime environment:
  Runtime::initialize(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #07: Stokes" << std::endl;

  // Specify the desired mesh refinement level, defaulted to 3.
  // Note that FEAT uses its own "Index" type rather than a wild mixture of int, uint, long
  // and such.
  Index level(3);

  // Now let's see if we have command line parameters: This tutorial supports passing
  // the refinement level as a command line parameter, to investigate the behavior of the L2/H1
  // errors of the discrete solution.
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[argc-1]).parse(ilevel) || (ilevel < 1))
    {
      // Failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level." << std::endl;
      std::cerr << "Note: The last argument must be a positive integer." << std::endl;
      // Abort our runtime environment
      Runtime::abort();
    }
    // If parsing was successful, use the given information and notify the user
    level = Index(ilevel);
    std::cout << "Refinement level: " << level << std::endl;
  }
  else
  {
    // No command line parameter given, so inform the user that defaults are used
    std::cout << "Refinement level (default): " << level << std::endl;
  }

  // call the tutorial's main function
  Tutorial07::main(level);

  // And finally, finalize our runtime environment. This function returns the 'EXIT_SUCCESS' return code,
  // so we can simply return this as the result of our main function to indicate a successful run.
  return Runtime::finalize();
}

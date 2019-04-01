// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/string.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/gpdv_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/assembly/oldroyd_assembler.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/null_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/solver/vanka.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/richardson.hpp>

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

using namespace FEAT;

namespace DbgAmaVanka
{
  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  struct BoundaryHitFunction
  {
    template<typename P_>
    bool operator()(const P_& p) const
    {
      return (p[0] < 0.001) || (p[1] < 0.001) || (p[1] > 0.999);
    }
  };

  template<int dim, int nsc>
  class OperatorK :
    public Assembly::BilinearOperator
  {
  public:
    static constexpr int BlockHeight = nsc;
    static constexpr int BlockWidth = dim;

    static constexpr TrafoTags trafo_config = TrafoTags::none;
    static constexpr SpaceTags test_config = SpaceTags::value;
    static constexpr SpaceTags trial_config = SpaceTags::grad;

    template<typename AsmTraits_>
    class Evaluator :
      public BilinearOperator::Evaluator<AsmTraits_>
    {
    public:
      typedef typename AsmTraits_::OperatorValueType OperatorValueType;
      typedef typename AsmTraits_::TestBasisData TestBasisData;
      typedef typename AsmTraits_::TrialBasisData TrialBasisData;

    public:
      explicit Evaluator(const OperatorK& DOXY(operat)) {}

      // 2D version with 4 components
      template<typename T_>
      static void eval(Tiny::Matrix<T_, 4, 2, 4, 2>& K, const Tiny::Vector<T_, 2, 2>& dx, const T_ s)
      {
        // sigma_1 [11] = dx u_1 = (dx, 0) : (u_1, u_2)
        K(0,0) = s * dx(0);
        K(0,1) = s * T_(0);

        // sigma_2 [12] = 1/2 * (dy u_1 + dx u_2) = 1/2 * (dy, dx) : (u_1, u_2)
        K(1,0) = s * dx(1) / T_(2);
        K(1,1) = s * dx(0) / T_(2);

        // sigma_3 [21] = 1/2 * (dy u_1 + dx u_2) = 1/2 * (dy, dx) : (u_1, u_2)
        K(2,0) = s * dx(1) / T_(2);
        K(2,1) = s * dx(0) / T_(2);

        // sigma_4 [22] = dy u_2 = (0, dy) : (u_1, u_2)
        K(3,0) = s * T_(0),
        K(3,1) = s * dx(1);
      }

      // 2D version with 3 components
      template<typename T_>
      static void eval(Tiny::Matrix<T_, 3, 2, 3, 2>& K, const Tiny::Vector<T_, 2, 2>& dx, const T_ s)
      {
        // sigma_1 [11] = dx u_1 = (dx, 0) : (u_1, u_2)
        K(0,0) = s * dx(0);
        K(0,1) = s * T_(0);

        // sigma_2 [22] = dy u_2 = (0, dy) : (u_1, u_2)
        K(1,0) = s * T_(0),
        K(1,1) = s * dx(1);

        // sigma_3 [12] = 1/2 * (dy u_1 + dx u_2) = 1/2 * (dy, dx) : (u_1, u_2)
        K(2,0) = s * dx(1) / T_(2);
        K(2,1) = s * dx(0) / T_(2);
      }

      OperatorValueType operator()(const TrialBasisData& phi, const TestBasisData& psi)
      {
        // call 3 or 4 component version
        OperatorValueType K;
        eval(K, phi.grad, -psi.value);
        return K;
      }
    }; // class OperatorK
  }; // class OperatorK


  template<int dim, int nsc>
  class OperatorR :
    public Assembly::BilinearOperator
  {
  public:
    static constexpr int BlockHeight = dim;
    static constexpr int BlockWidth = nsc;

    static constexpr TrafoTags trafo_config = TrafoTags::none;
    static constexpr SpaceTags test_config = SpaceTags::value;
    static constexpr SpaceTags trial_config = SpaceTags::grad;

    template<typename AsmTraits_>
    class Evaluator :
      public BilinearOperator::Evaluator<AsmTraits_>
    {
    public:
      typedef typename AsmTraits_::OperatorValueType OperatorValueType;
      typedef typename AsmTraits_::TestBasisData TestBasisData;
      typedef typename AsmTraits_::TrialBasisData TrialBasisData;

    public:
      explicit Evaluator(const OperatorR& DOXY(operat)) {}

      // 2D version with 4 components
      template<typename T_>
      static void eval(Tiny::Matrix<T_, 2, 4, 2, 4>& R, const Tiny::Vector<T_, 2, 2>& dx, const T_ u)
      {
        // u_1 = dx sigma_11 + dy sigma_12 = dx sigma_1 + dy sigma_2
        //     = (dx, dy, 0, 0) : (sigma_11, sigma_12, sigma_21, sigma_22)
        R(0,0) = u * dx(0);
        R(0,1) = u * dx(1);
        R(0,2) = u * T_(0);
        R(0,3) = u * T_(0);

        // u_2 = dx sigma_21 + dy sigma_22 = dx sigma_3 + dy sigma_2
        //     = (0, 0, dx, dy) : (sigma_11, sigma_12, sigma_21, sigma_22)
        R(1,0) = u * T_(0);
        R(1,1) = u * T_(0);
        R(1,2) = u * dx(0);
        R(1,3) = u * dx(1);
      }

      // 2D version with 3 components
      template<typename T_>
      static void eval(Tiny::Matrix<T_, 2, 3, 2, 3>& R, const Tiny::Vector<T_, 2, 2>& dx, const T_ u)
      {
        // u_1 = dx sigma_11 + dy sigma_12 = dx sigma_1 + dy sigma_3
        //     = (dx, 0, dy) : (sigma_11, sigma_22, sigma_12)
        R(0,0) = u * dx(0);
        R(0,1) = u * T_(0);
        R(0,2) = u * dx(1);

        // u_2 = dx sigma_21 + dy sigma_22 = dx sigma_3 + dy sigma_2
        //     = (0, dy, dx) : (sigma_11, sigma_22, sigma_12)
        R(1,0) = u * T_(0);
        R(1,1) = u * dx(1);
        R(1,2) = u * dx(0);
      }

      OperatorValueType operator()(const TrialBasisData& phi, const TestBasisData& psi)
      {
        // call 3 or 4 component version
        OperatorValueType R;
        eval(R, phi.grad, -psi.value);
        return R;
      }
    }; // class OperatorR
  }; // class OperatorR

  // 2D version with 4 stress components
  template<typename Mesh_>
  void add_stress_to_vtk(Geometry::ExportVTK<Mesh_>& exp, const LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, 4>& vs)
  {
    const std::size_t n = vs.size();
    std::vector<double> s11(n), s22(n), s12(n), s21(n);
    const auto* s = vs.elements();
    for(std::size_t i(0); i < n; ++i)
    {
      s11[i] = s[i][0];
      s12[i] = s[i][1];
      s21[i] = s[i][2];
      s22[i] = s[i][3];
    }
    exp.add_vertex_scalar("sigma_11", s11.data());
    exp.add_vertex_scalar("sigma_12", s12.data());
    exp.add_vertex_scalar("sigma_21", s21.data());
    exp.add_vertex_scalar("sigma_22", s22.data());
  }

  // 2D version with 3 stress components
  template<typename Mesh_>
  void add_stress_to_vtk(Geometry::ExportVTK<Mesh_>& exp, const LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, 3>& vs)
  {
    const std::size_t n = vs.size();
    std::vector<double> s11(n), s22(n), s12(n);
    const auto* s = vs.elements();
    for(std::size_t i(0); i < n; ++i)
    {
      s11[i] = s[i][0];
      s22[i] = s[i][1];
      s12[i] = s[i][2];
    }
    exp.add_vertex_scalar("sigma_11", s11.data());
    exp.add_vertex_scalar("sigma_22", s22.data());
    exp.add_vertex_scalar("sigma_12", s12.data());
  }

  // 3D version with 6 stress components
  template<typename Mesh_>
  void add_stress_to_vtk(Geometry::ExportVTK<Mesh_>& exp, const LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, 6>& vs)
  {
    const std::size_t n = vs.size();
    std::vector<double> s11(n), s22(n), s33(n), s12(n), s23(n), s31(n);
    const auto* s = vs.elements();
    for(std::size_t i(0); i < n; ++i)
    {
      s11[i] = s[i][0];
      s22[i] = s[i][1];
      s33[i] = s[i][2];
      s12[i] = s[i][3];
      s23[i] = s[i][4];
      s31[i] = s[i][5];
    }
    exp.add_vertex_scalar("sigma_11", s11.data());
    exp.add_vertex_scalar("sigma_22", s22.data());
    exp.add_vertex_scalar("sigma_33", s33.data());
    exp.add_vertex_scalar("sigma_12", s12.data());
    exp.add_vertex_scalar("sigma_23", s23.data());
    exp.add_vertex_scalar("sigma_31", s31.data());
  }

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  template<int dim_>
  void run_scalar(SimpleArgParser& args);
  template<int dim_>
  void run_2field(SimpleArgParser& args);
  template<int dim_, int nsc_>
  void run_3field(SimpleArgParser& args);

  void main(int argc, char* argv[])
  {
    SimpleArgParser args(argc, argv);

    args.support("level");
    args.support("vtk");
    args.support("scalar");
    args.support("3field");
    args.support("omega");
    args.support("cuda");
    args.support("mu");

    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      std::cerr << std::endl;
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << std::endl;
      Runtime::abort();
    }

    if(args.check("scalar") >= 0)
      run_scalar<2>(args);
    else if(args.check("3field") >= 0)
      run_3field<2, 3>(args);
    else
      run_2field<2>(args);
  }

  template<int dim>
  void run_scalar(SimpleArgParser& args)
  {
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> MeshPartType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
    typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;
    typedef LAFEM::UnitFilter<MemType, DataType, IndexType> FilterType;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    Index level = 4;
    double omega = 0.5;

    args.parse("level", level);
    args.parse("omega", omega);

    std::cout << "Level: " << level << std::endl;
    std::cout << "Omega: " << omega << std::endl;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);

    TrafoType trafo(mesh);
    SpaceType space(trafo);

    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    VectorType vec_sol = matrix.create_vector_r();
    VectorType vec_rhs = matrix.create_vector_l();

    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    Assembly::Common::LaplaceOperator laplace_operator;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, laplace_operator, space, cubature_factory);

    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);
    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, force_functional, space, cubature_factory);

    FilterType filter;

    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(boundary);
    unit_asm.assemble(filter, space);

    filter.filter_mat(matrix);
    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);


    auto new_vanka = Solver::new_amavanka(matrix, filter);
    new_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space));
    auto new_solver = Solver::new_richardson(matrix, filter, omega, new_vanka);

    new_solver->init();
    new_solver->set_plot_mode(Solver::PlotMode::all);
    new_solver->correct(vec_sol, vec_rhs);
    new_solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    String vtk_name(String("./dbg-amavanka-scalar-lvl") + stringify(level));

    if(args.check("vtk") >= 0)
    {
      args.parse("vtk", vtk_name);

      std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

      Geometry::ExportVTK<MeshType> exporter(mesh);
      exporter.add_vertex_scalar("sol", vec_sol.elements());
      exporter.write(vtk_name);
    }
  }

  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */

  template<int dim, int nsc>
  void run_3field(SimpleArgParser& args)
  {
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> MeshPartType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    typedef Space::Lagrange3::Element<TrafoType> SpaceStressType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> VectorVeloType;
    typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorPresType;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, nsc> VectorStressType;

    // tuple vector: (v,s,p)
    typedef LAFEM::TupleVector<VectorVeloType, VectorPresType, VectorStressType> SystemVectorType;

    //    [ A  B  R ]   [ v ]   [ f ]
    //    [ D  .  . ] * [ p ] = [ 0 ]
    //    [ K  .  M ]   [ s ]   [ 0 ]

    // matrices A, B and D for the stokes part
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, dim> MatrixTypeA;
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim,   1> MatrixTypeB;
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType,   1, dim> MatrixTypeD;

    // matrices M, K and R for the stress part
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, nsc, nsc> MatrixTypeM;
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, nsc, dim> MatrixTypeK;
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, nsc> MatrixTypeR;

    // null matrices to fill up the empty parts of the system matrix
    //typedef LAFEM::NullMatrix<MemType, DataType, IndexType, dim, dim> NullMatrixTypeVxV;
    typedef LAFEM::NullMatrix<MemType, DataType, IndexType,   1,   1> NullMatrixTypePxP;
    typedef LAFEM::NullMatrix<MemType, DataType, IndexType, nsc,   1> NullMatrixTypeSxP;
    typedef LAFEM::NullMatrix<MemType, DataType, IndexType,   1, nsc> NullMatrixTypePxS;

    // our system matrix
    typedef LAFEM::TupleMatrix<
      LAFEM::TupleMatrixRow<MatrixTypeA,       MatrixTypeB,       MatrixTypeR>,
      LAFEM::TupleMatrixRow<MatrixTypeD, NullMatrixTypePxP, NullMatrixTypePxS>,
      LAFEM::TupleMatrixRow<MatrixTypeK, NullMatrixTypeSxP,       MatrixTypeM>
    > SystemMatrixType;

    typedef LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, dim> FilterVeloType;
    typedef LAFEM::NoneFilter<MemType, DataType, IndexType> FilterPresType;
    typedef LAFEM::NoneFilterBlocked<MemType, DataType, IndexType, nsc> FilterStressType;

    typedef LAFEM::TupleFilter<FilterVeloType, FilterPresType, FilterStressType> SystemFilterType;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    Index level = 2;
    double omega = 0.5;
    double mu = 0.0;

    args.parse("level", level);
    args.parse("omega", omega);
    args.parse("mu", mu);

    std::cout << "Level: " << level << std::endl;
    std::cout << "Omega: " << omega << std::endl;
    std::cout << "Mu...: " << mu << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    // Create the mesh
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    // Create the boundary mesh-part
    BoundaryHitFunction boundary_hit_function;
    Geometry::HitTestFactory<BoundaryHitFunction, MeshType> boundary_factory(boundary_hit_function, mesh);
    MeshPartType boundary(boundary_factory);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Trafo and Finite Element Space initialisation

    // Create the trafo
    TrafoType trafo(mesh);

    // In the case of Stokes, there are two finite element spaces involved here:
    SpaceStressType space_stress(trafo);
    SpaceVeloType space_velo(trafo); // velocity space
    SpacePresType space_pres(trafo); // pressure space

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    std::cout << "Symbolic Matrix Assembly" << std::endl;

    // As usual, we create first an empty matrix.
    SystemMatrixType matrix;

    MatrixTypeA& matrix_a = matrix.template at<0,0>();
    MatrixTypeB& matrix_b = matrix.template at<0,1>();
    MatrixTypeD& matrix_d = matrix.template at<1,0>();

    MatrixTypeM& matrix_m = matrix.template at<2,2>();
    MatrixTypeR& matrix_r = matrix.template at<0,2>();
    MatrixTypeK& matrix_k = matrix.template at<2,0>();

    // assemble stokes matrix structures
    if(mu > 0.0)
      Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_a, space_velo);
    else
      Assembly::SymbolicAssembler::assemble_matrix_diag(matrix_a, space_velo);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_b, space_velo, space_pres);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_d, space_pres, space_velo);

    // assemble stress matrix structures
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_m, space_stress);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_k, space_stress, space_velo);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_r, space_velo, space_stress);

    // assemble null-matrix structures (necessary to determine the null-matrix sizes)
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix.template at<1,1>(), space_pres);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix.template at<2,1>(), space_stress, space_pres);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix.template at<1,2>(), space_pres, space_stress);

    // create vectors
    SystemVectorType vec_sol = matrix.create_vector_r();
    SystemVectorType vec_rhs = matrix.create_vector_l();

    VectorVeloType& vec_sol_v = vec_sol.template at<0>();
    VectorPresType& vec_sol_p = vec_sol.template at<1>();
    VectorStressType& vec_sol_s = vec_sol.template at<2>();

    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    std::cout << "Numeric Matrix Assembly" << std::endl;

    Cubature::DynamicFactory cubature("auto-degree:7");

    // assemble matrix block A
    //Assembly::Common::LaplaceOperatorBlocked<dim> operator_a; // gradient tensor
    Assembly::Common::DuDvOperatorBlocked<dim>  operator_a; // deformation tensor
    if(mu > 0.0)
      Assembly::BilinearOperatorAssembler::assemble_block_matrix1(matrix_a, operator_a, space_velo, cubature, mu);

    // assemble matrix blocks B and D
    Assembly::GradPresDivVeloAssembler::assemble(matrix_b, matrix_d, space_velo, space_pres, cubature);

    // assemble matrix block M
    //Assembly::Common::IdentityOperatorBlocked<nsc> operator_m;
    //Assembly::BilinearOperatorAssembler::assemble_block_matrix1(matrix_m, operator_m, space_stress, cubature);
    Assembly::OldroydAssembler::assemble_matrix(matrix_m, vec_sol_v, space_velo, space_stress, cubature, 1.0, 0.0, 0.0);

    // assemble matrix block K
    OperatorK<dim,nsc> operator_k;
    Assembly::BilinearOperatorAssembler::assemble_block_matrix2(matrix_k, operator_k, space_stress, space_velo, cubature);

    // assemble matrix block R
    OperatorR<dim,nsc> operator_r;
    Assembly::BilinearOperatorAssembler::assemble_block_matrix2(matrix_r, operator_r, space_velo, space_stress, cubature);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Filter Assembly" << std::endl;

    SystemFilterType filter;

    FilterVeloType& filter_v = filter.template at<0>();
    //FilterPresType& filter_p = filter.template at<1>();
    //FilterStressType& filter_s = filter.template at<2>();

    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // As usual, add the whole boundary to the assembler:
    unit_asm.add_mesh_part(boundary);

    Analytic::Common::ParProfileVector profile_function(0.0, 0.0, 0.0, 1.0, 1.0);
    unit_asm.assemble(filter_v, space_velo, profile_function);

    //Assembly::MeanFilterAssembler::assemble(filter_p, space_pres, cubature);

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    auto vec_sol_2 = vec_sol.clone();

    filter_v.filter_mat(matrix_a);

    filter_v.filter_offdiag_row_mat(matrix_b);
    filter_v.filter_offdiag_row_mat(matrix_r);

    auto new_vanka = Solver::new_amavanka(matrix, filter);
    new_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_velo));
    new_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_pres));
    new_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_stress));
    auto new_solver = Solver::new_richardson(matrix, filter, omega, new_vanka);

    new_solver->init();
    new_solver->set_plot_mode(Solver::PlotMode::all);
    new_solver->correct(vec_sol, vec_rhs);
    new_solver->done();

#ifdef FEAT_HAVE_CUDA
    if(args.check("cuda") >= 0)
    {
      std::cout << std::endl << "Running this stuff on CUDA..." << std::endl;
      typename SystemMatrixType::template ContainerTypeByMDI<Mem::CUDA, DataType, IndexType> matrix_cuda;
      typename SystemFilterType::template FilterTypeByMDI<Mem::CUDA, DataType, IndexType> filter_cuda;
      typename SystemVectorType::template ContainerTypeByMDI<Mem::CUDA, DataType, IndexType> vec_sol_cuda;
      typename SystemVectorType::template ContainerTypeByMDI<Mem::CUDA, DataType, IndexType> vec_rhs_cuda;

      matrix_cuda.convert(matrix);
      filter_cuda.convert(filter);
      vec_sol_cuda.convert(vec_sol_2);
      vec_rhs_cuda.convert(vec_rhs);

      auto cuda_vanka = Solver::new_amavanka(matrix_cuda, filter_cuda);
      cuda_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_velo));
      cuda_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_pres));
      cuda_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_stress));
      auto cuda_solver = Solver::new_richardson(matrix_cuda, filter_cuda, omega, cuda_vanka);

      cuda_solver->init();
      cuda_solver->set_plot_mode(Solver::PlotMode::all);
      cuda_solver->set_plot_name("CUDA Vanka");
      cuda_solver->correct(vec_sol_cuda, vec_rhs_cuda);
      cuda_solver->done();

      std::cout << std::endl << "CUDA Vanka Summary:" << std::endl;
      std::cout << cuda_solver->get_summary() << std::endl;
      std::cout << "Vanka Apply Time: " << stringify_fp_fix(cuda_vanka->time_apply(),3) << std::endl;
    }
#endif // FEAT_HAVE_CUDA

    std::cout << std::endl << "New Vanka Summary:" << std::endl;
    std::cout << new_solver->get_summary() << std::endl;
    std::cout << "Vanka Apply Time: " << stringify_fp_fix(new_vanka->time_apply(),3) << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Analyse velocity field

    std::cout << std::endl << "Velocity Field Analysis" << std::endl;

    Assembly::VelocityInfo<DataType, dim> velo_info = Assembly::VelocityAnalyser::compute(vec_sol_v, space_velo, cubature);

    std::cout << velo_info << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    String vtk_name(String("./dbg-amavanka-3field-lvl") + stringify(level));

    if(args.check("vtk") >= 0)
    {
      args.parse("vtk", vtk_name);

      std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(mesh);

      exporter.add_vertex_vector("velocity", vec_sol_v);

      add_stress_to_vtk(exporter, vec_sol_s);

      VectorPresType cell_pres;
      Assembly::DiscreteCellProjector::project(cell_pres, vec_sol_p, space_pres, cubature);

      // Now we can add the cell-projected pressure to our VTK exporter:
      exporter.add_cell_scalar("pressure", cell_pres.elements());

      // finally, write the VTK file
      exporter.write(vtk_name);
    }

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  } // int main(...)

  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */

  template<int dim>
  void run_2field(SimpleArgParser& args)
  {
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> MeshPartType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> VectorVeloType;
    typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorPresType;

    // tuple vector: (v,s,p)
    typedef LAFEM::TupleVector<VectorVeloType, VectorPresType> SystemVectorType;

    // matrices A, B and D for the stokes part
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, dim> MatrixTypeA;
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim,   1> MatrixTypeB;
    typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType,   1, dim> MatrixTypeD;

    //typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> SystemMatrixType;

    // our system matrix
    typedef LAFEM::TupleMatrix<
      LAFEM::TupleMatrixRow<MatrixTypeA, MatrixTypeB>,
      LAFEM::TupleMatrixRow<MatrixTypeD, LAFEM::NullMatrix<MemType, DataType, IndexType, 1, 1>>
    > SystemMatrixType;


    typedef LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, dim> FilterVeloType;
    typedef LAFEM::NoneFilter<MemType, DataType, IndexType> FilterPresType;

    typedef LAFEM::TupleFilter<FilterVeloType, FilterPresType> SystemFilterType;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    Index level = 2;
    double omega = 0.5;

    args.parse("level", level);
    args.parse("omega", omega);

    std::cout << "Level: " << level << std::endl;
    std::cout << "Omega: " << omega << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    // Create the mesh
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    // Create the boundary mesh-part
    BoundaryHitFunction boundary_hit_function;
    Geometry::HitTestFactory<BoundaryHitFunction, MeshType> boundary_factory(boundary_hit_function, mesh);
    MeshPartType boundary(boundary_factory);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Trafo and Finite Element Space initialisation

    // Create the trafo
    TrafoType trafo(mesh);

    // In the case of Stokes, there are two finite element spaces involved here:
    SpaceVeloType space_velo(trafo); // velocity space
    SpacePresType space_pres(trafo); // pressure space

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    std::cout << "Symbolic Matrix Assembly" << std::endl;

    // As usual, we create first an empty matrix.
    SystemMatrixType matrix;

    MatrixTypeA& matrix_a = matrix.template at<0,0>();
    MatrixTypeB& matrix_b = matrix.template at<0,1>();
    MatrixTypeD& matrix_d = matrix.template at<1,0>();

    // assemble stokes matrix structures
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_a, space_velo);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_b, space_velo, space_pres);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix_d, space_pres, space_velo);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix.template at<1,1>(), space_pres);

    // create vectors
    SystemVectorType vec_sol = matrix.create_vector_r();
    SystemVectorType vec_rhs = matrix.create_vector_l();

    VectorVeloType& vec_sol_v = vec_sol.template at<0>();
    VectorPresType& vec_sol_p = vec_sol.template at<1>();

    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    std::cout << "Numeric Matrix Assembly" << std::endl;

    Cubature::DynamicFactory cubature("auto-degree:7");

    // assemble matrix block A
    Assembly::Common::LaplaceOperatorBlocked<dim> operator_a; // gradient tensor
    //Assembly::Common::DuDvOperatorBlocked<dim>  operator_a; // deformation tensor
    Assembly::BilinearOperatorAssembler::assemble_block_matrix1(matrix_a, operator_a, space_velo, cubature);

    // assemble matrix blocks B and D
    Assembly::GradPresDivVeloAssembler::assemble(matrix_b, matrix_d, space_velo, space_pres, cubature);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Filter Assembly" << std::endl;

    SystemFilterType filter;

    FilterVeloType& filter_v = filter.template at<0>();
    //FilterPresType& filter_p = filter.template at<1>();

    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // As usual, add the whole boundary to the assembler:
    unit_asm.add_mesh_part(boundary);

    Analytic::Common::ParProfileVector profile_function(0.0, 0.0, 0.0, 1.0, 1.0);
    unit_asm.assemble(filter_v, space_velo, profile_function);

    //Assembly::MeanFilterAssembler::assemble(filter_p, space_pres, cubature);

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    filter_v.filter_mat(matrix_a);
    filter_v.filter_offdiag_row_mat(matrix_b);

    auto vec_sol_1 = vec_sol.clone();
    auto vec_sol_2 = vec_sol.clone();


    LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> saddle_matrix;
    saddle_matrix.block_a().clone(matrix.template at<0,0>());
    saddle_matrix.block_b().clone(matrix.template at<0,1>());
    saddle_matrix.block_d().clone(matrix.template at<1,0>());

    //auto new_vanka = std::make_shared<Solver::AmaVanka<SystemMatrixType, SystemFilterType>>(matrix, filter);
    //new_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_velo));
    //new_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_pres));
    auto new_vanka = Solver::new_amavanka(saddle_matrix, filter);
    auto new_solver = Solver::new_richardson(matrix, filter, omega, new_vanka);

    new_solver->init();
    new_solver->set_plot_mode(Solver::PlotMode::iter);
    new_solver->set_plot_name("New Vanka");
    new_solver->correct(vec_sol_1, vec_rhs);
    new_solver->done();

    auto old_vanka = Solver::new_vanka(saddle_matrix, filter, Solver::VankaType::block_full_add);
    auto old_solver = Solver::new_richardson(matrix, filter, omega, old_vanka);

    old_solver->init();
    old_solver->set_plot_mode(Solver::PlotMode::iter);
    old_solver->set_plot_name("Old Vanka");
    old_solver->correct(vec_sol_2, vec_rhs);
    old_solver->done();

#ifdef FEAT_HAVE_CUDA
    if(args.check("cuda") >= 0)
    {
      std::cout << std::endl << "Running this stuff on CUDA..." << std::endl;
      typename SystemMatrixType::template ContainerTypeByMDI<Mem::CUDA, DataType, IndexType> matrix_cuda;
      typename SystemFilterType::template FilterTypeByMDI<Mem::CUDA, DataType, IndexType> filter_cuda;
      typename SystemVectorType::template ContainerTypeByMDI<Mem::CUDA, DataType, IndexType> vec_sol_cuda;
      typename SystemVectorType::template ContainerTypeByMDI<Mem::CUDA, DataType, IndexType> vec_rhs_cuda;

      matrix_cuda.convert(matrix);
      filter_cuda.convert(filter);
      vec_sol_cuda.convert(vec_sol);
      vec_rhs_cuda.convert(vec_rhs);

      auto cuda_vanka = Solver::new_amavanka(matrix_cuda, filter_cuda);
      cuda_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_velo));
      cuda_vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_pres));
      auto cuda_solver = Solver::new_richardson(matrix_cuda, filter_cuda, omega, cuda_vanka);

      cuda_solver->init();
      cuda_solver->set_plot_mode(Solver::PlotMode::iter);
      cuda_solver->set_plot_name("CUDA Vanka");
      cuda_solver->correct(vec_sol_cuda, vec_rhs_cuda);
      cuda_solver->done();

      std::cout << std::endl << "CUDA Vanka Summary:" << std::endl;
      std::cout << cuda_solver->get_summary() << std::endl;
      std::cout << "Vanka Apply Time: " << stringify_fp_fix(cuda_vanka->time_apply(),3) << std::endl;
    }
#endif // FEAT_HAVE_CUDA

    std::cout << std::endl << "New Vanka Summary:" << std::endl;
    std::cout << new_solver->get_summary() << std::endl;
    std::cout << "Vanka Apply Time: " << stringify_fp_fix(new_vanka->time_apply(),3) << std::endl;

    std::cout << std::endl << "Old Vanka Summary:" << std::endl;
    std::cout << old_solver->get_summary() << std::endl;
    std::cout << "Vanka Apply Time: " << stringify_fp_fix(old_vanka->time_apply(),3) << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Analyse velocity field

    vec_sol.copy(vec_sol_1);

    std::cout << std::endl << "Velocity Field Analysis" << std::endl;

    Assembly::VelocityInfo<DataType, dim> velo_info = Assembly::VelocityAnalyser::compute(vec_sol_v, space_velo, cubature);

    std::cout << velo_info << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    String vtk_name(String("./dbg-amavanka-2field-lvl") + stringify(level));

    if(args.check("vtk") >= 0)
    {
      args.parse("vtk", vtk_name);

      std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(mesh);

      exporter.add_vertex_vector("velocity", vec_sol_v);

      VectorPresType cell_pres;
      Assembly::DiscreteCellProjector::project(cell_pres, vec_sol_p, space_pres, cubature);

      // Now we can add the cell-projected pressure to our VTK exporter:
      exporter.add_cell_scalar("pressure", cell_pres.elements());

      // finally, write the VTK file
      exporter.write(vtk_name);
    }

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  } // int main(...)
} // namespace DbgAmaVanka

// Here's our main function
int main(int argc, char* argv[])
{
  Runtime::initialise(argc, argv);
  DbgAmaVanka::main(argc, argv);
  return Runtime::finalise();
}

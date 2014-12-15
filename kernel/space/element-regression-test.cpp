#include <test_system/test_system.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/dirichlet_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/argyris/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/hermite3/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/space/crouzeix_raviart/element.hpp>
#include <kernel/space/bogner_fox_schmit/element.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/proto_solver.hpp>
#include <kernel/util/time_stamp.hpp>

namespace ElementRegression
{
  using namespace FEAST;
  using namespace FEAST::TestSystem;

  typedef double DataType;
  typedef Index IndexType;
  typedef Mem::Main MemType;
  typedef Algo::Generic AlgoType;

  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> FilterType;

  template<bool>
  struct H0Error
  {
    template<typename Vector_, typename Function_, typename Space_, typename Cubature_>
    static double compute(const Vector_&, const Function_&, const Space_&, Cubature_&)
    {
      return 0.0;
    }
  };

  template<bool>
  struct H1Error
  {
    template<typename Vector_, typename Function_, typename Space_, typename Cubature_>
    static double compute(const Vector_&, const Function_&, const Space_&, Cubature_&)
    {
      return 0.0;
    }
  };

  template<bool>
  struct H2Error
  {
    template<typename Vector_, typename Function_, typename Space_, typename Cubature_>
    static double compute(const Vector_&, const Function_&, const Space_&, Cubature_&)
    {
      return 0.0;
    }
  };

  template<>
  struct H0Error<true>
  {
    template<typename Vector_, typename Function_, typename Space_, typename Cubature_>
    static double compute(const Vector_& vector, const Function_& function, const Space_& space, Cubature_& cubature)
    {
      return double(Assembly::ScalarErrorComputerL2::compute(vector, function, space, cubature));
    }
  };

  template<>
  struct H1Error<true>
  {
    template<typename Vector_, typename Function_, typename Space_, typename Cubature_>
    static double compute(const Vector_& vector, const Function_& function, const Space_& space, Cubature_& cubature)
    {
      return double(Assembly::ScalarErrorComputerH1::compute(vector, function, space, cubature));
    }
  };

  template<>
  struct H2Error<true>
  {
    template<typename Vector_, typename Function_, typename Space_, typename Cubature_>
    static double compute(const Vector_& vector, const Function_& function, const Space_& space, Cubature_& cubature)
    {
      return double(Assembly::ScalarErrorComputerH2::compute(vector, function, space, cubature));
    }
  };

  template<typename Shape_>
  struct MeshGen;

  template<int n_>
  struct MeshGen<Shape::Hypercube<n_>>
  {
    static Geometry::ConformalMesh<Shape::Hypercube<n_>> make(Index level)
    {
      Geometry::RefinedUnitCubeFactory<Geometry::ConformalMesh<Shape::Hypercube<n_>>> factory(level);
      return Geometry::ConformalMesh<Shape::Hypercube<n_>>(factory);
    }
  };

  template<int n_>
  struct MeshGen<Shape::Simplex<n_>>
  {
    static Geometry::ConformalMesh<Shape::Simplex<n_>> make(Index level)
    {
      Geometry::RefinedUnitCubeFactory<Geometry::ConformalMesh<Shape::Hypercube<n_>>> factory(level);
      Geometry::ConformalMesh<Shape::Hypercube<n_>> hyper_mesh(factory);
      Geometry::ShapeConvertFactory<Geometry::ConformalMesh<Shape::Simplex<n_>>> convert(hyper_mesh);
      return Geometry::ConformalMesh<Shape::Simplex<n_>>(convert);
    }
  };

  template<int n_>
  struct MeshDisturb;

  template<>
  struct MeshDisturb<1>
  {
    template<typename Mesh_>
    static void apply(Mesh_& mesh, double dsh)
    {
      auto& vtx = mesh.get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        unsigned r0 = (unsigned(i) + 7u) * 3571u; // pseudo-random
        double t = double(1 - int(r0 & 2u)) * dsh;
        if((vtx[i][0] > 0.001) && (vtx[i][0] < 0.999))
          vtx[i][0] += t;
      }
    }
  };

  template<>
  struct MeshDisturb<2>
  {
    template<typename Mesh_>
    static void apply(Mesh_& mesh, double dsh)
    {
      const double fs = 2.0 * Math::pi<double>();
      auto& vtx = mesh.get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        // generate rotation matrix
        unsigned r0 = (unsigned(i) + 7u) * 3571u; // pseudo-random
        double t = fs * double(r0 % 3491u) / 3490.0;
        if((vtx[i][0] > 0.001) && (vtx[i][0] < 0.999))
          vtx[i][0] += dsh * Math::cos(t);
        if((vtx[i][1] > 0.001) && (vtx[i][1] < 0.999))
          vtx[i][1] += dsh * Math::sin(t);
      }
    }
  };

  template<>
  struct MeshDisturb<3>
  {
    template<typename Mesh_>
    static void apply(Mesh_& mesh, double dsh)
    {
      const double fs = 2.0 * Math::pi<double>();
      auto& vtx = mesh.get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        // generate rotation matrix
        unsigned r0 = (unsigned(i) + 7u) * 3571u; // pseudo-random
        unsigned r1 = (unsigned(i) + 5u) * 3433u; // pseudo-random
        double t1 = fs * double(r0 % 3491u) / 3490.0;
        double t2 = -1.0 + double(r1 % 2179u) / 1089.0;
        double t3 = (t2 < 0.0 ? -Math::sqrt(1.0 - t2*t2) : Math::sqrt(1.0 - t2*t2));
        if((vtx[i][0] > 0.001) && (vtx[i][0] < 0.999))
          vtx[i][0] += dsh * Math::cos(t1);
        if((vtx[i][1] > 0.001) && (vtx[i][1] < 0.999))
          vtx[i][1] += dsh * Math::sin(t1);
        if((vtx[i][2] > 0.001) && (vtx[i][2] < 0.999))
          vtx[i][2] += dsh * t3;
      }
    }
  };
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */

  template<typename Shape_, template<typename...> class Element_, bool h0_, bool h1_, bool h2_>
  class ElementRegressionBase :
    public TestSystem::BaseTest
  {
  protected:
    typedef Geometry::ConformalMesh<Shape_> MeshType;
    typedef Geometry::CellSubSet<Shape_> RegionType;
    //typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactory;
    typedef Geometry::BoundaryFactory<MeshType> RegionFactory;

    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Element_<TrafoType> SpaceType;

    typedef Assembly::Common::SineBubbleFunction SolFunction;

    Index _level;
    Cubature::DynamicFactory cubature_factory;
    SolFunction sol_func;
    double h0_ref, h1_ref, h2_ref;

  public:
    explicit ElementRegressionBase(String name, Index level, double h0, double h1, double h2) :
      TestSystem::BaseTest(name),
      _level(level),
      cubature_factory("auto-degree:" + stringify(Math::sqr(SpaceType::local_degree+1)+1)),
      h0_ref(h0),
      h1_ref(h1),
      h2_ref(h2)
    {
    }

    virtual void run() const override
    {
      // run the actual test
      run_error_test();

      // calculate convergence rates; this is not executed during the test run, but it is useful
      // for generating reference errors for the actual test runs.
      //calc_conv_rates();
    }

    virtual void calc_conv_rates() const
    {
      // print out the test name
      std::cout << this->_id << std::endl;

      static constexpr int nl = 6;

      Tiny::Matrix<double, nl, 3> errs;
      TimeStamp stamp;

      //for(Index i(0); i < Index(nl); ++i)
      for(Index i(0); i <= Math::min(_level, Index(nl)); ++i)
      {
        // run test
        errs[i] = run_level(i);

        // plot errors
        std::cout << i << ":";
        if(h0_) std::cout << " " << scientify(errs[i][0], 12);
        if(h1_) std::cout << " " << scientify(errs[i][1], 12);
        if(h2_) std::cout << " " << scientify(errs[i][2], 12);
        std::cout << " |";
        if(i > 0)
        {
          if(h0_) std::cout << " " << scientify(errs[i-1][0] / errs[i][0], 12);
          if(h1_) std::cout << " " << scientify(errs[i-1][1] / errs[i][1], 12);
          if(h2_) std::cout << " " << scientify(errs[i-1][2] / errs[i][2], 12);
        }
        std::cout << std::endl;
      }
      std::cout << "Elapsed Time: " << (TimeStamp().elapsed_string(stamp)) << std::endl;
    }

    virtual void run_error_test() const
    {
      // run test
      auto err = run_level(_level);

      // relative tolerance
      const double eps = 1E-2;

      // check h0-error
      if(h0_)
      {
        //std::cout << "H0-Error: " << scientify(err[0], 12) << std::endl;
        TEST_CHECK_EQUAL_WITHIN_EPS(err[0], h0_ref, eps*h0_ref);
      }
      if(h1_)
      {
        //std::cout << "H1-Error: " << scientify(err[1], 12) << std::endl;
        TEST_CHECK_EQUAL_WITHIN_EPS(err[1], h1_ref, eps*h1_ref);
      }
      if(h2_)
      {
        //std::cout << "H2-Error: " << scientify(err[2], 12) << std::endl;
        TEST_CHECK_EQUAL_WITHIN_EPS(err[2], h2_ref, eps*h2_ref);
      }
    }

    virtual void solve_system(const MatrixType& matrix, const FilterType& filter, VectorType& vec_sol, const VectorType& vec_rhs) const
    {
      // create two empty vector
      VectorType vec_def(vec_sol.clone(LAFEM::CloneMode::Layout));
      VectorType vec_cor(vec_sol.clone(LAFEM::CloneMode::Layout));

      // create a SSOR preconditioner
      LAFEM::SSORPrecond<MatrixType> precon(matrix);

      // create a PCG solver
      LAFEM::PCGSolver<AlgoType, MatrixType, FilterType> solver(matrix, filter, &precon);

      // configure solver
      solver.set_max_iter(1000);
      solver.set_tol_rel(1E-8);
      //solver.set_plot(true);

      // initialise solver
      TEST_CHECK_MSG(solver.init(), "Failed to initialise CG solver!");

      // compute defect
      matrix.template apply<AlgoType>(vec_def, vec_sol, vec_rhs, -DataType(1));
      filter.template filter_def<AlgoType>(vec_def);

      // apply solver on defect
      LAFEM::SolverStatus status = solver.apply(vec_cor, vec_def);
      TEST_CHECK_MSG(LAFEM::status_success(status), "Failed to solve linear system!");

      // update solution
      filter.template filter_cor<AlgoType>(vec_cor);
      vec_sol.template axpy<AlgoType>(vec_cor, vec_sol);

      // release solver
      solver.done();
    }

    Tiny::Vector<double, 3> run_level(Index level) const
    {
      // create mesh
      MeshType mesh(MeshGen<Shape_>::make(level));

      // slighty disturb the mesh
      MeshDisturb<Shape_::dimension>::apply(mesh, 0.05 / double(1 << level));

      // create boundary
      RegionFactory bnd_factory(mesh);
      RegionType boundary(bnd_factory);

      // create trafo
      TrafoType trafo(mesh);

      // create space
      SpaceType space(trafo);

      // create matrix structure
      MatrixType matrix;
      Assembly::SymbolicMatrixAssembler<>::assemble1(matrix, space);

      // create an empty filter
      FilterType filter(space.get_num_dofs());

      // create three vectors
      VectorType vec_sol(matrix.create_vector_r());
      VectorType vec_rhs(matrix.create_vector_r());

      // format everything
      matrix.format();
      vec_sol.format();
      vec_rhs.format();

      // assemble matrix, filter and vectors
      this->assemble_matrix(matrix, space);
      this->assemble_filter(filter, space, boundary);
      this->assemble_rhs(vec_rhs, space);
      this->assemble_sol(vec_sol, space);

      // impose filter onto system
      filter.template filter_mat<AlgoType>(matrix);
      filter.template filter_sol<AlgoType>(vec_sol);
      filter.template filter_rhs<AlgoType>(vec_rhs);

      // solve the system
      this->solve_system(matrix, filter, vec_sol, vec_rhs);

      // error vector
      Tiny::Vector<double, 3> err;

      // compute H0/H1/H2-errors
      err[0] = H0Error<h0_>::compute(vec_sol, this->sol_func, space, this->cubature_factory);
      err[1] = H1Error<h1_>::compute(vec_sol, this->sol_func, space, this->cubature_factory);
      err[2] = H2Error<h2_>::compute(vec_sol, this->sol_func, space, this->cubature_factory);

      // return error vector
      return err;
    }

    virtual void assemble_matrix(MatrixType&, SpaceType&) const = 0;
    virtual void assemble_filter(FilterType&, SpaceType&, RegionType&) const {}
    virtual void assemble_rhs(VectorType&, SpaceType&) const = 0;
    virtual void assemble_sol(VectorType&, SpaceType&) const {}
  };

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */

  template<typename Shape_, template<typename> class Element_, bool h0_, bool h1_, bool h2_>
  class ElementRegressionL2 :
    public ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_>
  {
  public:
    typedef ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_> BaseClass;

    typedef typename BaseClass::RegionType RegionType;
    typedef typename BaseClass::SpaceType SpaceType;

    explicit ElementRegressionL2(String name, Index level, double h0 = 0.0, double h1 = 0.0, double h2 = 0.0) :
      BaseClass(name, level, h0, h1, h2)
    {
    }

    virtual void assemble_matrix(MatrixType& matrix, SpaceType& space) const override
    {
      Assembly::Common::IdentityOperator operat;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, operat, space, this->cubature_factory);
    }

    virtual void assemble_rhs(VectorType& vec_rhs, SpaceType& space) const override
    {
      Assembly::Common::ForceFunctional<typename BaseClass::SolFunction> functional(this->sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, functional, space, this->cubature_factory);
    }
  };

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */

  template<typename Shape_, template<typename> class Element_, bool h0_, bool h1_, bool h2_>
  class ElementRegressionH1 :
    public ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_>
  {
  public:
    typedef ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_> BaseClass;

    typedef typename BaseClass::RegionType RegionType;
    typedef typename BaseClass::SpaceType SpaceType;

    explicit ElementRegressionH1(String name, Index level, double h0 = 0.0, double h1 = 0.0, double h2 = 0.0) :
      BaseClass(name, level, h0, h1, h2)
    {
    }

    virtual void assemble_matrix(MatrixType& matrix, SpaceType& space) const override
    {
      Assembly::Common::LaplaceOperator operat;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, operat, space, this->cubature_factory);
    }

    virtual void assemble_filter(FilterType& filter, SpaceType& space, RegionType& region) const override
    {
      Assembly::DirichletAssembler<SpaceType> dirichlet_asm(space);
      dirichlet_asm.add_cell_set(region);
      dirichlet_asm.assemble(filter);
    }

    virtual void assemble_rhs(VectorType& vec_rhs, SpaceType& space) const override
    {
      const DataType alpha = DataType(Shape_::dimension) * Math::sqr(Math::pi<DataType>());
      Assembly::Common::ForceFunctional<typename BaseClass::SolFunction> functional(this->sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, functional, space, this->cubature_factory, alpha);
    }
  };

  template<typename Trafo_>
  using Discontinuous0 = Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<0>>;

  template<typename Trafo_>
  using RannacherTurekStdNonPar = Space::RannacherTurek::Element<Trafo_>;

  /* ############################################################################################# */
  /* ############################################################################################# */
  /* ############################################################################################# */
  // Discontinuous-0 element (aka P0/Q0)

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Discontinuous0, true, false, false>
    l2_hy1_ldiscontinuous0_lvl5("L2:Discontinuous0:Hypercube<1>:5", 5, 2.019148792005e-002);

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Discontinuous0, true, false, false>
    l2_hy2_ldiscontinuous0_lvl4("L2:Discontinuous0:Hypercube<2>:4", 4, 4.018779996141e-002);

  // L2-Projection Hypercube<3>
  ElementRegressionL2<Shape::Hypercube<3>, Discontinuous0, true, false, false>
    l2_hy3_ldiscontinuous0_lvl3("L2:Discontinuous0:Hypercube<3>:3", 3, 6.885897034156e-002);

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Discontinuous0, true, false, false>
    l2_sx2_ldiscontinuous0_lvl3("L2:Discontinuous0:Simplex<2>:3", 3, 4.647768730126e-002);

  // L2-Projection Simplex<3>
  ElementRegressionL2<Shape::Simplex<3>, Discontinuous0, true, false, false>
    l2_sx3_ldiscontinuous0_lvl2("L2:Discontinuous0:Simplex<3>:2", 2, 6.018518308883e-002);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-1 element (aka P1/Q1)

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Lagrange1::Element, true, true, false>
    l2_hy1_lagrange1_lvl5("L2:Lagrange1:Hypercube<1>:5", 5, 2.665286837321e-004, 6.356129208899e-002);

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::Lagrange1::Element, true, true, false>
    l2_hy2_lagrange1_lvl4("L2:Lagrange1:Hypercube<2>:4", 4, 1.051542601586e-003, 1.271471192054e-001);

  // L2-Projection Hypercube<3>
  ElementRegressionL2<Shape::Hypercube<3>, Space::Lagrange1::Element, true, true, false>
    l2_hy3_lagrange1_lvl3("L2:Lagrange1:Hypercube<3>:3", 3, 3.670654860222e-003, 2.220142781953e-001);

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Lagrange1::Element, true, true, false>
    l2_sx2_lagrange1_lvl3("L2:Lagrange1:Simplex<2>:3", 3, 2.874906340752e-003, 2.375203126593e-001);

  // L2-Projection Simplex<3>
  ElementRegressionL2<Shape::Simplex<3>, Space::Lagrange1::Element, true, true, false>
    l2_sx3_lagrange1_lvl2("L2:Lagrange1:Simplex<3>:2", 2, 1.073303905424e-002, 4.969272469307e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-2 element (aka P2/Q2)

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Lagrange2::Element, true, true, true>
    l2_hy1_lagrange2_lvl4("L2:Lagrange2:Hypercube<1>:4", 4, 2.918342684669e-005, 3.411668271045e-003, 4.000134360015e-001);

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::Lagrange2::Element, true, true, true>
    l2_hy2_lagrange2_lvl3("L2:Lagrange2:Hypercube<2>:3", 3, 2.149226355212e-004, 1.420342358641e-002, 8.128053861467e-001);

  // L2-Projection Hypercube<3>
  ElementRegressionL2<Shape::Hypercube<3>, Space::Lagrange2::Element, true, true, true>
    l2_hy3_lagrange2_lvl2("L2:Lagrange2:Hypercube<3>:2", 2, 1.211653471543e-003, 4.968583665489e-002, 1.402115404765e+000);

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Lagrange2::Element, true, true, true>
    l2_sx2_lagrange2_lvl3("L2:Lagrange2:Simplex<2>:3", 3, 1.477450949611e-004, 1.202305256340e-002, 8.967033891044e-001);

  // L2-Projection Simplex<3>
  ElementRegressionL2<Shape::Simplex<3>, Space::Lagrange2::Element, true, true, true>
    l2_sx3_lagrange2_lvl2("L2:Lagrange2:Simplex<3>:2", 2, 6.920385959247e-004, 5.149687612368e-002, 2.013985124586e+000);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Ranancher-Turek element (aka Q1~)

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, RannacherTurekStdNonPar, true, true, false>
    l2_hy2_rannacher_turek_lvl4("L2:RannacherTurek:Hypercube<2>:4", 4, 1.911701836636e-003, 1.786971993692e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Crouzeix-Raviart element

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::CrouzeixRaviart::Element, true, true, false>
    l2_sx2_crouzeix_raviart_lvl4("L2:CrouzeixRaviart:Simplex<2>:4", 4, 5.736744214046e-004, 1.076734893721e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Bogner-Fox-Schmit element

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::BognerFoxSchmit::Element, true, true, true>
    l2_hy2_bogner_fox_schmit_lvl3("L2:BognerFoxSchmit:Hypercube<2>:3", 3, 9.371691544851e-005, 9.890219682180e-003, 6.424848635514e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Hermite-3 element

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Hermite3::Element, true, true, true>
    l2_hy1_hermite3_lvl5("L2:Hermite3:Hypercube<1>:5", 5, 6.183857974105e-008, 1.250686784444e-005, 2.633222452923e-003);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Argyris element

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Argyris::Element, true, true, true>
    l2_sx2_argyris_lvl3("L2:Argyris:Simplex<2>:3", 3, 9.132110987091e-008, 1.598712842387e-005, 2.195620770059e-003);

  /* ############################################################################################# */
  /* ############################################################################################# */
  /* ############################################################################################# */

  // Lagrange-1 element (aka P1/Q1)

  // H1-Projection Hypercube<1>
  ElementRegressionH1<Shape::Hypercube<1>, Space::Lagrange1::Element, true, true, false>
    h1_hy1_lagrange1_lvl5("H1:Lagrange1:Hypercube<1>:5", 5, 6.374300356367e-004, 6.343140191651e-002);

  // H1-Projection Hypercube<2>
  ElementRegressionH1<Shape::Hypercube<2>, Space::Lagrange1::Element, true, true, false>
    h1_hy2_lagrange1_lvl4("H1:Lagrange1:Hypercube<2>:4", 4, 1.942272640429e-003, 1.268651495575e-001);

  // H1-Projection Hypercube<3>
  ElementRegressionH1<Shape::Hypercube<3>, Space::Lagrange1::Element, true, true, false>
    h1_hy3_lagrange1_lvl3("H1:Lagrange1:Hypercube<3>:3", 3, 5.920955962244e-003, 2.207053505554e-001);

  // H1-Projection Simplex<2>
  ElementRegressionH1<Shape::Simplex<2>, Space::Lagrange1::Element, true, true, false>
    h1_sx2_lagrange1_lvl3("H1:Lagrange1:Simplex<2>:3", 3, 6.259087730394e-003, 2.332900934924e-001);

  // H1-Projection Simplex<3>
  ElementRegressionH1<Shape::Simplex<3>, Space::Lagrange1::Element, true, true, false>
    h1_sx3_lagrange1_lvl1("H1:Lagrange1:Simplex<3>:2", 2, 2.472424266085e-002, 4.803811239931e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-2 element (aka P2/Q2)

  // H1-Projection Hypercube<1>
  ElementRegressionH1<Shape::Hypercube<1>, Space::Lagrange2::Element, true, true, true>
    h1_hy1_lagrange2_lvl4("H1:Lagrange2:Hypercube<1>:4", 4, 3.183068629008e-005, 3.238329097007e-003, 3.967568427665e-001);

  // H1-Projection Hypercube<2>
  ElementRegressionH1<Shape::Hypercube<2>, Space::Lagrange2::Element, true, true, true>
    h1_hy2_lagrange2_lvl3("H1:Lagrange2:Hypercube<2>:3", 3, 2.511253226626e-004, 1.302093337879e-002, 8.012538477278e-001);

  // H1-Projection Hypercube<3>
  ElementRegressionH1<Shape::Hypercube<3>, Space::Lagrange2::Element, true, true, true>
    h1_hy3_lagrange2_lvl1("H1:Lagrange2:Hypercube<3>:2", 2, 1.675103396311e-003, 4.471054331889e-002, 1.379478609799e+000);

  // H1-Projection Simplex<2>
  ElementRegressionH1<Shape::Simplex<2>, Space::Lagrange2::Element, true, true, true>
    h1_sx2_lagrange2_lvl3("H1:Lagrange2:Simplex<2>:3", 3, 1.702956582454e-004, 1.202216232129e-002, 8.981919665319e-001);

  // H1-Projection Simplex<3>
  ElementRegressionH1<Shape::Simplex<3>, Space::Lagrange2::Element, true, true, true>
    h1_sx3_lagrange2_lvl1("H1:Lagrange2:Simplex<3>:2", 2, 1.008388594508e-003, 5.001897374726e-002, 1.934606569160e+000);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Ranancher-Turek element (aka Q1~)

  // H1-Projection Hypercube<2>
  ElementRegressionH1<Shape::Hypercube<2>, RannacherTurekStdNonPar, true, true, false>
    h1_hy2_rannacher_turek_lvl4("H1:RannacherTurek:Hypercube<2>:4", 4, 1.940017656827e-003, 1.787395012244e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Crouzeix-Raviart element

  // H1-Projection Simplex<2>
  ElementRegressionH1<Shape::Simplex<2>, Space::CrouzeixRaviart::Element, true, true, false>
    h1_sx2_crouzeix_raviart_lvl4("H1:CrouzeixRaviart:Simplex<2>:4", 4, 2.028469318821e-003, 1.463272513514e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Hermite-3 element

  // L2-Projection Hypercube<1>
  ElementRegressionH1<Shape::Hypercube<1>, Space::Hermite3::Element, true, true, true>
    h1_hy1_hermite3_lvl5("H1:Hermite3:Hypercube<1>:5", 5, 6.199872813170e-008, 1.243908289727e-005, 2.594740745804e-003);

} // namespace ElementRegression

#include <test_system/test_system.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/util/time_stamp.hpp>

#include <kernel/space/argyris/element.hpp>
#include <kernel/space/bogner_fox_schmit/element.hpp>
#include <kernel/space/crouzeix_raviart/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/hermite3/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>

namespace ElementRegression
{
  using namespace FEAST;
  using namespace FEAST::TestSystem;

  typedef double DataType;
  typedef Index IndexType;
  typedef Mem::Main MemType;

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

  template<typename Shape_, template<typename...> class Element_, bool h0_, bool h1_, bool h2_, typename... ElArgs_>
  class ElementRegressionBase :
    public TestSystem::BaseTest
  {
  protected:
    typedef Geometry::ConformalMesh<Shape_> MeshType;
    typedef Geometry::MeshPart<MeshType> RegionType;
    //typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactory;
    typedef Geometry::BoundaryFactory<MeshType> RegionFactory;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Element_<TrafoType, ElArgs_...> SpaceType;

    typedef Assembly::Common::SineBubbleFunction SolFunction;

    Index _level;
    Cubature::DynamicFactory cubature_factory;
    SolFunction sol_func;
    double h0_ref, h1_ref, h2_ref;

  public:
    explicit ElementRegressionBase(String name, Index level, double h0, double h1, double h2) :
      TestSystem::BaseTest(name + ":" + SpaceType::name() + ":" + Shape_::name() + ":" + stringify(level)),
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

      int nl_min(1);
      int nl_max = (Math::min(int(_level), nl));

      for(int i(nl_min); i <= nl_max; ++i)
      {
        // run test
        errs[i] = run_level(Index(i));

        // plot errors
        std::cout << i << ":";
        if(h0_) std::cout << " " << scientify(errs[i][0], 12);
        if(h1_) std::cout << " " << scientify(errs[i][1], 12);
        if(h2_) std::cout << " " << scientify(errs[i][2], 12);
        std::cout << " |";
        if(i > nl_min)
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

    virtual Tiny::Vector<double, 3> run_level(Index level) const = 0;
  };

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */

  template<typename Shape_, template<typename...> class Element_, bool h0_, bool h1_, bool h2_, typename... ElArgs_>
  class ElementRegressionInterpol :
    public ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_, ElArgs_...>
  {
  public:
    typedef ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_, ElArgs_...> BaseClass;

    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::RegionType RegionType;
    typedef typename BaseClass::TrafoType TrafoType;
    typedef typename BaseClass::SpaceType SpaceType;

    explicit ElementRegressionInterpol(Index level, double h0 = 0.0, double h1 = 0.0, double h2 = 0.0) :
      BaseClass("INT", level, h0, h1, h2)
    {
    }

    virtual Tiny::Vector<double, 3> run_level(Index level) const override
    {
      // create mesh
      MeshType mesh(MeshGen<Shape_>::make(level));

      // slighty disturb the mesh
      MeshDisturb<Shape_::dimension>::apply(mesh, 0.05 / double(1 << level));

      // create boundary
      Geometry::BoundaryFactory<MeshType> bnd_factory(mesh);
      RegionType boundary(bnd_factory);

      // create trafo
      TrafoType trafo(mesh);

      // create space
      SpaceType space(trafo);

      // interpolate our solution function
      VectorType vector(space.get_num_dofs());
      Assembly::Interpolator::project(vector, this->sol_func, space);

      // error vector
      Tiny::Vector<double, 3> err;

      // compute H0/H1/H2-errors
      err[0] = H0Error<h0_>::compute(vector, this->sol_func, space, this->cubature_factory);
      err[1] = H1Error<h1_>::compute(vector, this->sol_func, space, this->cubature_factory);
      err[2] = H2Error<h2_>::compute(vector, this->sol_func, space, this->cubature_factory);

      // return error vector
      return err;
    }
  };

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */

  template<typename Shape_, template<typename...> class Element_, bool h0_, bool h1_, bool h2_, typename... ElArgs_>
  class ElementRegressionSystem :
    public ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_, ElArgs_...>
  {
  public:
    typedef ElementRegressionBase<Shape_, Element_, h0_, h1_, h2_, ElArgs_...> BaseClass;

    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::RegionType RegionType;
    typedef typename BaseClass::TrafoType TrafoType;
    typedef typename BaseClass::SpaceType SpaceType;

    /// relative solver tolerance
    double sol_tol;

    explicit ElementRegressionSystem(String name, Index level, double h0, double h1, double h2, double stol) :
      BaseClass(name, level, h0, h1, h2),
      sol_tol(stol)
    {
    }

    virtual void solve_system(const MatrixType& matrix, const FilterType& filter, VectorType& vec_sol, const VectorType& vec_rhs) const
    {
      // create a SSOR preconditioner
      auto precon = Solver::new_ssor_precond(matrix, filter, DataType(1));

      // create a PCG solver
      Solver::PCG<MatrixType, FilterType> solver(matrix, filter, precon);

      // configure solver
      solver.set_max_iter(1000);
      solver.set_tol_rel(sol_tol);
      //solver.set_plot(true);

      // initialise solver
      solver.init();

      // apply solver
      Solver::Status status = solver.correct(vec_sol, vec_rhs);
      TEST_CHECK_MSG(Solver::status_success(status), "Failed to solve linear system!");

      // release solver
      solver.done();
    }

    virtual Tiny::Vector<double, 3> run_level(Index level) const override
    {
      // create mesh
      MeshType mesh(MeshGen<Shape_>::make(level));

      // slighty disturb the mesh
      MeshDisturb<Shape_::dimension>::apply(mesh, 0.05 / double(1 << level));

      // create boundary
      Geometry::BoundaryFactory<MeshType> bnd_factory(mesh);
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
      filter.filter_mat(matrix);
      filter.filter_sol(vec_sol);
      filter.filter_rhs(vec_rhs);

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

  template<typename Shape_, template<typename...> class Element_, bool h0_, bool h1_, bool h2_, typename... ElArgs_>
  class ElementRegressionL2 :
    public ElementRegressionSystem<Shape_, Element_, h0_, h1_, h2_, ElArgs_...>
  {
  public:
    typedef ElementRegressionSystem<Shape_, Element_, h0_, h1_, h2_, ElArgs_...> BaseClass;

    typedef typename BaseClass::RegionType RegionType;
    typedef typename BaseClass::SpaceType SpaceType;

    explicit ElementRegressionL2(Index level, double h0 = 0.0, double h1 = 0.0, double h2 = 0.0, double stol = 1E-8) :
      BaseClass("L2", level, h0, h1, h2, stol)
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

  template<typename Shape_, template<typename...> class Element_, bool h0_, bool h1_, bool h2_, typename... ElArgs_>
  class ElementRegressionH1 :
    public ElementRegressionSystem<Shape_, Element_, h0_, h1_, h2_, ElArgs_...>
  {
  public:
    typedef ElementRegressionSystem<Shape_, Element_, h0_, h1_, h2_, ElArgs_...> BaseClass;

    typedef typename BaseClass::RegionType RegionType;
    typedef typename BaseClass::SpaceType SpaceType;
    typedef typename BaseClass::MeshType MeshType;

    explicit ElementRegressionH1(Index level, double h0 = 0.0, double h1 = 0.0, double h2 = 0.0, double stol = 1E-8) :
      BaseClass("H1", level, h0, h1, h2, stol)
    {
    }

    virtual void assemble_matrix(MatrixType& matrix, SpaceType& space) const override
    {
      Assembly::Common::LaplaceOperator operat;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, operat, space, this->cubature_factory);
    }

    virtual void assemble_filter(FilterType& filter, SpaceType& space, RegionType& region) const override
    {
      Assembly::UnitFilterAssembler<MeshType> dirichlet_asm;
      dirichlet_asm.add_mesh_part(region);
      dirichlet_asm.assemble(filter, space);
    }

    virtual void assemble_rhs(VectorType& vec_rhs, SpaceType& space) const override
    {
      const DataType alpha = DataType(Shape_::dimension) * Math::sqr(Math::pi<DataType>());
      Assembly::Common::ForceFunctional<typename BaseClass::SolFunction> functional(this->sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, functional, space, this->cubature_factory, alpha);
    }
  };

  /* ############################################################################################# */
  /* ############################################################################################# */
  /* ############################################################################################# */
  // Discontinuous-0 element (aka P0/Q0)

  // Interpolation Hypercube<1>
  ElementRegressionInterpol<Shape::Hypercube<1>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    int_hy1_discontinuous0_lvl5(5, 2.019358459600e-002);

  // Interpolation Hypercube<2>
  ElementRegressionInterpol<Shape::Hypercube<2>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    int_hy2_discontinuous0_lvl4(4, 4.023170882191e-002);

  // Interpolation Hypercube<3>
  ElementRegressionInterpol<Shape::Hypercube<3>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    int_hy3_discontinuous0_lvl3(3, 6.920948509400e-002);

  // Interpolation Simplex<2>
  ElementRegressionInterpol<Shape::Simplex<2>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    int_sx2_discontinuous0_lvl3(3, 4.652869047574e-002);

  // Interpolation Simplex<3>
  ElementRegressionInterpol<Shape::Simplex<3>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    int_sx3_discontinuous0_lvl2(2, 6.042741607084e-002);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Discontinuous-1 element (aka P1dc)

  // Interpolation Hypercube<1>
  ElementRegressionInterpol<Shape::Hypercube<1>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    int_hy1_discontinuous1_lvl5(5, 3.903950715797e-004, 6.343140191651e-002);

  // Interpolation Hypercube<2>
  ElementRegressionInterpol<Shape::Hypercube<2>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    int_hy2_discontinuous1_lvl4(4, 2.506700785425e-003, 1.783826344447e-001);

  // Interpolation Hypercube<3>
  ElementRegressionInterpol<Shape::Hypercube<3>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    int_hy3_discontinuous1_lvl3(3, 1.099318729260e-002, 3.762276070837e-001);

  // Interpolation Simplex<2>
  ElementRegressionInterpol<Shape::Simplex<2>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    int_sx2_discontinuous1_lvl3(3, 7.145236132824e-003, 2.549649415216e-001);

  // Interpolation Simplex<3>
  ElementRegressionInterpol<Shape::Simplex<3>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    int_sx3_discontinuous1_lvl2(2, 2.334296314873e-002, 5.347979318175e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-1 element (aka P1/Q1)

  // Interpolation Hypercube<1>
  ElementRegressionInterpol<Shape::Hypercube<1>, Space::Lagrange1::Element, true, true, false>
    int_hy1_lagrange1_lvl5(5, 6.374300340530e-004, 6.343140191651e-002);

  // Interpolation Hypercube<2>
  ElementRegressionInterpol<Shape::Hypercube<2>, Space::Lagrange1::Element, true, true, false>
    int_hy2_lagrange1_lvl4(4, 3.392065460728e-003, 1.270820247318e-001);

  // Interpolation Hypercube<3>
  ElementRegressionInterpol<Shape::Hypercube<3>, Space::Lagrange1::Element, true, true, false>
    int_hy3_lagrange1_lvl3(3, 1.400674108905e-002, 2.261504131381e-001);

  // Interpolation Simplex<2>
  ElementRegressionInterpol<Shape::Simplex<2>, Space::Lagrange1::Element, true, true, false>
    int_sx2_lagrange1_lvl3(3, 7.145236132824e-003, 2.549649415216e-001);

  // Interpolation Simplex<3>
  ElementRegressionInterpol<Shape::Simplex<3>, Space::Lagrange1::Element, true, true, false>
    int_sx3_lagrange1_lvl2(2, 2.334296314873e-002, 5.347979318175e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-2 element (aka P2/Q2)

  // Interpolation Hypercube<1>
  ElementRegressionInterpol<Shape::Hypercube<1>, Space::Lagrange2::Element, true, true, true>
    int_hy1_lagrange2_lvl4(4, 3.183204796132e-005, 3.238402046955e-003, 3.967721647842e-001);

  // Interpolation Hypercube<2>
  ElementRegressionInterpol<Shape::Hypercube<2>, Space::Lagrange2::Element, true, true, true>
    int_hy2_lagrange2_lvl3(3, 2.545483838641e-004, 1.304807029180e-002, 8.020441638704e-001);

  // Interpolation Hypercube<3>
  ElementRegressionInterpol<Shape::Hypercube<3>, Space::Lagrange2::Element, true, true, true>
    int_hy3_lagrange2_lvl2(2, 1.736674508307e-003, 4.480090797344e-002, 1.381246248538e+000);

  // Interpolation Simplex<2>
  ElementRegressionInterpol<Shape::Simplex<2>, Space::Lagrange2::Element, true, true, true>
    int_sx2_lagrange2_lvl3(3, 1.806543043741e-004, 1.313481388559e-002, 9.791050171640e-001);

  // Interpolation Simplex<3>
  ElementRegressionInterpol<Shape::Simplex<3>, Space::Lagrange2::Element, true, true, true>
    int_sx3_lagrange2_lvl2(2, 1.090590184712e-003, 5.661280319448e-002, 2.129317058613e+000);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Ranancher-Turek element (aka Q1~)

  // Interpolation Hypercube<2>
  ElementRegressionInterpol<Shape::Hypercube<2>, Space::RannacherTurek::Element, true, true, false, Space::RannacherTurek::Variant::StdNonPar>
    int_hy2_rannacher_turek_lvl4(4, 2.500068277725e-003, 1.783313801003e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Crouzeix-Raviart element

  // Interpolation Simplex<2>
  ElementRegressionInterpol<Shape::Simplex<2>, Space::CrouzeixRaviart::Element, true, true, false>
    int_sx2_crouzeix_raviart_lvl4(4, 8.995898938549e-004, 1.100847779732e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Hermite-3 element

  // Interpolation Hypercube<1>
  ElementRegressionInterpol<Shape::Hypercube<1>, Space::Hermite3::Element, true, true, true>
    int_hy1_hermite3_lvl5(5, 1.185419582132e-007, 1.270547785014e-005, 2.568553585246e-003);

  // Interpolation Simplex<2>
  ElementRegressionInterpol<Shape::Simplex<2>, Space::Hermite3::Element, true, true, true>
    int_sx2_hermite3_lvl4(4, 7.698098703027e-007, 1.618142170926e-004, 3.377837471015e-002);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Argyris element

  // Interpolation Simplex<2>
  ElementRegressionInterpol<Shape::Simplex<2>, Space::Argyris::Element, true, true, true>
    int_sx2_argyris_lvl3(3, 2.682673846080e-008, 1.197864886810e-006, 8.035421521830e-005);

  /* ############################################################################################# */
  /* ############################################################################################# */
  /* ############################################################################################# */
  // Discontinuous-0 element (aka P0/Q0)

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    l2_hy1_discontinuous0_lvl5(5, 2.019148792005e-002);

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    l2_hy2_discontinuous0_lvl4(4, 4.018779996141e-002);

  // L2-Projection Hypercube<3>
  ElementRegressionL2<Shape::Hypercube<3>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    l2_hy3_discontinuous0_lvl3(3, 6.885897034156e-002);

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    l2_sx2_discontinuous0_lvl3(3, 4.647768730126e-002);

  // L2-Projection Simplex<3>
  ElementRegressionL2<Shape::Simplex<3>, Space::Discontinuous::Element, true, false, false, Space::Discontinuous::Variant::StdPolyP<0>>
    l2_sx3_discontinuous0_lvl2(2, 6.018518308883e-002);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Discontinuous-1 element (aka P0/Q0)

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    l2_hy1_discontinuous1_lvl5(5, 2.602535197460e-004, 6.343245888649e-002);

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    l2_hy2_discontinuous1_lvl4(4, 1.908779428796e-003, 1.784054134714e-001);

  // L2-Projection Hypercube<3>
  ElementRegressionL2<Shape::Hypercube<3>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    l2_hy3_discontinuous1_lvl3(3, 8.576331121695e-003, 3.754920164364e-001);

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    l2_sx2_discontinuous1_lvl3(3, 2.269755384588e-003, 2.152403144585e-001);

  // L2-Projection Simplex<3>
  ElementRegressionL2<Shape::Simplex<3>, Space::Discontinuous::Element, true, true, false, Space::Discontinuous::Variant::StdPolyP<1>>
    l2_sx3_discontinuous1_lvl2(2, 6.111473872536e-003, 3.568249688054e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-1 element (aka P1/Q1)

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Lagrange1::Element, true, true, false>
    l2_hy1_lagrange1_lvl5(5, 2.665286837321e-004, 6.356129208899e-002);

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::Lagrange1::Element, true, true, false>
    l2_hy2_lagrange1_lvl4(4, 1.051542601586e-003, 1.271471192054e-001);

  // L2-Projection Hypercube<3>
  ElementRegressionL2<Shape::Hypercube<3>, Space::Lagrange1::Element, true, true, false>
    l2_hy3_lagrange1_lvl3(3, 3.670654860222e-003, 2.220142781953e-001);

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Lagrange1::Element, true, true, false>
    l2_sx2_lagrange1_lvl3(3, 2.874906340752e-003, 2.375203126593e-001);

  // L2-Projection Simplex<3>
  ElementRegressionL2<Shape::Simplex<3>, Space::Lagrange1::Element, true, true, false>
    l2_sx3_lagrange1_lvl2(2, 1.073303905424e-002, 4.969272469307e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-2 element (aka P2/Q2)

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Lagrange2::Element, true, true, true>
    l2_hy1_lagrange2_lvl4(4, 2.918342684669e-005, 3.411668271045e-003, 4.000134360015e-001);

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::Lagrange2::Element, true, true, true>
    l2_hy2_lagrange2_lvl3(3, 2.149226355212e-004, 1.420342358641e-002, 8.128053861467e-001);

  // L2-Projection Hypercube<3>
  ElementRegressionL2<Shape::Hypercube<3>, Space::Lagrange2::Element, true, true, true>
    l2_hy3_lagrange2_lvl2(2, 1.211653471543e-003, 4.968583665489e-002, 1.402115404765e+000);

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Lagrange2::Element, true, true, true>
    l2_sx2_lagrange2_lvl3(3, 1.477450949611e-004, 1.202305256340e-002, 8.967033891044e-001);

  // L2-Projection Simplex<3>
  ElementRegressionL2<Shape::Simplex<3>, Space::Lagrange2::Element, true, true, true>
    l2_sx3_lagrange2_lvl2(2, 6.920385959247e-004, 5.149687612368e-002, 2.013985124586e+000);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Ranancher-Turek element (aka Q1~)

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::RannacherTurek::Element, true, true, false, Space::RannacherTurek::Variant::StdNonPar>
    l2_hy2_rannacher_turek_lvl4(4, 1.911701836636e-003, 1.786971993692e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Crouzeix-Raviart element

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::CrouzeixRaviart::Element, true, true, false>
    l2_sx2_crouzeix_raviart_lvl4(4, 5.736744214046e-004, 1.076734893721e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Bogner-Fox-Schmit element

  // L2-Projection Hypercube<2>
  ElementRegressionL2<Shape::Hypercube<2>, Space::BognerFoxSchmit::Element, true, true, true>
    l2_hy2_bogner_fox_schmit_lvl3(3, 9.371691544851e-005, 9.890219682180e-003, 6.424848635514e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Hermite-3 element

  // L2-Projection Hypercube<1>
  ElementRegressionL2<Shape::Hypercube<1>, Space::Hermite3::Element, true, true, true>
    l2_hy1_hermite3_lvl5(5, 6.165546112841e-08, 1.248051929957e-05, 2.626748998399e-03, 1E-10);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Argyris element

  // L2-Projection Simplex<2>
  ElementRegressionL2<Shape::Simplex<2>, Space::Argyris::Element, true, true, true>
    l2_sx2_argyris_lvl3(3, 2.679004360690e-09, 3.692277739529e-07, 5.378755288880e-05, 1E-12);

  /* ############################################################################################# */
  /* ############################################################################################# */
  /* ############################################################################################# */

  // Lagrange-1 element (aka P1/Q1)

  // H1-Projection Hypercube<1>
  ElementRegressionH1<Shape::Hypercube<1>, Space::Lagrange1::Element, true, true, false>
    h1_hy1_lagrange1_lvl5(5, 6.374300356367e-004, 6.343140191651e-002);

  // H1-Projection Hypercube<2>
  ElementRegressionH1<Shape::Hypercube<2>, Space::Lagrange1::Element, true, true, false>
    h1_hy2_lagrange1_lvl4(4, 1.942272640429e-003, 1.268651495575e-001);

  // H1-Projection Hypercube<3>
  ElementRegressionH1<Shape::Hypercube<3>, Space::Lagrange1::Element, true, true, false>
    h1_hy3_lagrange1_lvl3(3, 5.920955962244e-003, 2.207053505554e-001);

  // H1-Projection Simplex<2>
  ElementRegressionH1<Shape::Simplex<2>, Space::Lagrange1::Element, true, true, false>
    h1_sx2_lagrange1_lvl3(3, 6.259087730394e-003, 2.332900934924e-001);

  // H1-Projection Simplex<3>
  ElementRegressionH1<Shape::Simplex<3>, Space::Lagrange1::Element, true, true, false>
    h1_sx3_lagrange1_lvl1(2, 2.472424266085e-002, 4.803811239931e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Lagrange-2 element (aka P2/Q2)

  // H1-Projection Hypercube<1>
  ElementRegressionH1<Shape::Hypercube<1>, Space::Lagrange2::Element, true, true, true>
    h1_hy1_lagrange2_lvl4(4, 3.183068629008e-005, 3.238329097007e-003, 3.967568427665e-001);

  // H1-Projection Hypercube<2>
  ElementRegressionH1<Shape::Hypercube<2>, Space::Lagrange2::Element, true, true, true>
    h1_hy2_lagrange2_lvl3(3, 2.511253226626e-004, 1.302093337879e-002, 8.012538477278e-001);

  // H1-Projection Hypercube<3>
  ElementRegressionH1<Shape::Hypercube<3>, Space::Lagrange2::Element, true, true, true>
    h1_hy3_lagrange2_lvl1(2, 1.675103396311e-003, 4.471054331889e-002, 1.379478609799e+000);

  // H1-Projection Simplex<2>
  ElementRegressionH1<Shape::Simplex<2>, Space::Lagrange2::Element, true, true, true>
    h1_sx2_lagrange2_lvl3(3, 1.702956582454e-004, 1.202216232129e-002, 8.981919665319e-001);

  // H1-Projection Simplex<3>
  ElementRegressionH1<Shape::Simplex<3>, Space::Lagrange2::Element, true, true, true>
    h1_sx3_lagrange2_lvl1(2, 1.008388594508e-003, 5.001897374726e-002, 1.934606569160e+000);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Ranancher-Turek element (aka Q1~)

  // H1-Projection Hypercube<2>
  ElementRegressionH1<Shape::Hypercube<2>, Space::RannacherTurek::Element, true, true, false, Space::RannacherTurek::Variant::StdNonPar>
    h1_hy2_rannacher_turek_lvl4(4, 1.940017656827e-003, 1.787395012244e-001);

  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  /* ********************************************************************************************* */
  // Crouzeix-Raviart element

  // H1-Projection Simplex<2>
  ElementRegressionH1<Shape::Simplex<2>, Space::CrouzeixRaviart::Element, true, true, false>
    h1_sx2_crouzeix_raviart_lvl4(4, 2.028469318821e-003, 1.463272513514e-001);

} // namespace ElementRegression

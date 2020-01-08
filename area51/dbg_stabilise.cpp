// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
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
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/solver/umfpack.hpp>

namespace Andicore
{
  using namespace FEAT;

  //typedef Shape::Hypercube<1> ShapeType;  // 1D
  //typedef Shape::Triangle ShapeType;      // 2D, same as Shape::Simplex<2>
  typedef Shape::Quadrilateral ShapeType;   // 2D, same as Shape::Hypercube<2>
  //typedef Shape::Tetrahedron ShapeType;   // 3D, same as Shape::Simplex<3>
  //typedef Shape::Hexahedron ShapeType;    // 3D, same as Shape::Hypercube<3>

  static constexpr int dim = ShapeType::dimension;

  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;
  //typedef Space::Lagrange2::Element<TrafoType> SpaceType;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
  typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> BlockedVectorType;
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> FilterType;


  template<typename Mesh_>
  void disturb_mesh(Mesh_& mesh, double dsh)
  {
    auto& vtx = mesh.get_vertex_set();

    Random rng;

    for(Index i(0); i < vtx.get_num_vertices(); ++i)
    {
      // generate rotation matrix
      double t(0);
      rng >> t;
      t *= 2.0 * Math::pi<double>();
      if((vtx[i][0] > 0.001) && (vtx[i][0] < 0.999))
        vtx[i][0] += dsh * Math::cos(t);
      if((vtx[i][1] > 0.001) && (vtx[i][1] < 0.999))
        vtx[i][1] += dsh * Math::sin(t);
    }
  }

  class ConvectFunction :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 2;
    typedef Analytic::Image::Vector<2> ImageType;

    static constexpr bool can_value = true;
    static constexpr bool can_grad = false;
    static constexpr bool can_hess = false;

    template<typename EvalTraits_>
    class Evaluator :
      public Analytic::Function::Evaluator<EvalTraits_>
    {
    public:
      typedef typename EvalTraits_::DataType DataType;
      typedef typename EvalTraits_::PointType PointType;
      typedef typename EvalTraits_::ValueType ValueType;
      typedef typename EvalTraits_::HessianType HessianType;

    public:
      const DataType _pi2;

      explicit Evaluator(const ConvectFunction&) :
        _pi2(2.0*Math::pi<DataType>())
      {
      }

      ValueType value(const PointType& p)
      {
        ValueType conv;
        conv[0] =  (1.0 - Math::cos(_pi2*p[0])) * Math::sin(_pi2*p[1]);
        conv[1] = -(1.0 - Math::cos(_pi2*p[1])) * Math::sin(_pi2*p[0]);
        return conv;
      }
    };
  };

  template<typename SolFunction_>
  class AndicoreFunctional :
    public Assembly::LinearFunctional
  {
  public:
    const DataType _nu, _beta, _theta;
    const SolFunction_& _sol_func;
    explicit AndicoreFunctional( const SolFunction_& sol_func, const DataType nu, const DataType beta, const DataType theta) :
      _nu(nu),
      _beta(beta),
      _theta(theta),
      _sol_func(sol_func)
    {
    }

    static constexpr TrafoTags trafo_config = TrafoTags::img_point;
    static constexpr SpaceTags test_config = SpaceTags::value;

    template<typename AsmTraits_>
    class Evaluator :
      public Assembly::LinearFunctional::Evaluator<AsmTraits_>
    {
    public:
      typedef typename AsmTraits_::TrafoData TrafoData;
      typedef typename AsmTraits_::TestBasisData TestBasisData;
      typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;
      typedef typename AsmTraits_::DataType DataType;
      typedef Analytic::EvalTraits<DataType, SolFunction_> AnalyticEvalTraits;
      typedef typename AnalyticEvalTraits::ValueType ValueType;
      typedef typename SolFunction_::template Evaluator<AnalyticEvalTraits> SolEvaluator;

    protected:
      const DataType _pi2, _nu, _beta, _theta;
      SolEvaluator _sol_eval;
      ValueType _force_value;

    public:
      explicit Evaluator(const AndicoreFunctional& functional) :
        _pi2(2.0*Math::pi<DataType>()),
        _nu(functional._nu),
        _beta(functional._beta),
        _theta(functional._theta),
        _sol_eval(functional._sol_func),
        _force_value(DataType(0))
      {
      }

      void set_point(const TrafoData& tau)
      {
        typename AnalyticEvalTraits::ValueType    value = _sol_eval.value(tau.img_point);
        typename AnalyticEvalTraits::GradientType grad  = _sol_eval.gradient(tau.img_point);
        typename AnalyticEvalTraits::HessianType  hess  = _sol_eval.hessian(tau.img_point);

        // evaluate convection
        typename AnalyticEvalTraits::GradientType conv;
        conv[0] =  (1.0 - Math::cos(_pi2*tau.img_point[0])) * Math::sin(_pi2*tau.img_point[1]);
        conv[1] = -(1.0 - Math::cos(_pi2*tau.img_point[1])) * Math::sin(_pi2*tau.img_point[0]);

        _force_value = -_nu * hess.trace() + _beta * dot(conv, grad) + _theta*value;
      }

      ValueType operator()(const TestBasisData& psi) const
      {
        return _force_value * psi.value;
      }
    }; // class AndicoreFunctional<...>::Evaluator<...>
  }; // class AndicoreFunctional<...>


  void main(int argc, char** argv)
  {
    SimpleArgParser args(argc, argv);

    args.support("level");
    args.support("nu");
    args.support("beta");
    args.support("theta");
    args.support("upsam");
    args.support("gamma");

    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      std::cerr << std::endl;
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << std::endl;
      Runtime::abort();
    }

    DataType nu = 1E-5;   // diffusion
    DataType beta = 1.0;  // convection
    DataType theta = 0.0; // reaction
    DataType upsam = 0.1; // stabilisation
    DataType gamma = 0.01; // EOJ stabilisation
    Index level = 5;

    args.parse("level", level);
    args.parse("nu", nu);
    args.parse("beta", beta);
    args.parse("theta", theta);
    args.parse("upsam", upsam);
    args.parse("gamma", gamma);

    std::cout << String("Level").pad_back(20, '.') << ": " << level << std::endl;
    std::cout << String("nu").pad_back(20, '.') << ": " << nu << std::endl;
    std::cout << String("beta").pad_back(20, '.') << ": " << beta << std::endl;
    std::cout << String("theta").pad_back(20, '.') << ": " << theta << std::endl;
    std::cout << String("upsam").pad_back(20, '.') << ": " << upsam << std::endl;
    std::cout << String("gamma").pad_back(20, '.') << ": " << gamma << std::endl;

    // ********************************************************************************************

    Cubature::DynamicFactory cubature(String("auto-degree:") + stringify(2*(SpaceType::local_degree)+3));

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    // disturb mesh
    disturb_mesh(mesh, 0.15 / double(1 << level));

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);

    TrafoType trafo(mesh);
    SpaceType space(trafo);



    ConvectFunction conv_function;
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;
    AndicoreFunctional<decltype(sol_function)> rhs_functional(sol_function, nu, beta, theta);


    MatrixType matrix;

    //Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
    Assembly::SymbolicAssembler::assemble_matrix_ext1(matrix, space);

    VectorType vec_sol_1 = matrix.create_vector_r();
    VectorType vec_sol_2 = matrix.create_vector_r();
    VectorType vec_sol_3 = matrix.create_vector_r();
    VectorType vec_rhs = matrix.create_vector_l();

    BlockedVectorType vec_conv(matrix.rows());
    Assembly::Interpolator::project(vec_conv, conv_function, space);

    vec_sol_1.format();
    vec_sol_2.format();
    vec_sol_3.format();
    vec_rhs.format();

    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, rhs_functional, space, cubature);

    FilterType filter;
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(boundary);
    unit_asm.assemble(filter, space);

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol_1);
    filter.filter_sol(vec_sol_2);


    Assembly::BurgersAssembler<DataType, IndexType, 2> burgers;
    Assembly::TraceAssembler<TrafoType> trace_asm(trafo);
    trace_asm.compile_all_facets(true, false);


#ifdef FEAT_HAVE_UMFPACK

    auto solver = Solver::new_umfpack(matrix);
    solver->init_symbolic();


    // unstabilised
    std::cout << "Solving Unstabilised System..." << std::endl;
    burgers.nu = nu;
    burgers.beta = beta;
    burgers.theta = theta;
    matrix.format();
    burgers.assemble_scalar_matrix(matrix, vec_conv, space, cubature);
    filter.filter_mat(matrix);

    solver->init_numeric();
    Solver::solve(*solver, vec_sol_1, vec_rhs, matrix, filter);
    solver->done_numeric();

    Assembly::ScalarErrorInfo<DataType> errors_1 = Assembly::ScalarErrorComputer<1>::
      compute(vec_sol_1, sol_function, space, cubature);
    std::cout << errors_1 << std::endl;


    // streamline diffusion
    std::cout << "Solving Streamline Diffusion Stabilised System..." << std::endl;
    burgers.sd_delta = upsam;
    burgers.sd_nu = nu;
    burgers.set_sd_v_norm(vec_conv);
    matrix.format();
    burgers.assemble_scalar_matrix(matrix, vec_conv, space, cubature);
    filter.filter_mat(matrix);

    solver->init_numeric();
    Solver::solve(*solver, vec_sol_2, vec_rhs, matrix, filter);
    solver->done_numeric();

    Assembly::ScalarErrorInfo<DataType> errors_2 = Assembly::ScalarErrorComputer<1>::
      compute(vec_sol_2, sol_function, space, cubature);
    std::cout << errors_2 << std::endl;


    // EOJ stabilisation
    std::cout << "Solving Jump Stabilised System..." << std::endl;
    burgers.nu = nu;
    burgers.beta = beta;
    burgers.theta = theta;
    burgers.sd_nu = 0.0;
    burgers.sd_delta = 0.0;
    matrix.format();
    burgers.assemble_scalar_matrix(matrix, vec_conv, space, cubature);
    trace_asm.assemble_jump_stabil_operator_matrix(matrix, space, cubature, gamma, DataType(2), DataType(2));
    filter.filter_mat(matrix);

    solver->init_numeric();
    Solver::solve(*solver, vec_sol_3, vec_rhs, matrix, filter);
    solver->done_numeric();

    Assembly::ScalarErrorInfo<DataType> errors_3 = Assembly::ScalarErrorComputer<1>::
      compute(vec_sol_3, sol_function, space, cubature);
    std::cout << errors_3 << std::endl;

    solver->done_symbolic();

#else // no UMFPACK
    std::cout << "You need to compile with UMFPACK support." << std::endl;
#endif // FEAT_HAVE_UMFPACK

    if(args.check("vtk") >= 0)
    {
      String vtk_name(String("./dbg-stabilise-lvl") + stringify(level));
      args.parse("vtk", vtk_name);

      std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

      Geometry::ExportVTK<MeshType> exporter(mesh);
      exporter.add_vertex_vector("conv", vec_conv);
      exporter.add_vertex_scalar("sol_1", vec_sol_1.elements());
      exporter.add_vertex_scalar("sol_2", vec_sol_2.elements());
      exporter.add_vertex_scalar("sol_3", vec_sol_3.elements());
      exporter.add_vertex_scalar("rhs", vec_rhs.elements());
      exporter.write(vtk_name);
    }

    std::cout << "Finished!" << std::endl;
  }
} // namespace Andicore

// Here's our main function
int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  Andicore::main(argc, argv);
  return FEAT::Runtime::finalise();
}

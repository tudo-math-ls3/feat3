// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/string.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/pcg.hpp>

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial01
{
  typedef Shape::Quadrilateral ShapeType;
  typedef double DataType;
  typedef Mem::Main MemType;

  static double PI = Math::pi<double>();

  template<typename DT_>
  struct NeumannSineRhs
  {
    static DT_ eval(DT_ x, DT_ y)
    {
      return -PI*(Math::sin(PI * x) + Math::sin(PI * y));
    }
  };

  typedef Analytic::StaticWrapperFunction<2, NeumannSineRhs, true, false, false> NeumannRhsFunction;

  void main(Index level)
  {
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> BoundaryType;
    typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactoryType;
    typedef Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;
    MeshFactoryType mesh_factory(level);
    MeshType mesh(mesh_factory);

    std::cout << "Creating Boundary..." << std::endl;
    BoundaryFactoryType boundary_factory(mesh);
    BoundaryType boundary(boundary_factory);

    std::cout << "Creating Trafo..." << std::endl;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    TrafoType trafo(mesh);

    std::cout << "Creating Space..." << std::endl;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;
    SpaceType space(trafo);

    typedef LAFEM::DenseVector<MemType, DataType> VectorType;
    typedef LAFEM::SparseMatrixCSR<MemType, DataType> MatrixType;
    typedef LAFEM::NoneFilter<MemType, DataType, Index> FilterType;
    //typedef LAFEM::MeanFilter<MemType, DataType, Index> FilterType;

    FilterType filter;

    std::cout << "Allocating and initialising vectors and matrix..." << std::endl;
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
    VectorType vec_sol(space.get_num_dofs());
    VectorType vec_rhs(space.get_num_dofs());

    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    std::cout << "Assembling system matrix..." << std::endl;
    matrix.format();

    Assembly::Common::LaplaceOperator laplace_operator;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, laplace_operator, space, cubature_factory);

    std::cout << "Assembling right-hand-side vector..." << std::endl;
    vec_rhs.format();
    {
      Analytic::Common::SineBubbleFunction<2> sol_function;
      Assembly::Common::ForceFunctional<decltype(sol_function)> force_functional(sol_function);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, force_functional, space, cubature_factory,
        2.0 * Math::sqr(Math::pi<DataType>()));
    }

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    {
      Assembly::TraceAssembler<TrafoType> trace_asm(trafo);
      trace_asm.add_mesh_part(boundary);
      trace_asm.compile_facets();
      NeumannRhsFunction neu_func;
      Assembly::Common::ForceFunctional<decltype(neu_func)> neu_force(neu_func);
      trace_asm.assemble_functional_vector(vec_rhs, neu_force, space, cubature_factory);
    }

    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    vec_sol.format(4.0 / Math::sqr(PI));
    std::cout << "Assembling boundary conditions..." << std::endl;

    std::cout << "Solving linear system..." << std::endl;
    //auto precond = Solver::new_ssor_precond(matrix, filter);
    //auto precond = Solver::new_ilu_precond(matrix, filter);
    auto solver = Solver::new_pcg(matrix, filter/*, precond*/);

    solver->set_plot_mode(Solver::PlotMode::iter);
    solver->init();
    solver->set_max_iter(1000);

    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    solver->done();

    {
      Analytic::Common::SineBubbleFunction<2> sol_function;

      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::
        compute(vec_sol,sol_function, space, cubature_factory);

      std::cout << errors << std::endl;
    }

    // First of all, build the filename string
    String vtk_name(String("./dbg-trace-2-lvl") + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // add the vertex-projection of our solution and rhs vectors
    exporter.add_vertex_scalar("sol", vec_sol.elements());
    exporter.add_vertex_scalar("rhs", vec_rhs.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  } // int main(...)
} // namespace Tutorial01

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialize the FEAT runtime environment:
  Runtime::initialize(argc, argv);

  // Specify the desired mesh refinement level, defaulted to 3.
  // Note that FEAT uses its own "Index" type rather than a wild mixture of int, uint, long
  // and such.
  Index level(3);

  // Now let's see if we have command line parameters: This tutorial supports passing
  // the refinement level as a command line parameter, to investigate the behaviour of the L2/H1
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
  Tutorial01::main(level);

  // And finally, finalize our runtime environment. This function returns the 'EXIT_SUCCESS' return code,
  // so we can simply return this as the result of our main function to indicate a successful run.
  return Runtime::finalize();
}

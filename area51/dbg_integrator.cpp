#include <kernel/util/string.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/stop_watch.hpp>
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
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/domain_assembler_helpers.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/pcg.hpp>

using namespace FEAT;

namespace Tutorial01
{
  typedef Shape::Quadrilateral ShapeType;   // 2D, same as Shape::Hypercube<2>
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> FilterType;


  // Here's our tutorial's main function.
  void main(int argc, char** argv)
  {
    Index level = 4;
    std::size_t cores = 2;

    if(argc > 1) String(argv[1]).parse(level);
    if(argc > 2) String(argv[2]).parse(cores);

    std::cout << "Level: " << level << std::endl;
    std::cout << "Cores: " << cores << std::endl;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);

    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);

    MeshPartType boundary(boundary_factory);

    // Let's create a trafo object now. Its only parameter is the mesh that it is defined on.
    TrafoType trafo(mesh);

    // Create the desire finite element space. Its only parameter is the trafo that it is defined on.
    SpaceType space(trafo);


    // create matrix
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
    matrix.format();

    String cub_name = "auto-degree:3";

    Cubature::DynamicFactory cubature_factory(cub_name);
    Assembly::Common::LaplaceOperator laplace_op;

    StopWatch watch_1, watch_2, watch_3, watch_4, watch_5;

    watch_1.start();
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, laplace_op, space, cubature_factory, -4.0);
    watch_1.stop();

    std::cout << "Time 1: " << watch_1.elapsed_string().pad_front(10) << std::endl;

    // create integral assembler
    if(true)
    {
      Assembly::DomainAssembler<TrafoType> integrator(trafo);

      integrator.set_max_worker_threads(cores);
      integrator.set_threading_strategy(Assembly::ThreadingStrategy::layered);
      integrator.compile_all_elements();

      Assembly::BilinearOperatorMatrixAssemblyJob1<decltype(laplace_op), MatrixType, SpaceType>
        laplace_job(laplace_op, matrix, space, cub_name, 1.0);

      watch_2.start();
      integrator.assemble_master(laplace_job);
      watch_2.stop();
      std::cout << "Time 2: " << watch_2.elapsed_string().pad_front(10) << std::endl;

      watch_3.start();
      //integrator.assemble(laplace_job);
      assemble_bilinear_operator_matrix_1(integrator, matrix, laplace_op, space, cub_name);
      watch_3.stop();
      std::cout << "Time 3: " << watch_3.elapsed_string().pad_front(10) << std::endl;
    }


    // create integral assembler
    if(true)
    {
      Assembly::DomainAssembler<TrafoType> integrator(trafo);

      integrator.set_max_worker_threads(cores);
      integrator.set_threading_strategy(Assembly::ThreadingStrategy::colored);
      integrator.compile_all_elements();

      Assembly::BilinearOperatorMatrixAssemblyJob1<decltype(laplace_op), MatrixType, SpaceType>
        laplace_job(laplace_op, matrix, space, cub_name, 1.0);

      watch_4.start();
      integrator.assemble_master(laplace_job);
      watch_4.stop();
      std::cout << "Time 4: " << watch_4.elapsed_string().pad_front(10) << std::endl;

      watch_5.start();
      //integrator.assemble(laplace_job);
      assemble_bilinear_operator_matrix_1(integrator, matrix, laplace_op, space, cub_name);
      watch_5.stop();
      std::cout << "Time 5: " << watch_5.elapsed_string().pad_front(10) << std::endl;
    }

    if(true)
    {
      LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 2, 2> matrix_b;
      Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_b, space);
      matrix_b.format();

      Assembly::Common::LaplaceOperatorBlocked<2> laplace_op_b;

      Assembly::DomainAssembler<TrafoType> integrator(trafo);

      integrator.compile_all_elements();


      /*Assembly::BilinearOperatorMatrixAssemblyJob1<decltype(laplace_op_b), decltype(matrix_b), SpaceType>
        laplace_job(laplace_op_b, matrix_b, space, cub_name, 1.0);*/

      watch_2.start();
      //integrator.assemble_master(laplace_job);
      assemble_bilinear_operator_matrix_1(integrator, matrix_b, laplace_op_b, space, cub_name, 1.0);
      assemble_bilinear_operator_matrix_2(integrator, matrix_b, laplace_op_b, space, space, cub_name, -1.0);
      watch_2.stop();
      std::cout << "Time 2: " << watch_2.elapsed_string().pad_front(10) << std::endl;

      DataType vmax = 0.0;
      auto* val = matrix_b.val();
      for(Index i(0); i < matrix.used_elements(); ++i)
        vmax = Math::max(vmax, Math::abs(val[i].norm_frobenius_sqr()));
      vmax = Math::sqrt(vmax);

      std::cout << "MAX-ERR = " << stringify_fp_sci(vmax) << std::endl;

    }

    DataType vmax = 0.0;
    DataType* val = matrix.val();
    for(Index i(0); i < matrix.used_elements(); ++i)
      vmax = Math::max(vmax, Math::abs(val[i]));

    std::cout << "MAX-ERR = " << stringify_fp_sci(vmax) << std::endl;
  } // void main(...)
} // namespace Tutorial01

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialize the FEAT runtime environment:
  Runtime::initialize(argc, argv);

  // call the tutorial's main function
  Tutorial01::main(argc, argv);

  // And finally, finalize our runtime environment. This function returns the 'EXIT_SUCCESS' return code,
  // so we can simply return this as the result of our main function to indicate a successful run.
  return Runtime::finalize();
}

#include <kernel/base_header.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/string.hpp>                          // for String
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>            // for RefinedUnitCubeFactory
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-2 Element (aka "Q2")
#include <kernel/space/lagrange3/element.hpp>              // the Lagrange-3 Element (aka "Q3")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp>    // for Assembly::assemble_***
#include <kernel/assembly/common_operators.hpp>

#include <kernel/util/simple_arg_parser.hpp>               // NEW: for SimpleArgParser
#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;



void run_csr_io(int level)
{
  typedef Shape::Quadrilateral ShapeType;   // 2D, same as Shape::Hypercube<2>
  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  // typedef Geometry::MeshPart<MeshType> MeshPartType;

  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // typedef Space::Lagrange1::Element<TrafoType> SpaceType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceType;
  typedef double DataType;
  typedef Index IndexType;
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;
  Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);

  MeshType mesh(mesh_factory);
  TrafoType trafo(mesh);
  SpaceType space(trafo);
  MatrixType matrix;
  Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
  matrix.format();
  Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
  domain_assembler.compile_all_elements();
  String cubature_name =
    "auto-degree:5";          // automatic cubature rule for 5th degree polynomials
  Assembly::Common::LaplaceOperator laplace_operator;
  Assembly::assemble_bilinear_operator_matrix_1(
    domain_assembler, matrix, laplace_operator, space, cubature_name);

  MatrixType mat_1;
  MatrixType mat_2;


  const String file_name = "bench_example_file";
  auto write_out_test = [&file_name, &matrix](){matrix.write_out(FileMode::fm_mtx, file_name + ".mtx");};
  auto write_out_test_bin = [&file_name, &matrix](){matrix.write_out(FileMode::fm_binary, file_name + ".bin");};

  auto read_from_test = [&file_name, &mat_1]() {mat_1.read_from(FileMode::fm_mtx, file_name + ".mtx");};
  auto read_from_test_bin = [&file_name, &mat_2]() {mat_2.read_from(FileMode::fm_binary, file_name + ".bin");};

  std::cout << "CSR Write out MTX\n";
  run_bench(write_out_test, 1., 1.);
  std::cout << "CSR Write out binary\n";
  run_bench(write_out_test_bin, 1., 1.);
  std::cout << "CSR Read from MTX\n";
  run_bench(read_from_test, 1., 1.);
  std::cout << "CSR Read from binary\n";
  run_bench(read_from_test_bin, 1., 1.);
}





int main(int argc, char** argv)
{

  Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  SimpleArgParser parser(argc, argv);
  parser.support("lvl");
  int lvl = 1;
  if(parser.parse("lvl", lvl) < 0)
  {
    std::cout << "Could not parse lvl\n";
    Runtime::abort();
  }
  if(lvl < 1)
  {
    std::cout << "Required to provide positive integer as lvl\n";
    Runtime::abort();
  }
  run_csr_io(lvl);
  return 0;
}
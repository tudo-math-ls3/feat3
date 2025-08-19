// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#ifdef FEAT_HAVE_OMP
#include <omp.h>
#else
inline static int omp_get_max_threads()
{
  return 1;
}
#endif

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

#ifdef BENCH_FP_32
typedef float DataType;
const char * fp_string = "fp_32";
#else
typedef double DataType;
const char * fp_string = "fp_64";
#endif

#ifdef BECNH_IT_32
typedef std::uint32_t IndexType;
const char * ix_string = "int_32";
#else
typedef std::uint64_t IndexType;
const char * ix_string = "int_64";
#endif

template <typename MatrixType_>
void run(const MatrixType_& matrix, PreferredBackend backend)
{
  Backend::set_preferred_backend(PreferredBackend::generic);
  auto vec_rhs = matrix.create_vector_r();
  auto vec_sol = matrix.create_vector_l();

  vec_rhs.format();
  vec_sol.format();

  Index used_elements(matrix.template used_elements<Perspective::pod>());

  for (Index i (0) ; i < vec_rhs.size() ; ++i)
  {
    decltype(vec_rhs(i)) temp{DataType(i%100) / DataType(100)};
    vec_rhs(i, temp);
  }


  double flops = double(used_elements);
  flops *= 2;

  double bytes = double(matrix.bytes()) + 2 * double(vec_rhs.bytes());
  std::printf("Flops per run: %.0f MFlops, Bytes per run: %.0fMB\n", flops/1E+6, bytes/(1024*1024));

  Backend::set_preferred_backend(backend);
  auto func = [&] () { matrix.apply(vec_sol, vec_rhs); };
  FEAT_APPLICATION_MARKER_START("matrix-apply");
  run_bench(func, flops, bytes);
  FEAT_APPLICATION_MARKER_STOP("matrix-apply");

  std::cout<<"control norm: "<< vec_sol.norm2()<<"\n";
}

int main(int argc, char ** argv)
{
  FEAT_APPLICATION_MARKER_REGISTER("matrix-apply");
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  SimpleArgParser parser(argc, argv);
  parser.support("mtx-file");
  parser.support("backend");
  parser.support("block-height");
  parser.support("block-width");

  if(parser.check("mtx-file") <= 0)
  {
    XABORTM("mtx-file argument required");
  }
  String mtx_file = parser.query("mtx-file")->second.front();

  PreferredBackend backend = PreferredBackend::generic;
  String backend_string = "generic";
  if(parser.check("backend") > 0)
  {
    backend_string = parser.query("backend")->second.front();
    if(backend_string.compare_no_case("generic") == 0)
    {
      backend = PreferredBackend::generic;
    }
    else if(backend_string.compare_no_case("cuda") == 0)
    {
      backend = PreferredBackend::cuda;
    }
    else if(backend_string.compare_no_case("mkl") == 0)
    {
      backend = PreferredBackend::mkl;
    }
    else
    {
      XABORTM("No known backend option " + backend_string);
    }
  }

  int height{0}, width{0};
  height = parser.parse_default("block-height", int(1));
  width = parser.parse_default("block-width", int(1));


  // read in mtx file
  LAFEM::SparseMatrixCSR<DataType, IndexType> in_matrix(LAFEM::FileMode::fm_mtx, mtx_file);

  std::printf("Num threads %i, Dataype %s, IndexType %s, backend %s, real rows %i, real cols %i \n", omp_get_max_threads(), fp_string, ix_string, backend_string.data(), int(in_matrix.rows()), int(in_matrix.columns()));
  if(height == 1 && width == 1)
    std::printf("Plain CRS: Block Height: %i,  Block Widht: %i \n", height, width);
  else
    std::printf("BCRS: Block Height: %i,  Block Widht: %i \n", height, width);

  if(height == 1 && width == 1)
  {
    run(in_matrix, backend);
  }
  else if(height == 2 && width == 1)
  {
    LAFEM::SparseMatrixBCSR<DataType, IndexType, 2, 1> matrix;
    matrix.convert(in_matrix);
    run(in_matrix, backend);
  }
  else if(height == 1 && width == 2)
  {
    LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, 2> matrix;
    matrix.convert(in_matrix);
    run(in_matrix, backend);
  }
  else if(height == 2 && width == 2)
  {
    LAFEM::SparseMatrixBCSR<DataType, IndexType, 2, 2> matrix;
    matrix.convert(in_matrix);
    run(in_matrix, backend);
  }
  else if(height == 3 && width == 1)
  {
    LAFEM::SparseMatrixBCSR<DataType, IndexType, 3, 1> matrix;
    matrix.convert(in_matrix);
    run(in_matrix, backend);
  }
  else if(height == 1 && width == 3)
  {
    LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, 3> matrix;
    matrix.convert(in_matrix);
    run(in_matrix, backend);
  }
  else if(height == 3 && width == 3)
  {
    LAFEM::SparseMatrixBCSR<DataType, IndexType, 3, 3> matrix;
    matrix.convert(in_matrix);
    run(in_matrix, backend);
  }
  else
  {
    std::printf("Combination height %i, width %i not implemented.\n", int(height), int(width));
    Runtime::abort();
  }

  return 0;
}
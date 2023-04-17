// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/solver/pcg.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/runtime.hpp>

#include <iostream>

using namespace FEAT;


template <typename DT_, typename IT_>
void run(PreferredBackend backend)
{
  Runtime::set_preferred_backend(PreferredBackend::generic);

  std::cout << String(100, '=') << std::endl;
  std::cout << "Backend..: " << backend << std::endl;
  std::cout << "DataType.: " << Type::Traits<DT_>::name() << std::endl;
  std::cout << "IndexType: " << Type::Traits<IT_>::name() << std::endl;

  // create matrix
  LAFEM::PointstarFactoryFE<DT_, IT_> factory(Index(2000));
  LAFEM::SparseMatrixCSR<DT_, IT_> matrix = factory.matrix_csr();

  std::cout << "NumDofs..: " << matrix.rows() << std::endl;
  std::cout << "NumNZEs..: " << matrix.used_elements() << std::endl;

  // create vectors
  LAFEM::DenseVector<DT_, IT_> vec_ref = factory.vector_q2_bubble();
  LAFEM::DenseVector<DT_, IT_> vec_rhs = vec_ref.clone(LAFEM::CloneMode::Layout);
  LAFEM::DenseVector<DT_, IT_> vec_sol = vec_ref.clone(LAFEM::CloneMode::Layout);

  // compute RHS and format solution
  matrix.apply(vec_rhs, vec_ref);
  vec_sol.format();

  LAFEM::NoneFilter<DT_, IT_> filter;

  // switch backend
  Runtime::set_preferred_backend(backend);

  // create PCG solver
  auto solver = Solver::new_pcg(matrix, filter);
  solver->set_max_iter(250);
  solver->set_tol_rel(Math::sqrt(Math::eps<DT_>()));
  solver->set_plot_mode(Solver::PlotMode::none);

  solver->init();

  TimeStamp stamp_1;
  solver->apply(vec_sol, vec_rhs);
  TimeStamp stamp_2;

  solver->done();

  Runtime::set_preferred_backend(PreferredBackend::generic);

  // compute error to reference solution
  vec_ref.axpy(vec_sol, vec_ref, -DT_(1));
  DT_ err = vec_ref.norm2();

  const unsigned long long n   = matrix.rows();
  const unsigned long long nze = matrix.used_elements();
  const unsigned long long k   = solver->get_num_iter();

  // compute FLOP/s and BYTE/s
  // 2x DOT in startup; 3x AXPY, 3x DOT, 1x SpMV in loop
  unsigned long long flop = 2ull*n + k * (6ull*n + nze);

  // 2x COPY + 2x DOT in startup
  // 1x COPY + 1x NORM2 + 2x DOT + 3x AXPY + 1x SpMV in loop
  unsigned long long byte = (8ull*n + k*(15ull*n + nze)) * sizeof(DT_) + k*(n+nze)*sizeof(IT_);

  double sec = stamp_2.elapsed(stamp_1);
  std::cout << "NumIters.: " << k << std::endl;
  std::cout << "RefError.: " << stringify_fp_sci(err) << std::endl;
  std::cout << "Runtime..: " << stamp_2.elapsed_string(stamp_1) << std::endl;
  std::cout << "GFLOP/s..: " << double(flop) / sec * 1E-9 << std::endl;
  std::cout << "GBYTE/s..: " << double(byte) / sec * 0.931E-9 << std::endl;
}

int main(int argc, char ** argv)
{
  Runtime::initialize(argc, argv);
#ifdef FEAT_HAVE_CUDA
#ifdef FEAT_HAVE_HALFMATH
  run<Half, unsigned int>(PreferredBackend::cuda);
#endif
  run<float, unsigned int>(PreferredBackend::cuda);
  run<double, unsigned int>(PreferredBackend::cuda);
#endif
  run<double, unsigned int>(PreferredBackend::generic);
  Runtime::finalize();
}

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;


template <typename SM_>
void run(PreferredBackend backend)
{
  Backend::set_preferred_backend(PreferredBackend::generic);
  typedef typename SM_::DataType DT_;
  typedef typename SM_::IndexType IT_;

  std::vector<IT_> num_of_nodes;
  num_of_nodes.push_back(4000);
  num_of_nodes.push_back(4000);

  // generate FE matrix A
  SparseMatrixBanded<DT_, IT_> bm(PointstarStructureFE::template value<DT_>(1, num_of_nodes));
  for (Index i(0) ; i < bm.get_elements_size().at(0) ; ++i)
    bm.val()[i] = DT_((i%4) + 1);
  SM_ sys;
  sys.convert(bm);
  Index size(sys.rows());

  Backend::set_preferred_backend(backend);
  std::cout<<backend<<" "<<SM_::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<"\n";
  std::cout<<"vector size: "<<size<<" used elements: "<<sys.used_elements()<<"\n";
  DenseVector<DT_, IT_> x(size);
  for (Index i (0) ; i < x.size() ; ++i)
    x(i, DT_(i%100) / DT_(100));
  DenseVector<DT_, IT_> r(size, DT_(4711));

  double flops(double(sys.used_elements()));
  flops *= 2;

  double bytes(double(sys.used_elements()));
  bytes *= double(sizeof(DT_));
  bytes += double(sys.used_elements() * sizeof(IT_));
  bytes += double(size * sizeof(DT_));

  auto func = [&] () { sys.apply(r, x); };
  run_bench(func, flops, bytes);

  MemoryPool::synchronize();
  std::cout<<"control norm: "<<x.norm2()<<"\n";
}

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
#ifdef FEAT_HAVE_CUDA
  run<SparseMatrixCSR<double, Index> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<double, unsigned int> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<float, Index> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<float, unsigned int> >(PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
  run<SparseMatrixCSR<Half, Index> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<Half, unsigned int> >(PreferredBackend::cuda);
#endif
#endif
  //run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, Index> >();
  //run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned int> >();
#ifdef FEAT_HAVE_MKL
  //run<Algo::MKL, SparseMatrixCSR<Mem::Main, double, unsigned long> >();
#endif
/*#ifdef FEAT_HAVE_CUDA
  run<Algo::CUDA, SparseMatrixBanded<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixBanded<Mem::CUDA, double, unsigned int> >();
#endif
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, unsigned int> >();*/
  return 0;
}

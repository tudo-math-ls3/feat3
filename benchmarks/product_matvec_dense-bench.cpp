// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template <typename DM_>
void run(PreferredBackend backend)
{
  Backend::set_preferred_backend(PreferredBackend::generic);
  typedef typename DM_::DataType DT_;

  //Index size(4096);
  Index size(16384);
  //Index size(32768);
  //Index size(65536);

  DenseMatrix<DT_, Index> x(size, size);
  DenseVector<DT_, Index> y(size);
  for (Index i(0) ; i < size ; ++i)
  {
    y.elements()[i]= DT_(-(i%100) * DT_(0.1) * (i%10));
  }
  for (Index i(0) ; i < size * size ; ++i)
  {
    x.elements()[i] = DT_((i%100) * DT_(0.5) * (i%10));
  }

  DenseVector<DT_, Index> r(size, 4711.);

  Backend::set_preferred_backend(backend);
  std::cout<<backend<<" "<<DM_::name()<<" "<<Type::Traits<DT_>::name()<<" rows/cols: " << size << "\n";

  double flops(double(x.used_elements()));
  flops *= 2;

  double bytes(double(x.used_elements()));
  bytes *= double(sizeof(DT_));
  bytes += double(size * sizeof(DT_));


  auto func = [&] () { x.apply(r, y); };
  run_bench(func, flops, bytes);

  MemoryPool::synchronize();
  std::cout<<"control norm: "<<r.norm2()<<"\n";
}

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
/*  run<DenseMatrix<Half, Index> >(PreferredBackend::generic);
  run<DenseMatrix<float, Index> >(PreferredBackend::generic);
  run<DenseMatrix<double, Index> >(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
  run<DenseMatrix<float, Index> >(PreferredBackend::mkl);
  run<DenseMatrix<double, Index> >(PreferredBackend::mkl);
#endif*/
#ifdef FEAT_HAVE_CUDA
#ifdef FEAT_HAVE_HALFMATH
  run<DenseMatrix<FEAT::Half, Index> >(PreferredBackend::cuda);
#endif
  run<DenseMatrix<float, Index> >(PreferredBackend::cuda);
  run<DenseMatrix<double, Index> >(PreferredBackend::cuda);
#endif
  return 0;
}

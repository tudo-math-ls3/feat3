// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template <typename VT_>
void run(PreferredBackend backend)
{
  Backend::set_preferred_backend(PreferredBackend::generic);
  typedef typename VT_::DataType DT_;
  typedef typename VT_::IndexType IT_;

  Index size(500000000ul);
  std::cout<<backend<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<"\n";
  std::cout<<"vector size: "<<size<<"\n";
  VT_ x(size, DT_(1.234));
  VT_ y(size, DT_(4711));

  Backend::set_preferred_backend(backend);
  double flops = double(size);
  flops *= 2;

  double bytes = double(size);
  bytes *= 2;
  bytes *= sizeof(DT_);

  auto func = [&] () { x.dot(y); };
  run_bench(func, flops, bytes);
}

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
#ifdef FEAT_HAVE_CUDA
#ifdef FEAT_HAVE_HALFMATH
  run<DenseVector<Half, Index> >(PreferredBackend::cuda);
#endif
  run<DenseVector<float, Index> >(PreferredBackend::cuda);
  run<DenseVector<double, Index> >(PreferredBackend::cuda);
#endif
/*  run<Algo::Generic, DenseVector<Mem::Main, float, Index> >();
  run<Algo::Generic, DenseVector<Mem::Main, double, Index> >();
#ifdef FEAT_HAVE_MKL
  run<Algo::MKL, DenseVector<Mem::Main, float, Index> >();
  run<Algo::MKL, DenseVector<Mem::Main, double, Index> >();
#endif*/
  return 0;
}

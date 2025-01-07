// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_matrix.hpp>
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

  //Index size(64);
  Index size(128);
  if (backend == PreferredBackend::cuda)
    size *= 8;
  size *= 16/sizeof(DT_);

  DenseMatrix<DT_, Index> x(size, size), y(size, size);
  for (Index i(0) ; i < size ; ++i)
  {
    for (Index j(0) ; j < size ; ++j)
    {
      x(i, j, DT_(i%100) * DT_(0.5) * DT_(i%10));
      y(i, j, -DT_(i%100) * DT_(0.1) * DT_(i%10));
    }
  }

  DM_ r(size, size, 4711.);

  DT_ alpha = DT_(1.);
  DT_ beta = DT_(0.);

  Backend::set_preferred_backend(backend);

  std::cout<<backend<<" "<<DM_::name()<<" "<<Type::Traits<DT_>::name()<<" rows/cols: " << size << "\n";

  double flops = 2. * double(x.rows() * x.rows() * r.columns() + r.columns());
  double bytes = 2. * double(x.rows() * x.columns() * y.columns() + r.columns() * r.rows());
  bytes *= sizeof(DT_);

  switch (backend)
  {
    case PreferredBackend::generic :
      {
        auto func = [&] () { Arch::ProductMatMat::dense_generic<DT_>(r.elements(), alpha, beta, x.elements(), y.elements(), r.elements(), r.rows(), r.columns(), x.columns()); };
        run_bench(func, flops, bytes);
        break;
      }

#ifdef FEAT_HAVE_MKL
    case PreferredBackend::mkl :
      {
        auto func = [&] () { Arch::ProductMatMat::dense_mkl(r.elements(), alpha, beta, x.elements(), y.elements(), r.elements(), r.rows(), r.columns(), x.columns()); };
        run_bench(func, flops, bytes);
        break;
      }
#endif

#ifdef FEAT_HAVE_CUDA
    case PreferredBackend::cuda :
      {
        auto func = [&] () { Arch::ProductMatMat::dense_cuda<DT_>(r.elements(), alpha, beta, x.elements(), y.elements(), r.elements(), r.rows(), r.columns(), x.columns()); };
        run_bench(func, flops, bytes);
        break;
      }
#endif

    default:
      throw InternalError("unsupported arch detected!");
  }

  MemoryPool::synchronize();
  std::cout<<"control norm: "<<double(r.norm_frobenius())<<"\n";
}

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  /*run<DenseMatrix<Half, Index> >(PreferredBackend::generic);
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
  //Util::cuda_reset_algos();
  run<DenseMatrix<float, Index> >(PreferredBackend::cuda);
  //Util::cuda_reset_algos();
  run<DenseMatrix<double, Index> >(PreferredBackend::cuda);
#endif
  return 0;
}

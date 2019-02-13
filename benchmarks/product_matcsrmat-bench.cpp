// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template <typename DT_, typename IT_>
void run(PreferredBackend backend)
{
  Runtime::set_preferred_backend(PreferredBackend::generic);

  Index size(128);
  if (backend == PreferredBackend::cuda)
    size *= 8;
  size *= 16/sizeof(DT_);
  size=8192;
  size *= 4;

  DT_ alpha(1.);
  DT_ beta(0.);

  DenseMatrix<DT_, Index> r(size, size, DT_(0)), y(size, size);
  for (Index i(0) ; i < size ; ++i)
  {
    for (Index j(0) ; j < size ; ++j)
    {
      y(i, j, DT_((1 + i*j) % 100));
    }
  }

  SparseMatrixFactory<DT_, IT_> x_fac(size, size);
  for (IT_ row(0) ; row < x_fac.rows() ; ++row)
  {
    for (IT_ col(0) ; col < x_fac.columns() ; ++col)
    {
      if(row == col)
      {
        x_fac.add(row, col, DT_(2));
      }
      else if((row == col+1) || (row+1 == col))
      {
        x_fac.add(row, col, DT_(-1));
      }
      else if((row == col+Index(Math::sqrt(double(size)))) || (row+Index(Math::sqrt(double(size))) == col))
      {
        x_fac.add(row, col, DT_(-1));
      }
    }
  }
  SparseMatrixCSR<DT_, IT_> x(x_fac.make_csr());

  Runtime::set_preferred_backend(backend);

  std::cout<<backend<<" "<<DenseMatrix<DT_, IT_>::name()<<" "<<SparseMatrixCSR<DT_, IT_>::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<" rows/cols: " << size << std::endl;

  //roughly 5 csr entries per row
  double flops = 2. * double(5 * x.rows() * r.columns());
  double bytes = 2. * double(5*2 * x.columns() * y.columns() + r.columns() * r.rows());
  bytes *= sizeof(DT_);

  switch (backend)
  {
    case PreferredBackend::generic :
      {
        auto func = [&] () { Arch::ProductMatMat::dsd_generic<DT_>(r.elements(), alpha, beta, x.val(), x.col_ind(), x.row_ptr(), x.used_elements(), y.elements(), r.rows(), r.columns(), x.columns()); };
        run_bench(func, flops, bytes);
        break;
      }

#ifdef FEAT_HAVE_CUDA
    case PreferredBackend::cuda :
      {
        auto func = [&] () { Arch::ProductMatMat::dsd_cuda<DT_>(r.elements(), alpha, beta, x.val(), x.col_ind(), x.row_ptr(), x.used_elements(), y.elements(), r.rows(), r.columns(), x.columns()); };
        run_bench(func, flops, bytes);
        break;
      }
#endif

    default:
      throw InternalError("unsupported arch detected!");
  }

  MemoryPool::synchronize();
  std::cout<<"control norm: "<<r.norm_frobenius()<<std::endl;
}

int main(int argc, char ** argv)
{
  Runtime::initialize(argc, argv);
  /*run<DenseMatrix<Half, Index> >(PreferredBackend::generic);
  run<DenseMatrix<float, Index> >(PreferredBackend::generic);
  run<DenseMatrix<double, Index> >(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
  run<DenseMatrix<float, Index> >(PreferredBackend::mkl);
  run<DenseMatrix<double, Index> >(PreferredBackend::mkl);
#endif*/
#ifdef FEAT_HAVE_CUDA
#ifdef FEAT_HAVE_HALFMATH
  run<FEAT::Half, unsigned int>(PreferredBackend::cuda);
#endif
  run<float, unsigned int>(PreferredBackend::cuda);
  run<double, unsigned int>(PreferredBackend::cuda);
#endif
  Runtime::finalize();
}

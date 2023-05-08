// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template <typename DT_, typename IT_>
void run(PreferredBackend backend, const String filename, bool transpose)
{
  DT_ alpha(1.);
  DT_ beta(0.);

  Backend::set_preferred_backend(PreferredBackend::generic);

  SparseMatrixCSR<DT_, IT_> x(FileMode::fm_csr, filename);
  if (transpose)
    x=x.transpose();
  std::cout<<"csr loaded "<<x.rows()<< " " <<x.columns()<<std::endl;


  DenseMatrix<DT_, Index> r(x.rows(), 5000, DT_(0)), y(x.columns(), 5000, DT_(2.345));
  for (Index i(0) ; i < y.rows() ; ++i)
  {
    for (Index j(0) ; j < y.columns() ; ++j)
    {
      y(i, j, DT_((1 + i*j) % 100));
    }
  }

  Backend::set_preferred_backend(backend);


  std::cout<<backend<<" "<<DenseMatrix<DT_, IT_>::name()<<" "<<SparseMatrixCSR<DT_, IT_>::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<" rows/cols/dense: " << x.rows()<<" "<<x.columns()<< " "<<r.columns() << std::endl;

  double flops = 2. * double(x.used_elements() * y.columns());
  double bytes = 2. * double(2 * x.used_elements() * x.columns() * y.columns() + r.columns() * r.rows());
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
  if (!(argc == 2 || argc == 3))
  {
    throw InternalError("this benchmarks need the path to a csr (binary) matrix file as its single command line parameter. A following parameter t for transpose is optional");
  }
  String filename = String(argv[1]);
  bool transpose(false);
  if (argc == 3)
  {
    if (String(argv[2]) == "t")
      transpose=true;
    else if (String(argv[2]) == "n")
      transpose=false;
    else
      throw InternalError("second parameter " + String(argv[2]) + " not known! Only t or n for transpose option are accepted.");
  }
  /*run<DenseMatrix<Half, Index> >(PreferredBackend::generic, filename, transpose);
  run<DenseMatrix<float, Index> >(PreferredBackend::generic, filename, transpose);
  run<DenseMatrix<double, Index> >(PreferredBackend::generic, filename, transpose);
#ifdef FEAT_HAVE_MKL
  run<DenseMatrix<float, Index> >(PreferredBackend::mkl, filename, transpose);
  run<DenseMatrix<double, Index> >(PreferredBackend::mkl, filename, transpose);
#endif*/
#ifdef FEAT_HAVE_CUDA
#ifdef FEAT_HAVE_HALFMATH
  run<FEAT::Half, unsigned int>(PreferredBackend::cuda, filename, transpose);
#endif
  run<float, unsigned int>(PreferredBackend::cuda, filename, transpose);
  run<double, unsigned int>(PreferredBackend::cuda, filename, transpose);
#endif
  Runtime::finalize();
}

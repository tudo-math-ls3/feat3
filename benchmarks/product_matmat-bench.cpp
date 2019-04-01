// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template<typename Algo_, typename DT_>
struct ProductMatMatBench;

template<typename DT_>
struct ProductMatMatBench<Algo::Generic, DT_>
{
  static void f(DenseMatrix<Mem::Main, DT_, Index> & r, const DenseMatrix<Mem::Main, DT_, Index> & x,
    DenseMatrix<Mem::Main, DT_, Index> & y)
  {
    Arch::ProductMatMat<Mem::Main>::dense_generic(r.elements(), x.elements(),
        y.elements(), r.rows(), r.columns(), x.columns());
  }
};

template<typename DT_>
struct ProductMatMatBench<Algo::MKL,  DT_>
{
  static void f(DenseMatrix<Mem::Main, DT_, Index> & r, const DenseMatrix<Mem::Main, DT_, Index> & x,
    DenseMatrix<Mem::Main, DT_, Index> & y)
  {
    Arch::ProductMatMat<Mem::Main>::dense_mkl(r.elements(), x.elements(),
        y.elements(), r.rows(), r.columns(), x.columns());
  }
};

template<typename DT_>
struct ProductMatMatBench<Algo::CUDA,  DT_>
{
  static void f(DenseMatrix<Mem::CUDA, DT_, Index> & r, const DenseMatrix<Mem::CUDA, DT_, Index> & x,
    DenseMatrix<Mem::CUDA, DT_, Index> & y)
  {
    Arch::ProductMatMat<Mem::CUDA>::dense(r.elements(), x.elements(),
        y.elements(), r.rows(), r.columns(), x.columns());
  }
};


template <typename Algo_, typename DM_>
void run()
{
  typedef typename DM_::DataType DT_;
  typedef typename DM_::MemType Mem_;

  Index size(2000);
  DenseMatrix<Mem::Main, DT_, Index> x_local(size, size), y_local(size, size);
  for (Index i(0) ; i < size ; ++i)
  {
    for (Index j(0) ; j < size ; ++j)
    {
      x_local(i, j, DT_(DT_(i%100) * DT_(0.5) * DT_(i%10)));
      y_local(i, j, DT_(-DT_(i%100) * DT_(0.1) * DT_(i%10)));
    }
  }

  DM_ x;
  x.convert(x_local);
  DM_ y;
  y.convert(x_local);
  DM_ r(size, size, 4711);

  std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<DM_::name()<<" "<<Type::Traits<DT_>::name()<<" rows/cols: " << size << std::endl;

  double flops = 2. * x.rows() * x.rows() * r.columns() + r.columns();
  double bytes = 2. * x.rows() * x.columns() * y.columns() + r.columns() * r.rows();
  bytes *= sizeof(DT_);


  auto func = [&] () { ProductMatMatBench<Algo_, DT_>::f(r, x, y); };
  run_bench<Mem_>(func, flops, bytes);

  std::cout<<"control norm: "<<r.norm_frobenius()<<std::endl;
}

int main(int argc, char ** argv)
{
  Runtime::initialise(argc, argv);
  run<Algo::Generic, DenseMatrix<Mem::Main, float, Index> >();
  run<Algo::Generic, DenseMatrix<Mem::Main, double, Index> >();
#ifdef FEAT_HAVE_MKL
  run<Algo::MKL, DenseMatrix<Mem::Main, float, Index> >();
  run<Algo::MKL, DenseMatrix<Mem::Main, double, Index> >();
#endif
#ifdef FEAT_HAVE_CUDA
  run<Algo::CUDA, DenseMatrix<Mem::CUDA, float, Index> >();
  run<Algo::CUDA, DenseMatrix<Mem::CUDA, double, Index> >();
#endif
  Runtime::finalise();
}

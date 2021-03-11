// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template<typename Algo_, typename DT_, typename IT_>
struct ProductMatVecBench;

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::Generic, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    DenseMatrix<Mem::Main, DT_, IT_> & A)
  {
    Arch::Apply<Mem::Main>::dense_generic(x.elements(), DT_(1), DT_(0), b.elements(), A.elements(), b.elements(), A.rows(), A.columns());
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::MKL, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    DenseMatrix<Mem::Main, DT_, IT_> & A)
  {
    Arch::Apply<Mem::Main>::dense(x.elements(), DT_(1), DT_(0), b.elements(), A.elements(), b.elements(), A.rows(), A.columns());
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::CUDA, DT_, IT_>
{
  static void f(DenseVector<Mem::CUDA, DT_, IT_> & x, const DenseVector<Mem::CUDA, DT_, IT_> & b,
    DenseMatrix<Mem::CUDA, DT_, IT_> & A)
  {
    Arch::Apply<Mem::CUDA>::dense(x.elements(), DT_(1), DT_(0), b.elements(), A.elements(), b.elements(), A.rows(), A.columns());
  }
};

template <typename Mem_, typename Algo_, typename DT_>
void run()
{
  using IT_ = Index;
  DenseMatrix<Mem::Main, DT_, IT_> sys_main(4096, 4096);
  for (Index i(0) ; i < sys_main.used_elements() ; ++i)
  {
    sys_main.elements()[i] = DT_(sys_main.used_elements() / (i % sys_main.rows() + 1));
  }
  DenseMatrix<Mem_, DT_, IT_> sys;
  sys.convert(sys_main);

  Index size(sys.rows());
  std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<Type::Traits<DT_>::name()<<std::endl;
  std::cout<<"vector size: "<<size<<" used elements: "<<sys.used_elements()<<std::endl;
  DenseVector<Mem::Main, DT_, IT_> bmain(size);
  for (Index i (0) ; i < bmain.size() ; ++i)
    bmain(i, DT_(i%100) / DT_(100));
  DenseVector<Mem_, DT_, IT_> b;
  b.convert(bmain);
  DenseVector<Mem_, DT_, IT_> x(size, DT_(4711));

  double flops(double(sys.used_elements()));
  flops *= 2;

  double bytes(double(sys.used_elements()));
  bytes *= double(sizeof(DT_));
  bytes += double(sys.used_elements() * sizeof(IT_));
  bytes += double(size * sizeof(DT_));

  auto func = [&] () { ProductMatVecBench<Algo_, DT_, IT_>::f(x, b, sys); };
  run_bench<Mem_>(func, flops, bytes);

  std::cout<<"control norm: "<<x.norm2()<<std::endl;
}

int main(int argc, char ** argv)
{
  Runtime::initialize(argc, argv);
#ifdef FEAT_HAVE_CUDA
  run<Mem::CUDA, Algo::CUDA, double>();
  run<Mem::CUDA, Algo::CUDA, float>();
#endif
#ifdef FEAT_HAVE_MKL
  run<Mem::Main, Algo::MKL, double>();
  run<Mem::Main, Algo::MKL, float>();
#endif
  run<Mem::Main, Algo::Generic, double>();
  run<Mem::Main, Algo::Generic, float>();
  Runtime::finalize();
}

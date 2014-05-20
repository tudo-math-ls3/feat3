#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>

#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::Benchmark;

template <typename Algo_, typename VT_>
void run()
{
  typedef typename VT_::DataType DT_;
  typedef typename VT_::IndexType IT_;
  typedef typename VT_::MemType Mem_;

  Index size(5000000ul);
  std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<std::endl;
  std::cout<<"vector size: "<<size<<std::endl;
  DenseVector<Mem_, DT_, IT_> x(size, DT_(1.234));
  DenseVector<Mem_, DT_, IT_> y(size, DT_(4711));
  DT_ s(23);

  double flops(size);
  flops *= 2;

  double bytes(size);
  bytes *= 3;
  bytes *= sizeof(DT_);

  auto func = [&] () { y.template axpy<Algo_>(x, y, s); };
  run_bench<Mem_>(func, flops, bytes);
}

int main(int argc, char ** argv)
{
  run<Algo::CUDA, DenseVector<Mem::CUDA, double, Index> >();
  run<Algo::Generic, DenseVector<Mem::Main, double, Index> >();
  run<Algo::MKL, DenseVector<Mem::Main, double, Index> >();
}

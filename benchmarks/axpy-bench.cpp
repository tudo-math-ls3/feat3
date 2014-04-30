#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/type_traits.hpp>

#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;

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

  std::vector<double> times;
  const Index iters(25);
  for (Index i(0) ; i < 25 ; ++i)
  {
    TimeStamp at, bt;
    at.stamp();
    for (Index j(0) ; j < iters ; ++j)
    {
      y.template axpy<Algo_>(x, y, s);
    }
    MemoryPool<Mem_>::synchronize();
    bt.stamp();
    times.push_back(bt.elapsed(at));
  }

  double mean(0);
  for (auto & time : times)
    mean += time;
  mean /= DT_(times.size());
  std::cout<<"TOE: "<<std::fixed<<mean<<std::endl;
  double flops(size);
  flops *= 2;
  flops *= iters;
  flops /= mean;
  flops /= 1000; // kilo
  flops /= 1000; // mega
  flops /= 1000; // giga
  std::cout<<"GFlop/s: "<<flops<<std::endl;
  double bytes(size);
  bytes *= 3;
  bytes *= sizeof(DT_);
  bytes *= iters;
  bytes /= mean;
  bytes /= 1024; // kilo
  bytes /= 1024; // mega
  bytes /= 1024; // giga
  std::cout<<"GByte/s: "<<bytes<<std::endl;
  std::cout<<"=============================================="<<std::endl;
}

int main(int argc, char ** argv)
{
  run<Algo::CUDA, DenseVector<Mem::CUDA, double, Index> >();
  run<Algo::Generic, DenseVector<Mem::Main, double, Index> >();
  run<Algo::MKL, DenseVector<Mem::Main, double, Index> >();
}

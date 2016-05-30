#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template<typename Algo_, typename DT_, typename IT_>
class AxpyBench;

template<typename DT_, typename IT_>
class AxpyBench<Algo::Generic, DT_, IT_>
{
  public:
  static void f(const DenseVector<Mem::Main, DT_, IT_> & x, DenseVector<Mem::Main, DT_, IT_> & y, DT_ s)
  {
    Arch::Axpy<Mem::Main>::dv_generic(y.elements(), s, x.elements(), y.elements(), x.size());
  }
};

template<typename DT_, typename IT_>
class AxpyBench<Algo::MKL, DT_, IT_>
{
  public:
  static void f(const DenseVector<Mem::Main, DT_, IT_> & x, DenseVector<Mem::Main, DT_, IT_> & y, DT_ s)
  {
    Arch::Axpy<Mem::Main>::dv_mkl(y.elements(), s, x.elements(), y.elements(), x.size());
  }
};

template<typename DT_, typename IT_>
class AxpyBench<Algo::CUDA, DT_, IT_>
{
  public:
  static void f(const DenseVector<Mem::CUDA, DT_, IT_> & x, DenseVector<Mem::CUDA, DT_, IT_> & y, DT_ s)
  {
    Arch::Axpy<Mem::CUDA>::dv(y.elements(), s, x.elements(), y.elements(), x.size());
  }
};


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

  double flops = double(size);
  flops *= 2;

  double bytes = double(size);
  bytes *= 3;
  bytes *= sizeof(DT_);

  auto func = [&] () { AxpyBench<Algo_, DT_, IT_>::f(x, y, s); };
  run_bench<Mem_>(func, flops, bytes);
}

int main(int argc, char ** argv)
{
  Runtime::initialise(argc, argv);
#ifdef FEAT_BACKENDS_CUDA
  run<Algo::CUDA, DenseVector<Mem::CUDA, double, Index> >();
#endif
  run<Algo::Generic, DenseVector<Mem::Main, double, Index> >();
#ifdef FEAT_BACKENDS_MKL
  run<Algo::MKL, DenseVector<Mem::Main, double, Index> >();
#endif
  Runtime::finalise();
}

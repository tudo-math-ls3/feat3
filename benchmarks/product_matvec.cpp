#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/type_traits.hpp>

#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename Algo_, typename SM_>
void run()
{
  typedef typename SM_::DataType DT_;
  typedef typename SM_::IndexType IT_;
  typedef typename SM_::MemType Mem_;
  typedef PointstarFactoryFE<DT_> PSF_;

  PSF_ factory(2500);
  SM_ sys(factory.matrix_csr());
  Index size(sys.rows());
  std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<SM_::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<std::endl;
  std::cout<<"vector size: "<<size<<" used elements: "<<sys.used_elements()<<std::endl;
  //DenseVector<Mem_, DT_, IT_> b(factory.vector_q2_bubble());
  DenseVector<Mem_, DT_, IT_> b(size, DT_(1.234));
  DenseVector<Mem_, DT_, IT_> x(size);

  std::vector<double> times;
  Index iters(25);
  for (Index i(0) ; i < 10 ; ++i)
  {
    TimeStamp at, bt;
    at.stamp();
    for (Index i(0) ; i < iters ; ++i)
    {
      sys.template apply<Algo_>(x, b);
    }
    MemoryPool<Mem_>::synchronize();
    bt.stamp();
    times.push_back(bt.elapsed(at));
  }

  double mean(0);
  for (auto & time : times)
    mean += time;
  mean /= DT_(times.size());
  std::cout<<"control norm: "<<x.template norm2<Algo_>()<<std::endl;
  std::cout<<"TOE: "<<std::fixed<<mean<<std::endl;
  double flops(sys.used_elements());
  flops *= 2;
  flops *= iters;
  flops /= mean;
  flops /= 1000; // kilo
  flops /= 1000; // mega
  flops /= 1000; // giga
  std::cout<<"GFlop/s: "<<flops<<std::endl;
  std::cout<<"=============================================="<<std::endl;
}

int main(int argc, char ** argv)
{
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, unsigned int> >();
  run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, unsigned int> >();
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned int> >();
  run<Algo::MKL, SparseMatrixCSR<Mem::Main, double> >();
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, unsigned int> >();
}

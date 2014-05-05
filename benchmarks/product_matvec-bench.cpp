#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
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

  std::vector<IT_> num_of_nodes;
  num_of_nodes.push_back(1200);
  num_of_nodes.push_back(1200);

  // generate FE matrix A
  SparseMatrixBanded<Mem::Main, DT_, IT_> bm(PointstarStructureFE<Algo::Generic>::template value<DT_>(1, num_of_nodes));
  for (Index i(0) ; i < bm.get_elements_size().at(0) ; ++i)
    bm.val()[i] = DT_((i%4) + 1);
  SM_ sys;
  sys.convert(bm);
  Index size(sys.rows());
  std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<SM_::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<std::endl;
  std::cout<<"vector size: "<<size<<" used elements: "<<sys.used_elements()<<std::endl;
  DenseVector<Mem_, DT_, IT_> b(size, DT_(1.234));
  DenseVector<Mem_, DT_, IT_> x(size, DT_(4711));

  std::vector<double> times;
  Index iters(25);
  for (Index i(0) ; i < 10 ; ++i)
  {
    TimeStamp at, bt;
    at.stamp();
    for (Index j(0) ; j < iters ; ++j)
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
  double bytes(sys.used_elements());
  bytes *= 2;
  bytes += size;
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
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, unsigned int> >();
  //run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, unsigned int> >();
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned int> >();
  run<Algo::MKL, SparseMatrixCSR<Mem::Main, double> >();
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, unsigned int> >();
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, unsigned int> >();
}

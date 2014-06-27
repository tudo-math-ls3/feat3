#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>

#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::Benchmark;

template <typename Algo_, typename SM_>
void run()
{
  typedef typename SM_::DataType DT_;
  typedef typename SM_::IndexType IT_;
  typedef typename SM_::MemType Mem_;

  std::vector<IT_> num_of_nodes;
  num_of_nodes.push_back(1300);
  num_of_nodes.push_back(1300);

  // generate FE matrix A
  SparseMatrixBanded<Mem::Main, DT_, IT_> bm(PointstarStructureFE<Algo::Generic>::template value<DT_>(1, num_of_nodes));
  for (Index i(0) ; i < bm.get_elements_size().at(0) ; ++i)
    bm.val()[i] = DT_((i%4) + 1);
  SM_ sys;
  sys.convert(bm);
  Index size(sys.rows());
  std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<SM_::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<std::endl;
  std::cout<<"vector size: "<<size<<" used elements: "<<sys.used_elements()<<std::endl;
  DenseVector<Mem_, DT_, IT_> b(size);
  for (Index i (0) ; i < b.size() ; ++i)
    b(i, DT_(i%100) / DT_(100));
  DenseVector<Mem_, DT_, IT_> x(size, DT_(4711));

  double flops(double(sys.used_elements()));
  flops *= 2;

  double bytes(double(sys.used_elements()));
  bytes *= 2;
  bytes += size;
  bytes *= sizeof(DT_);

  auto func = [&] () { sys.template apply<Algo_>(x, b); };
  run_bench<Mem_>(func, flops, bytes);

  std::cout<<"control norm: "<<x.template norm2<Algo_>()<<std::endl;
}

int main(int /*argc*/, char ** /*argv*/)
{
#ifdef FEAST_BACKENDS_CUDA
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, unsigned int> >();
  //run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, unsigned int> >();
#endif
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned int> >();
#ifdef FEAST_BACKENDS_MKL
  run<Algo::MKL, SparseMatrixCSR<Mem::Main, double> >();
#endif
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, unsigned int> >();
#ifdef FEAST_BACKENDS_CUDA
  run<Algo::CUDA, SparseMatrixBanded<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixBanded<Mem::CUDA, double, unsigned int> >();
#endif
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, unsigned int> >();
}

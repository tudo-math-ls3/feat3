// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template<typename Algo_, SparseLayoutId, typename DT_, typename IT_>
struct ProductMatVecBench;

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::Generic, SparseLayoutId::lt_csr, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    SparseMatrixCSR<Mem::Main, DT_, IT_> & A)
  {
    Arch::Apply<Mem::Main>::csr_generic(x.elements(), DT_(1), b.elements(), DT_(0), x.elements(), A.val(), A.col_ind(), A.row_ptr(),
        A.rows(), A.columns(), A.used_elements(), false);
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::MKL, SparseLayoutId::lt_csr, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    SparseMatrixCSR<Mem::Main, DT_, IT_> & A)
  {
    Arch::Apply<Mem::Main>::csr_mkl(x.elements(), DT_(1), b.elements(), DT_(0), x.elements(), A.val(), A.col_ind(), A.row_ptr(),
        A.rows(), A.columns(), A.used_elements(), false);
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::CUDA, SparseLayoutId::lt_csr, DT_, IT_>
{
  static void f(DenseVector<Mem::CUDA, DT_, IT_> & x, const DenseVector<Mem::CUDA, DT_, IT_> & b,
    SparseMatrixCSR<Mem::CUDA, DT_, IT_> & A)
  {
    Arch::Apply<Mem::CUDA>::csr(x.elements(), DT_(1), b.elements(), DT_(0), x.elements(), A.val(), A.col_ind(), A.row_ptr(),
        A.rows(), A.columns(), A.used_elements(), false);
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::Generic, SparseLayoutId::lt_ell, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    SparseMatrixELL<Mem::Main, DT_, IT_> & A)
  {
    Arch::Apply<Mem::Main>::ell_generic(x.elements(), DT_(1), b.elements(), DT_(0), x.elements(), A.val(), A.col_ind(), A.cs(), A.cl(),
        A.C(), A.rows());
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::CUDA, SparseLayoutId::lt_ell, DT_, IT_>
{
  static void f(DenseVector<Mem::CUDA, DT_, IT_> & x, const DenseVector<Mem::CUDA, DT_, IT_> & b,
    SparseMatrixELL<Mem::CUDA, DT_, IT_> & A)
  {
    Arch::Apply<Mem::CUDA>::ell(x.elements(), DT_(1), b.elements(), DT_(0), x.elements(), A.val(), A.col_ind(), A.cs(), A.cl(),
        A.C(), A.rows());
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::Generic, SparseLayoutId::lt_banded, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    SparseMatrixBanded<Mem::Main, DT_, IT_> & A)
  {
    Arch::Apply<Mem::Main>::banded_generic(x.elements(), DT_(1), b.elements(), DT_(0), x.elements(), A.val(), A.offsets(), A.num_of_offsets(), A.rows(), A.columns());
  }
};

template<typename DT_, typename IT_>
struct ProductMatVecBench<Algo::CUDA, SparseLayoutId::lt_banded, DT_, IT_>
{
  static void f(DenseVector<Mem::CUDA, DT_, IT_> & x, const DenseVector<Mem::CUDA, DT_, IT_> & b,
    SparseMatrixBanded<Mem::CUDA, DT_, IT_> & A)
  {
    Arch::Apply<Mem::CUDA>::banded(x.elements(), DT_(1), b.elements(), DT_(0), x.elements(), A.val(), A.offsets(), A.num_of_offsets(), A.rows(), A.columns());
  }
};


template <typename Algo_, typename SM_>
void run()
{
  typedef typename SM_::DataType DT_;
  typedef typename SM_::IndexType IT_;
  typedef typename SM_::MemType Mem_;

  std::vector<IT_> num_of_nodes;
  num_of_nodes.push_back(2000);
  num_of_nodes.push_back(2000);

  // generate FE matrix A
  SparseMatrixBanded<Mem::Main, DT_, IT_> bm(PointstarStructureFE::template value<DT_>(1, num_of_nodes));
  for (Index i(0) ; i < bm.get_elements_size().at(0) ; ++i)
    bm.val()[i] = DT_((i%4) + 1);
  SM_ sys;
  sys.convert(bm);
  Index size(sys.rows());
  std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<SM_::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<std::endl;
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

  auto func = [&] () { ProductMatVecBench<Algo_, SM_::layout_id, DT_, IT_>::f(x, b, sys); };
  run_bench<Mem_>(func, flops, bytes);

  std::cout<<"control norm: "<<x.norm2()<<std::endl;
}

int main(int argc, char ** argv)
{
  Runtime::initialize(argc, argv);
#ifdef FEAT_HAVE_CUDA
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double, unsigned int> >();
  //run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, unsigned int> >();
#endif
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned int> >();
#ifdef FEAT_HAVE_MKL
  run<Algo::MKL, SparseMatrixCSR<Mem::Main, double, unsigned long> >();
#endif
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixELL<Mem::Main, double, unsigned int> >();
#ifdef FEAT_HAVE_CUDA
  run<Algo::CUDA, SparseMatrixBanded<Mem::CUDA, double, Index> >();
  run<Algo::CUDA, SparseMatrixBanded<Mem::CUDA, double, unsigned int> >();
#endif
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, Index> >();
  run<Algo::Generic, SparseMatrixBanded<Mem::Main, double, unsigned int> >();
  Runtime::finalize();
}

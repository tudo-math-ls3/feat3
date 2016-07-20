#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/util/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template<typename Algo_, SparseLayoutId, typename DT_, typename IT_>
struct BlockProductMatVecBench;

template<typename DT_, typename IT_>
struct BlockProductMatVecBench<Algo::Generic, SparseLayoutId::lt_csr, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    SparseMatrixCSR<Mem::Main, DT_, IT_> & A)
  {
    Arch::ProductMatVec<Mem::Main>::csr_generic(x.elements(), A.val(), A.col_ind(), A.row_ptr(),
                                              b.elements(), A.rows(), A.columns(), A.used_elements());
  }

  template <int BlockHeight_, int BlockWidth_>
  static void f(DenseVectorBlocked<Mem::Main, DT_, IT_, BlockHeight_> & x, const DenseVectorBlocked<Mem::Main, DT_, IT_, BlockWidth_> & b,
    SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_> & A)
  {
    Arch::ProductMatVec<Mem::Main>::template csrb_generic<DT_, IT_, BlockHeight_, BlockWidth_>(x.template elements<Perspective::pod>(), A.template val<Perspective::pod>(), A.col_ind(), A.row_ptr(),
                                              b.template elements<Perspective::pod>(), A.rows(), A.columns(), A.used_elements());
  }
};

template<typename DT_, typename IT_>
struct BlockProductMatVecBench<Algo::MKL, SparseLayoutId::lt_csr, DT_, IT_>
{
  static void f(DenseVector<Mem::Main, DT_, IT_> & x, const DenseVector<Mem::Main, DT_, IT_> & b,
    SparseMatrixCSR<Mem::Main, DT_, IT_> & A)
  {
    Arch::ProductMatVec<Mem::Main>::csr_mkl(x.elements(), A.val(), A.col_ind(), A.row_ptr(),
                                              b.elements(), A.rows(), A.columns(), A.used_elements());
  }

  template <int BlockHeight_, int BlockWidth_>
  static void f(DenseVectorBlocked<Mem::Main, DT_, IT_, BlockHeight_> & x, const DenseVectorBlocked<Mem::Main, DT_, IT_, BlockWidth_> & b,
    SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_> & A)
  {
    Arch::ProductMatVec<Mem::Main>::csrb_mkl(x.template elements<Perspective::pod>(), A.template val<Perspective::pod>(), A.col_ind(), A.row_ptr(),
                                              b.template elements<Perspective::pod>(), A.rows(), A.columns(), A.used_elements(), BlockHeight_);
  }
};

template<typename DT_, typename IT_>
struct BlockProductMatVecBench<Algo::CUDA, SparseLayoutId::lt_csr, DT_, IT_>
{
  static void f(DenseVector<Mem::CUDA, DT_, IT_> & x, const DenseVector<Mem::CUDA, DT_, IT_> & b,
    SparseMatrixCSR<Mem::CUDA, DT_, IT_> & A)
  {
    Arch::ProductMatVec<Mem::CUDA>::csr(x.elements(), A.val(), A.col_ind(), A.row_ptr(),
                                              b.elements(), A.rows(), A.columns(), A.used_elements());
  }

  template <int BlockHeight_, int BlockWidth_>
  static void f(DenseVectorBlocked<Mem::CUDA, DT_, IT_, BlockHeight_> & x, const DenseVectorBlocked<Mem::CUDA, DT_, IT_, BlockWidth_> & b,
    SparseMatrixBCSR<Mem::CUDA, DT_, IT_, BlockHeight_, BlockWidth_> & A)
  {
    Arch::ProductMatVec<Mem::CUDA>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(x.template elements<Perspective::pod>(), A.template val<Perspective::pod>(), A.col_ind(), A.row_ptr(),
                                              b.template elements<Perspective::pod>(), A.rows(), A.columns(), A.used_elements());
  }
};


template <typename Algo_, typename SM_, int Blocksize_>
void run()
{
  typedef typename SM_::DataType DT_;
  typedef typename SM_::IndexType IT_;
  typedef typename SM_::MemType Mem_;

  const Index m = 600;
  const Index d = 2;
  PointstarFactoryFD<DT_, IT_> psf(m, d);

  // create 5-point star CSR matrix
  SparseMatrixCSR<Mem::Main, DT_, IT_> init_mat(psf.matrix_csr());
  auto layout = init_mat.layout();

  SparseMatrixBCSR<Mem::Main, DT_, IT_, Blocksize_, Blocksize_> sys_main(layout);
  for (Index i(0) ; i < sys_main.used_elements() ; ++i)
  {
    sys_main.val()[i] = init_mat.val()[i];
  }
  init_mat.clear();

  {
    SparseMatrixBCSR<Mem_, DT_, IT_, Blocksize_, Blocksize_> sys;
    sys.convert(sys_main);


    Index size(sys.rows());
    Index used_elements(sys.template used_elements<Perspective::pod>());
    std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<sys.name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<" Blocksize: " << stringify(Blocksize_) << std::endl;
    std::cout<<"vector size: "<<size<<", Blocksize: " << stringify(Blocksize_) << ", used elements: "<<used_elements<<std::endl;
    DenseVectorBlocked<Mem::Main, DT_, IT_, Blocksize_> bmain(size);
    for (Index i (0) ; i < bmain.size() ; ++i)
    {
      auto temp = bmain(i);
      temp = DT_(i%100) / DT_(100);
      bmain(i, temp);
    }
    DenseVectorBlocked<Mem_, DT_, IT_, Blocksize_> b;
    b.convert(bmain);
    DenseVectorBlocked<Mem_, DT_, IT_, Blocksize_> x(size, DT_(4711));

    double flops = double(used_elements);
    flops *= 2;

    double bytes = double(used_elements);
    bytes *= sizeof(DT_);
    bytes += used_elements * sizeof(IT_);
    bytes += size * Blocksize_ * sizeof(DT_);

    auto func = [&] () { BlockProductMatVecBench<Algo_, SM_::layout_id, DT_, IT_>::f(x, b, sys); };
    run_bench<Mem_>(func, flops, bytes);

    std::cout<<"control norm: "<<x.norm2()<<std::endl;
  }

  {
    SparseMatrixCSR<Mem_, DT_, IT_> sys;
    sys.convert(sys_main);


    Index size(sys.rows());
    Index used_elements(sys.template used_elements<Perspective::pod>());
    std::cout<<Mem_::name()<<" "<<Algo_::name()<<" "<<sys.name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<< std::endl;
    std::cout<<"vector size: "<<size<<", used elements: "<<used_elements<<std::endl;
    DenseVector<Mem::Main, DT_, IT_> bmain(size);
    for (Index i (0) ; i < bmain.size() ; ++i)
    {
      bmain(i, DT_(i%100) / DT_(100));
    }
    DenseVector<Mem_, DT_, IT_> b;
    b.convert(bmain);
    DenseVector<Mem_, DT_, IT_> x(size, DT_(4711));

    double flops = double(used_elements);
    flops *= 2;

    double bytes = double(used_elements);
    bytes *= sizeof(DT_);
    bytes += used_elements * sizeof(IT_);
    bytes += size * sizeof(DT_);

    auto func = [&] () { BlockProductMatVecBench<Algo_, SM_::layout_id, DT_, IT_>::f(x, b, sys); };
    run_bench<Mem_>(func, flops, bytes);

    std::cout<<"control norm: "<<x.norm2()<<std::endl;
  }
}

int main(int argc, char ** argv)
{
  Runtime::initialise(argc, argv);
#ifdef FEAT_HAVE_CUDA
  run<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double, unsigned int>, 2>();
#endif
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, Index>, 2>();
  run<Algo::Generic, SparseMatrixCSR<Mem::Main, double, unsigned int> , 2>();
#ifdef FEAT_HAVE_MKL
  run<Algo::MKL, SparseMatrixCSR<Mem::Main, double, unsigned long>, 2>();
#endif
  Runtime::finalise();
}

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;

template <typename SM_, int dim, int Blockheight_, int Blockwidth_ = Blockheight_>
void run(PreferredBackend backend)
{
  Backend::set_preferred_backend(PreferredBackend::generic);
  typedef typename SM_::DataType DT_;
  typedef typename SM_::IndexType IT_;

  const Index m = dim == 2 ? 600 : 100;
  const Index d = dim;
  PointstarFactoryFD<DT_, IT_> psf(m, d);

  // create 5-point star CSR matrix
  SparseMatrixCSR<DT_, IT_> init_mat(psf.matrix_csr());
  auto layout = init_mat.layout();

  SparseMatrixBCSR<DT_, IT_, Blockheight_, Blockwidth_> sys_main(layout);
  for (Index i(0) ; i < sys_main.used_elements() ; ++i)
  {
    sys_main.val()[i] = init_mat.val()[i];
  }
  init_mat.clear();

  {
    Backend::set_preferred_backend(PreferredBackend::generic);
    SparseMatrixBCSR<DT_, IT_, Blockheight_, Blockwidth_> sys;
    sys.convert(sys_main);


    Index size(sys.rows());
    Index used_elements(sys.template used_elements<Perspective::pod>());
    std::cout<<backend<<" "<<sys.name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<< " dimension " << d <<" Blocksize: " << stringify(Blockheight_) << ", " << stringify(Blockwidth_) << "\n";
    std::cout<<"vector size: "<<size<<", Blockheight: " << stringify(Blockheight_) << ", Blockheight: " << stringify(Blockheight_);
    std::cout << ", used elements: "<<used_elements<<"\n";
    DenseVectorBlocked<DT_, IT_, Blockheight_> b(size);
    for (Index i (0) ; i < b.size() ; ++i)
    {
      typename std::remove_const<decltype(b(i))>::type  temp;
      temp = DT_(i%100) / DT_(100);
      b(i, temp);
    }
    DenseVectorBlocked<DT_, IT_, Blockwidth_> x(size, DT_(4711));

    double flops = double(used_elements);
    flops *= 2;

    double bytes = double(used_elements);
    bytes *= DT_(sizeof(DT_));
    bytes += DT_(used_elements * sizeof(IT_));
    bytes += DT_(size * Blockheight_ * sizeof(DT_));

    Backend::set_preferred_backend(backend);
    auto func = [&] () { sys.apply(b, x); };
    run_bench(func, flops, bytes);

    std::cout<<"control norm: "<<x.norm2()<<"\n";
  }

  if constexpr(Blockheight_ == Blockwidth_)
  {
    Backend::set_preferred_backend(PreferredBackend::generic);
    SparseMatrixCSR<DT_, IT_> sys;
    sys.convert(sys_main);


    Index size(sys.rows());
    Index used_elements(sys.template used_elements<Perspective::pod>());
    std::cout<<backend<<" "<<sys.name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<< " dimension " << d << "\n";
    std::cout<<"vector size: "<<size<<", used elements: "<<used_elements<<"\n";
    DenseVector<DT_, IT_> b(size);
    for (Index i (0) ; i < b.size() ; ++i)
    {
      b(i, DT_(i%100) / DT_(100));
    }
    DenseVector<DT_, IT_> x(size, DT_(4711));

    double flops = double(used_elements);
    flops *= 2;

    double bytes = double(used_elements);
    bytes *= DT_(sizeof(DT_));
    bytes += DT_(used_elements * sizeof(IT_));
    bytes += DT_(size * sizeof(DT_));
    Backend::set_preferred_backend(backend);

    auto func = [&] () { sys.apply(b, x); };
    run_bench(func, flops, bytes);

    std::cout<<"control norm: "<<x.norm2()<<"\n";
  }
}

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
#ifdef FEAT_HAVE_CUDA
  run<SparseMatrixCSR<double, unsigned int>, 2, 2, 2 >(PreferredBackend::cuda);
  run<SparseMatrixCSR<double, unsigned long>, 2, 2, 2 >(PreferredBackend::cuda);
  run<SparseMatrixCSR<double, unsigned int>, 3, 3, 3 >(PreferredBackend::cuda);
  run<SparseMatrixCSR<double, unsigned long>, 3, 3, 3 >(PreferredBackend::cuda);
  run<SparseMatrixCSR<double, unsigned int>, 3, 3, 1 >(PreferredBackend::cuda);
  run<SparseMatrixCSR<double, unsigned long>, 3, 3, 1 >(PreferredBackend::cuda);
#endif
  run<SparseMatrixCSR<double, unsigned int>, 2, 2, 2 >(PreferredBackend::generic);
  run<SparseMatrixCSR<double, unsigned long>, 2, 2, 2 >(PreferredBackend::generic);
  run<SparseMatrixCSR<double, unsigned int>, 3, 3, 3 >(PreferredBackend::generic);
  run<SparseMatrixCSR<double, unsigned long>, 3, 3, 3 >(PreferredBackend::generic);
  run<SparseMatrixCSR<double, unsigned int>, 3, 3, 1 >(PreferredBackend::generic);
  run<SparseMatrixCSR<double, unsigned long>, 3, 3, 1 >(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
  run<SparseMatrixCSR<double, unsigned long>, 2, 2, 2 >(PreferredBackend::mkl);
#endif
  return 0;
}

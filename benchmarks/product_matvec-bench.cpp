// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;


template <typename SM_>
void run(PreferredBackend backend)
{
  Backend::set_preferred_backend(PreferredBackend::generic);
  typedef typename SM_::DataType DT_;
  typedef typename SM_::IndexType IT_;

  std::vector<IT_> num_of_nodes;
  num_of_nodes.push_back(4000);
  num_of_nodes.push_back(4000);

  // generate FE matrix A
  SparseMatrixBanded<DT_, IT_> bm(PointstarStructureFE::template value<DT_>(1, num_of_nodes));
  {
    Memory::TypedView<DT_> vb = bm.val_view_w();
    for (Index i(0) ; i < bm.get_elements_size().at(0) ; ++i)
      vb[i] = DT_((i%4) + 1);
  }
  SM_ sys;
  sys.convert(bm);
  Index size(sys.num_rows());

  std::cout<<backend<<" "<<SM_::name()<<" "<<Type::Traits<DT_>::name()<<" "<<Type::Traits<IT_>::name()<<"\n";
  std::cout<<"vector size: "<<size<<" used elements: "<<sys.num_nzes()<<"\n";
  DenseVector<DT_, IT_> x(size);
  {
    Memory::TypedView<DT_> vx = x.elements_view_w();
    for (Index i (0) ; i < x.size() ; ++i)
      vx[i] = DT_(i%100) / DT_(100);
  }
  DenseVector<DT_, IT_> r(size, DT_(4711));

  Backend::set_preferred_backend(backend);

  double flops(double(sys.num_nzes()));
  flops *= 2;

  double bytes(double(sys.num_nzes()));
  bytes *= double(sizeof(DT_));
  bytes += double(sys.num_nzes() * sizeof(IT_));
  bytes += double(size * sizeof(DT_));

  auto func = [&] () { sys.apply(r, x); };
  run_bench(func, flops, bytes);

  Runtime::synchronize_devices();
  std::cout<<"control norm: "<<double(x.norm2())<<"\n";
}

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  run<SparseMatrixCSR<double, Index> >(PreferredBackend::generic);
  run<SparseMatrixCSR<float, Index> >(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
  run<SparseMatrixCSR<double, Index> >(PreferredBackend::mkl);
  run<SparseMatrixCSR<float, Index> >(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_CUDA
  run<SparseMatrixCSR<double, Index> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<double, unsigned int> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<float, Index> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<float, unsigned int> >(PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
  run<SparseMatrixCSR<Half, Index> >(PreferredBackend::cuda);
  run<SparseMatrixCSR<Half, unsigned int> >(PreferredBackend::cuda);
#endif
#endif
  return 0;
}

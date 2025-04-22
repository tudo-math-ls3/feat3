// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/util/type_traits.hpp>
#include <benchmarks/benchmark.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/cuda_util.hpp>

#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Benchmark;


template <typename DT_, typename IT_, bool synch_>
void run([[maybe_unused]] SimpleArgParser& args)
{
  #ifdef FEAT_HAVE_CUDA
  Index lvl = 0;
  args.parse("level", lvl);
  Index num_inner_loop = 100;
  args.parse("n-inner-loop", num_inner_loop);
  Backend::set_preferred_backend(PreferredBackend::generic);
  typedef DT_ DataType;
  typedef IT_ IndexType;

  Index node_size = 10 << lvl;

  std::vector<IT_> num_of_nodes;
  num_of_nodes.push_back(node_size);
  num_of_nodes.push_back(node_size);

  // generate FE matrix A
  SparseMatrixCSR<DT_, IT_> mat(PointstarStructureFE::template value<DT_>(1, num_of_nodes));
  for (Index i(0) ; i < mat.get_elements_size().at(0) ; ++i)
    mat.val()[i] = DT_((i%5) + 1);

  DenseVector<DataType, IndexType> rhs = mat.create_vector_r();
  std::cout << "Run bench " << (synch_ ? String("Synched") : String("Unsynched")) << "\n";
  std::cout<<"vector size: "<< rhs.size() <<" used elements: "<< mat.used_elements()<<"\n";
  for (Index i (0) ; i < rhs.size() ; ++i)
    rhs(i, DT_(i%100) / DT_(100));
  DenseVector<DT_, IT_> sol = mat.create_vector_r();
  sol.format(DT_(4711));

  Backend::set_preferred_backend(PreferredBackend::cuda);

  // implicit transfer to device mem
  mat.apply(sol, rhs);

  double flops(double(mat.used_elements()));
  flops *= 2;
  flops *= num_inner_loop;

  double bytes(double(mat.used_elements()));
  bytes *= double(sizeof(DT_));
  bytes += double(mat.used_elements() * sizeof(IT_));
  bytes += double(rhs.size() * sizeof(DT_));

  bytes *= num_inner_loop;

  if constexpr (synch_)
  {
    auto func = [&] () {
      for(Index i = 0; i < num_inner_loop; ++i)
      {
        mat.apply(sol, rhs);
      }
    };
    run_bench(func, flops, bytes);
  }
  else
  {
    auto func = [&] () {
      Runtime::SyncGuard sync_guard;
      for(Index i = 0; i < num_inner_loop; ++i)
      {
        mat.apply(sol, rhs);
      }
    };
    run_bench(func, flops, bytes);
  }

  MemoryPool::synchronize();
  std::cout<<"control norm: "<<double(sol.norm2())<<"\n";
  #endif
}

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  SimpleArgParser args(argc, argv);
#ifdef FEAT_HAVE_CUDA
  run<double, Index, true>(args);
  run<double, Index, false>(args);
#endif
  return 0;
}
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/product_matvec.hpp>
#include <kernel/lafem/richardson.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class RichardsonTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  RichardsonTest()
    : TaggedTest<Arch_, DT_, Algo_>("richardson_test")
  {
  }

  virtual void run() const
  {
    Index size(1025);
    DenseVector<Arch_, DT_> x(size, DT_(1));
    DenseVector<Arch_, DT_> b(size);
    DenseVector<Mem::Main, DT_> ref_local(size, DT_(42));
    DenseVector<Arch_, DT_> ref(size);
    copy(ref, ref_local);

    SparseMatrixCOO<Mem::Main, DT_> csys(size, size);
    for (Index i(0) ; i < size ; ++i)
      csys(i, i, DT_(4));
    for (Index i(1) ; i < size ; ++i)
      csys(i - 1, i, DT_(-1));
    for (Index i(0) ; i < size - 1; ++i)
      csys(i + 1, i, DT_(-1));
    SparseMatrixCSR<Arch_, DT_> sys(csys);

    JacobiPreconditioner<Algo_, SparseMatrixCSR<Arch_, DT_>, DenseVector<Arch_, DT_> > jac(sys, DT_(0.7));

    ProductMatVec<Algo_>::value(b, sys, ref);

    Richardson<Algo_>::value(x, sys, b, jac, 1000, DT_(1e-16));

    DenseVector<Mem::Main, DT_> sol(size);
    copy(sol, x);
    for (Index i(0) ; i < size ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(sol(i), ref_local(i), 1e-10);
  }
};
RichardsonTest<Mem::Main, Algo::Generic, float> richardson_test_float;
RichardsonTest<Mem::Main, Algo::Generic, double> richardson_test_double;
#ifdef FEAST_GMP
RichardsonTest<Mem::Main, Algo::Generic, mpf_class> richardson_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
RichardsonTest<Mem::CUDA, Algo::CUDA, float> cuda_richardson_test_float;
RichardsonTest<Mem::CUDA, Algo::CUDA, double> cuda_richardson_test_double;
#endif

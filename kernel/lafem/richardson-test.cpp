#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/richardson.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Algo_,
  typename SM_>
class RichardsonTest
  : public TaggedTest<typename SM_::MemType, typename SM_::DataType, Algo_>
{

public:

  typedef typename SM_::DataType DT_;
  typedef typename SM_::MemType Mem_;
  RichardsonTest()
    : TaggedTest<Mem_, DT_, Algo_>("richardson_test " + SM_::type_name())
  {
  }

  virtual void run() const
  {
    //PointstarFactoryFD<DT_> factory(3, 2);
    PointstarFactoryFE<DT_> factory(13);
    SM_ sys(factory.matrix_csr());

    Index size(sys.rows());
    DenseVector<Mem_, DT_> x(size, DT_(1));
    DenseVector<Mem_, DT_> ref(factory.vector_q2_bubble());
    DenseVector<Mem::Main, DT_> ref_local(ref);
    DenseVector<Mem_, DT_> b(size);
    sys.template apply<Algo_>(b, ref);

    JacobiPreconditioner<Algo_, SM_, DenseVector<Mem_, DT_> > jac(sys, DT_(0.7));


    Richardson<Algo_>::value(x, sys, b, jac, 1000, DT_(1e-10));

    DenseVector<Mem::Main, DT_> sol(size);
    sol.copy(x);
    for (Index i(0) ; i < size ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(sol(i), ref_local(i), 1e-9);
  }
};
RichardsonTest<Algo::Generic, SparseMatrixCOO<Mem::Main, double> > coo_richardson_test_double;
RichardsonTest<Algo::Generic, SparseMatrixCSR<Mem::Main, double> > csr_richardson_test_double;
RichardsonTest<Algo::Generic, SparseMatrixELL<Mem::Main, double> > ell_richardson_test_double;
#ifdef FEAST_BACKENDS_MKL
RichardsonTest<Algo::MKL, SparseMatrixCOO<Mem::Main, double> > mkl_coo_richardson_test_double;
RichardsonTest<Algo::MKL, SparseMatrixCSR<Mem::Main, double> > mkl_csr_richardson_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
RichardsonTest<Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double> > cuda_csr_richardson_test_double;
RichardsonTest<Algo::CUDA, SparseMatrixELL<Mem::CUDA, double> > cuda_ell_richardson_test_double;
#endif

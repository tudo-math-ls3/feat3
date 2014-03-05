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
  typename PSF_,
  typename Algo_,
  typename SM_>
class RichardsonTest
  : public TaggedTest<typename SM_::MemType, typename SM_::DataType, Algo_>
{

public:

  typedef typename SM_::DataType DT_;
  typedef typename SM_::MemType Mem_;
  RichardsonTest()
    : TaggedTest<Mem_, DT_, Algo_>("richardson_test " + SM_::name())
  {
  }

  virtual void run() const
  {
    PSF_ factory(13);
    SM_ sys;
    sys.convert(factory.matrix_csr());

    Index size(sys.rows());
    DenseVector<Mem_, DT_> x(size, DT_(1));
    DenseVector<Mem_, DT_> ref;
    ref.convert(factory.vector_q2_bubble());
    DenseVector<Mem_, DT_> b(size);
    sys.template apply<Algo_>(b, ref);

    JacobiPreconditioner<Algo_, SM_, DenseVector<Mem_, DT_> > jac(sys, DT_(0.7));


    Richardson<Algo_>::value(x, sys, b, jac, 5000, DT_(1e-12));

    DenseVector<Mem::Main, DT_> sol(size);
    sol.copy(x);
    for (Index i(0) ; i < size ; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(sol(i), ref(i), 1e-8);
  }
};
RichardsonTest<PointstarFactoryFD<double>, Algo::Generic, SparseMatrixCOO<Mem::Main, double> > coo_fd_richardson_test_double;
RichardsonTest<PointstarFactoryFE<double>, Algo::Generic, SparseMatrixCOO<Mem::Main, double> > coo_fe_richardson_test_double;
RichardsonTest<PointstarFactoryFD<double>, Algo::Generic, SparseMatrixCSR<Mem::Main, double> > csr_fd_richardson_test_double;
RichardsonTest<PointstarFactoryFE<double>, Algo::Generic, SparseMatrixCSR<Mem::Main, double> > csr_fe_richardson_test_double;
RichardsonTest<PointstarFactoryFD<double>, Algo::Generic, SparseMatrixELL<Mem::Main, double> > ell_fd_richardson_test_double;
RichardsonTest<PointstarFactoryFE<double>, Algo::Generic, SparseMatrixELL<Mem::Main, double> > ell_fe_richardson_test_double;
#ifdef FEAST_BACKENDS_MKL
RichardsonTest<PointstarFactoryFD<double>, Algo::MKL, SparseMatrixCOO<Mem::Main, double> > mkl_coo_fd_richardson_test_double;
RichardsonTest<PointstarFactoryFE<double>, Algo::MKL, SparseMatrixCOO<Mem::Main, double> > mkl_coo_fe_richardson_test_double;
RichardsonTest<PointstarFactoryFD<double>, Algo::MKL, SparseMatrixCSR<Mem::Main, double> > mkl_csr_fd_richardson_test_double;
RichardsonTest<PointstarFactoryFE<double>, Algo::MKL, SparseMatrixCSR<Mem::Main, double> > mkl_csr_fe_richardson_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
RichardsonTest<PointstarFactoryFD<double>, Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double> > cuda_csr_fd_richardson_test_double;
RichardsonTest<PointstarFactoryFE<double>, Algo::CUDA, SparseMatrixCSR<Mem::CUDA, double> > cuda_csr_fe_richardson_test_double;
RichardsonTest<PointstarFactoryFD<double>, Algo::CUDA, SparseMatrixELL<Mem::CUDA, double> > cuda_ell_fd_richardson_test_double;
RichardsonTest<PointstarFactoryFE<double>, Algo::CUDA, SparseMatrixELL<Mem::CUDA, double> > cuda_ell_fe_richardson_test_double;
#endif

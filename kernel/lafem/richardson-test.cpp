#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/product.hpp>
#include <kernel/lafem/richardson.hpp>
#include <kernel/lafem/preconditioner.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class RichardsonTest
  : public TaggedTest<Arch_, DT_>
{

public:

  RichardsonTest()
    : TaggedTest<Arch_, DT_>("richardson_test")
  {
  }

  virtual void run() const
  {
    Index size(1025);
    DenseVector<Arch_, DT_> x(size);
    DenseVector<Arch_, DT_> b(size);
    DenseVector<Arch_, DT_> ref(size);
    for (Index i(0) ; i < size ; ++i)
    {
      ref(i, DT_(42));
      x(i, DT_(1));
    }

    SparseMatrixCOO<Arch_, DT_> csys(size, size);
    for (Index i(0) ; i < size ; ++i)
      csys(i, i, DT_(4));
    for (Index i(1) ; i < size ; ++i)
      csys(i - 1, i, DT_(-1));
    for (Index i(0) ; i < size - 1; ++i)
      csys(i + 1, i, DT_(-1));
    SparseMatrixCSR<Arch_, DT_> sys(csys);

    JacobiPreconditioner<BType_, SparseMatrixCSR<Arch_, DT_>, DenseVector<Arch_, DT_> > jac(sys, DT_(0.7));

    Product<Arch_, BType_>::value(b, sys, ref);

    Richardson<BType_>::value(x, sys, b, jac, 1000, DT_(1e-16));

    TEST_CHECK_EQUAL(x, ref);
  }
};
RichardsonTest<Archs::CPU, Archs::Generic, float> richardson_test_float;
RichardsonTest<Archs::CPU, Archs::Generic, double> richardson_test_double;

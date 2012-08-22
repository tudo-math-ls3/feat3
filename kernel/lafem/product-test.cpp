#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename BType_,
  typename DT_>
class CSRProductTest
  : public TaggedTest<Arch_, DT_>
{

public:

  CSRProductTest()
    : TaggedTest<Arch_, DT_>("product_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 3e2 ; size*=2)
    {
      SparseMatrixCOO<Arch_, DT_> ac(size, size + 2);
      DenseVector<Arch_, DT_> b(size + 2);
      DenseVector<Arch_, DT_> c(size, 4711);
      DenseVector<Arch_, DT_> ref(size);
      for (unsigned long row(0) ; row < ac.rows() ; ++row)
      {
        for (unsigned long col(0) ; col < ac.columns() ; ++col)
        {
          if (col % 5 == 0)
            ac(row, col, row * 5 / (col + 1));
        }
      }
      SparseMatrixCSR<Arch_, DT_> a(ac);
      for (Index i(0) ; i < size ; ++i)
      {
        b(i, DT_(size*2 - i/((i+3)/2) * 1.23));
      }
      for (Index i(0) ; i < size ; ++i)
      {
        DT_ sum(0);
        for (Index j(0) ; j < a.columns() ; ++j)
          sum += ac(i, j) * b(j);
        ref(i, sum);
      }

      Product<Arch_, BType_>::value(c, a, b);
      TEST_CHECK_EQUAL(c, ref);
    }
  }
};
CSRProductTest<Archs::CPU, Archs::Generic, float> csrv_product_test_float;
CSRProductTest<Archs::CPU, Archs::Generic, double> csr_product_test_double;

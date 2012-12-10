#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/product_matvec.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename SM_>
class ProductTest
  : public TaggedTest<Arch_, DT_>
{

public:

  ProductTest()
    : TaggedTest<Arch_, DT_>("product_matvec_test: " + SM_::name())
  {
  }

  virtual void run() const
  {
    const DT_ eps = std::pow(std::numeric_limits<DT_>::epsilon(), DT_(0.8));
    const DT_ pi = DT_(2) * std::acos(DT_(0));

    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      DenseVector<Mem::Main, DT_> b_local(size + 2);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (unsigned long row(0) ; row < a_local.rows() ; ++row)
      {
        for (unsigned long col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }
      for (Index i(0) ; i < size ; ++i)
      {
        b_local(i, std::sin(pi * DT_(i) / DT_(size-1)));
      }

      SM_ a(a_local);
      DenseVector<Arch_, DT_> b(size + 2);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size, 4711);

      ProductMatVec<Algo_>::value(c, a, b);
      copy(result_local, c);

      DT_ dev(DT_(0));
      for (Index i(0) ; i < size ; ++i)
      {
        DT_ sum(result_local(i));
        for (Index j(0) ; j < a_local.columns() ; ++j)
          sum -= a_local(i, j) * b_local(j);
        dev = std::max(dev, std::abs(sum));
      }

      TEST_CHECK(dev <= eps);
    }
  }
};
ProductTest<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> >csr_product_matvec_test_float;
ProductTest<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> >csr_product_matvec_test_double;
#ifdef FEAST_BACKENDS_CUDA
ProductTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixCSR<Mem::CUDA, float> >cuda_csr_product_matvec_test_float;
ProductTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixCSR<Mem::CUDA, double> >cuda_csr_product_matvec_test_double;
#endif
ProductTest<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> >ell_product_matvec_test_float;
ProductTest<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> >ell_product_matvec_test_double;
#ifdef FEAST_BACKENDS_CUDA
ProductTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixELL<Mem::CUDA, float> >cuda_ell_product_matvec_test_float;
ProductTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixELL<Mem::CUDA, double> >cuda_ell_product_matvec_test_double;
#endif

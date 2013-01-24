#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/algorithm.hpp>
#include <kernel/lafem/product.hpp>
#include <kernel/lafem/norm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVComponentProductTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVComponentProductTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_component_product_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> ref2(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        b_local(i, DT_(size*2 - i));
        ref(i, a_local(i) * b_local(i));
        ref2(i, a_local(i) * a_local(i));
      }

      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size);

      Product<Algo_>::value(c, a, b);
      copy(result_local, c);
      TEST_CHECK_EQUAL(result_local, ref);

      Product<Algo_>::value(a, a, b);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(a, a_local);
      Product<Algo_>::value(b, a, b);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(b, b_local);
      Product<Algo_>::value(a, a, a);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }
};
DVComponentProductTest<Mem::Main, Algo::Generic, float> dv_component_product_test_float;
DVComponentProductTest<Mem::Main, Algo::Generic, double> dv_component_product_test_double;
#ifdef FEAST_BACKENDS_MKL
DVComponentProductTest<Mem::Main, Algo::MKL, float> mkl_dv_component_product_test_float;
DVComponentProductTest<Mem::Main, Algo::MKL, double> mkl_dv_component_product_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVComponentProductTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_component_product_test_float;
DVComponentProductTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_component_product_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename SM_>
class ProductMatVecTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  ProductMatVecTest()
    : TaggedTest<Arch_, DT_, Algo_>("product_matvec_test: " + SM_::type_name())
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

      Product<Algo_>::value(c, a, b);
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
ProductMatVecTest<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> > csr_product_matvec_test_float;
ProductMatVecTest<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> > csr_product_matvec_test_double;
#ifdef FEAST_BACKENDS_MKL
ProductMatVecTest<Mem::Main, Algo::MKL, float, SparseMatrixCSR<Mem::Main, float> > mkl_csr_product_matvec_test_float;
ProductMatVecTest<Mem::Main, Algo::MKL, double, SparseMatrixCSR<Mem::Main, double> > mkl_csr_product_matvec_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
ProductMatVecTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixCSR<Mem::CUDA, float> > cuda_csr_product_matvec_test_float;
ProductMatVecTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixCSR<Mem::CUDA, double> > cuda_csr_product_matvec_test_double;
#endif
ProductMatVecTest<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> > ell_product_matvec_test_float;
ProductMatVecTest<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> > ell_product_matvec_test_double;
#ifdef FEAST_BACKENDS_CUDA
ProductMatVecTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixELL<Mem::CUDA, float> > cuda_ell_product_matvec_test_float;
ProductMatVecTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixELL<Mem::CUDA, double> > cuda_ell_product_matvec_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVDotProductTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVDotProductTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_dot_product_test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = std::pow(std::numeric_limits<DT_>::epsilon(), DT_(0.8));

    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      const DT_ den(DT_(1) / DT_(size));
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i+1) * den);    // a[i] = (i+1) / n
        b_local(i, DT_(1) / DT_(i+1)); // b[i] = 1 / (i+1)
      }

      DenseVector<Arch_, DT_> a(a_local);
      DenseVector<Arch_, DT_> b(b_local);

      // a*b = 1
      DT_ ref(DT_(1));
      DT_ c = Product<Algo_>::value(a, b);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = Product<Algo_>::value(b, a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = Product<Algo_>::value(b, b);
      ref = Norm2<Algo_>::value(b);
      ref *= ref;
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
    }
  }
};
DVDotProductTest<Mem::Main, Algo::Generic, float> dv_dot_product_test_float;
DVDotProductTest<Mem::Main, Algo::Generic, double> dv_dot_product_test_double;
#ifdef FEAST_BACKENDS_MKL
DVDotProductTest<Mem::Main, Algo::MKL, float> mkl_dv_dot_product_test_float;
DVDotProductTest<Mem::Main, Algo::MKL, double> mkl_dv_dot_product_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVDotProductTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_dot_product_test_float;
DVDotProductTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_dot_product_test_double;
#endif

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVComponentProductTest2
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVComponentProductTest2()
    : TaggedTest<Arch_, DT_, Algo_>("product: dv_component_product_test")
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

      c.template product<Algo_>(a, b);
      copy(result_local, c);
      TEST_CHECK_EQUAL(result_local, ref);

      a.template product<Algo_>(a, b);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(a, a_local);
      b.template product<Algo_>(a, b);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(b, b_local);
      a.template product<Algo_>(a, a);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }
};
DVComponentProductTest2<Mem::Main, Algo::Generic, float> dv_component_product_test_float2;
DVComponentProductTest2<Mem::Main, Algo::Generic, double> dv_component_product_test_double2;
#ifdef FEAST_GMP
//DVComponentProductTest2<Mem::Main, Algo::Generic, mpf_class> dv_component_product_test_mpf_class2;
#endif
#ifdef FEAST_BACKENDS_MKL
DVComponentProductTest2<Mem::Main, Algo::MKL, float> mkl_dv_component_product_test_float2;
DVComponentProductTest2<Mem::Main, Algo::MKL, double> mkl_dv_component_product_test_double2;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVComponentProductTest2<Mem::CUDA, Algo::CUDA, float> cuda_dv_component_product_test_float2;
DVComponentProductTest2<Mem::CUDA, Algo::CUDA, double> cuda_dv_component_product_test_double2;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename SM_>
class ProductMatVecTest2
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  ProductMatVecTest2()
    : TaggedTest<Arch_, DT_, Algo_>("product: product_matvec_test: " + SM_::type_name())
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::Limits<DT_>::epsilon(), DT_(0.8));
    const DT_ pi = Math::Limits<DT_>::pi();

    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      DenseVector<Mem::Main, DT_> b_local(size + 2);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }
      for (Index i(0) ; i < size ; ++i)
      {
        b_local(i, Math::sin(pi * DT_(i) / DT_(size-1)));
      }

      SM_ a(a_local);
      DenseVector<Arch_, DT_> b(size + 2);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size, 4711);

      c.template product<Algo_>(a, b);
      copy(result_local, c);

      DT_ dev(DT_(0));
      for (Index i(0) ; i < size ; ++i)
      {
        DT_ sum(result_local(i));
        for (Index j(0) ; j < a_local.columns() ; ++j)
          sum -= a_local(i, j) * b_local(j);
        dev = Math::max(dev, Math::abs(sum));
      }

      TEST_CHECK(dev <= eps);
    }
  }
};
ProductMatVecTest2<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> > csr_product_matvec_test_float2;
ProductMatVecTest2<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> > csr_product_matvec_test_double2;
#ifdef FEAST_GMP
//ProductMatVecTest2<Mem::Main, Algo::Generic, mpf_class, SparseMatrixCSR<Mem::Main, mpf_class> > csr_product_matvec_test_mpf_class2;
#endif
#ifdef FEAST_BACKENDS_MKL
ProductMatVecTest2<Mem::Main, Algo::MKL, float, SparseMatrixCSR<Mem::Main, float> > mkl_csr_product_matvec_test_float2;
ProductMatVecTest2<Mem::Main, Algo::MKL, double, SparseMatrixCSR<Mem::Main, double> > mkl_csr_product_matvec_test_double2;
#endif
#ifdef FEAST_BACKENDS_CUDA
ProductMatVecTest2<Mem::CUDA, Algo::CUDA, float, SparseMatrixCSR<Mem::CUDA, float> > cuda_csr_product_matvec_test_float2;
ProductMatVecTest2<Mem::CUDA, Algo::CUDA, double, SparseMatrixCSR<Mem::CUDA, double> > cuda_csr_product_matvec_test_double2;
#endif
ProductMatVecTest2<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> > ell_product_matvec_test_float2;
ProductMatVecTest2<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> > ell_product_matvec_test_double2;
#ifdef FEAST_BACKENDS_CUDA
ProductMatVecTest2<Mem::CUDA, Algo::CUDA, float, SparseMatrixELL<Mem::CUDA, float> > cuda_ell_product_matvec_test_float2;
ProductMatVecTest2<Mem::CUDA, Algo::CUDA, double, SparseMatrixELL<Mem::CUDA, double> > cuda_ell_product_matvec_test_double2;
#endif
ProductMatVecTest2<Mem::Main, Algo::Generic, float, SparseMatrixCOO<Mem::Main, float> > coo_product_matvec_test_float2;
ProductMatVecTest2<Mem::Main, Algo::Generic, double, SparseMatrixCOO<Mem::Main, double> > coo_product_matvec_test_double2;
#ifdef FEAST_BACKENDS_MKL
ProductMatVecTest2<Mem::Main, Algo::MKL, float, SparseMatrixCOO<Mem::Main, float> > mkl_coo_product_matvec_test_float2;
ProductMatVecTest2<Mem::Main, Algo::MKL, double, SparseMatrixCOO<Mem::Main, double> > mkl_coo_product_matvec_test_double2;
#endif
#ifdef FEAST_BACKENDS_CUDA
//ProductMatVecTest2<Mem::CUDA, Algo::CUDA, float, SparseMatrixCOO<Mem::CUDA, float> > cuda_coo_product_matvec_test_float2;
//ProductMatVecTest2<Mem::CUDA, Algo::CUDA, double, SparseMatrixCOO<Mem::CUDA, double> > cuda_coo_product_matvec_test_double2;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVDotProductTest2
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVDotProductTest2()
    : TaggedTest<Arch_, DT_, Algo_>("product: dv_dot_product_test")
  {
  }

  virtual void run() const
  {
    const DT_ eps = Math::pow(Math::Limits<DT_>::epsilon(), DT_(0.8));

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
      DT_ c = a.template product<Algo_>(b);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.template product<Algo_>(a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.template product<Algo_>(b);
      ref = b.template norm2<Algo_>();
      ref *= ref;
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
    }
  }
};
DVDotProductTest2<Mem::Main, Algo::Generic, float> dv_dot_product_test_float2;
DVDotProductTest2<Mem::Main, Algo::Generic, double> dv_dot_product_test_double2;
#ifdef FEAST_GMP
//DVDotProductTest2<Mem::Main, Algo::Generic, mpf_class> dv_dot_product_test_mpf_class2;
#endif
#ifdef FEAST_BACKENDS_MKL
DVDotProductTest2<Mem::Main, Algo::MKL, float> mkl_dv_dot_product_test_float2;
DVDotProductTest2<Mem::Main, Algo::MKL, double> mkl_dv_dot_product_test_double2;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVDotProductTest2<Mem::CUDA, Algo::CUDA, float> cuda_dv_dot_product_test_float2;
DVDotProductTest2<Mem::CUDA, Algo::CUDA, double> cuda_dv_dot_product_test_double2;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVScaleTest2
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVScaleTest2()
    : TaggedTest<Arch_, DT_, Algo_>("product: dv_scale_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DT_ s(DT_(4.321));
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        ref(i, a_local(i) * s);
      }

      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);

      b.template product<Algo_>(a, s);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      a.template product<Algo_>(a, s);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);
    }
  }
};
DVScaleTest2<Mem::Main, Algo::Generic, float> dv_scale_test_float;
DVScaleTest2<Mem::Main, Algo::Generic, double> dv_scale_test_double;
#ifdef FEAST_GMP
//DVScaleTest2<Mem::Main, Algo::Generic, mpf_class> dv_scale_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_MKL
DVScaleTest2<Mem::Main, Algo::MKL, float> mkl_dv_scale_test_float;
DVScaleTest2<Mem::Main, Algo::MKL, double> mkl_dv_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVScaleTest2<Mem::CUDA, Algo::CUDA, float> cuda_dv_scale_test_float;
DVScaleTest2<Mem::CUDA, Algo::CUDA, double> cuda_dv_scale_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename SM_>
class SMScaleTest2
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  SMScaleTest2()
    : TaggedTest<Arch_, DT_, Algo_>("product: scale_test: " + SM_::type_name())
  {
  }

  virtual void run() const
  {
    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      SparseMatrixCOO<Mem::Main, DT_> ref_local(size, size + 2);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));

          if(row == col)
            ref_local(row, col, DT_(2) * s);
          else if((row == col+1) || (row+1 == col))
            ref_local(row, col, DT_(-1) * s);
        }
      }

      SM_ a(a_local);
      SM_ b(a.clone());

      b.template product<Algo_>(a, s);
      SparseMatrixCOO<Mem::Main, DT_> b_local(b);
      TEST_CHECK_EQUAL(b_local, ref_local);

      a.template product<Algo_>(a, s);
      SparseMatrixCOO<Arch_, DT_> a_coo(a);
      a_local = a_coo;
      TEST_CHECK_EQUAL(a_local, ref_local);
    }
  }
};
SMScaleTest2<Mem::Main, Algo::Generic, float, SparseMatrixCOO<Mem::Main, float> > sm_coo_scale_test_float;
SMScaleTest2<Mem::Main, Algo::Generic, double, SparseMatrixCOO<Mem::Main, double> > sm_coo_scale_test_double;
#ifdef FEAST_GMP
//SMScaleTest2<Mem::Main, Algo::Generic, mpf_class, SparseMatrixCOO<Mem::Main, mpf_class> > sm_coo_scale_test_mpf_class;
#endif
SMScaleTest2<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> > sm_csr_scale_test_float;
SMScaleTest2<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> > sm_csr_scale_test_double;
#ifdef FEAST_GMP
//SMScaleTest2<Mem::Main, Algo::Generic, mpf_class, SparseMatrixCSR<Mem::Main, mpf_class> > sm_csr_scale_test_mpf_class;
#endif
SMScaleTest2<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> > sm_ell_scale_test_float;
SMScaleTest2<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> > sm_ell_scale_test_double;
#ifdef FEAST_GMP
//SMScaleTest2<Mem::Main, Algo::Generic, mpf_class, SparseMatrixELL<Mem::Main, mpf_class> > sm_ell_scale_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_MKL
SMScaleTest2<Mem::Main, Algo::MKL, float, SparseMatrixCOO<Mem::Main, float> > mkl_sm_coo_scale_test_float;
SMScaleTest2<Mem::Main, Algo::MKL, double, SparseMatrixCOO<Mem::Main, double> > mkl_sm_coo_scale_test_double;
SMScaleTest2<Mem::Main, Algo::MKL, float, SparseMatrixCSR<Mem::Main, float> > mkl_sm_csr_scale_test_float;
SMScaleTest2<Mem::Main, Algo::MKL, double, SparseMatrixCSR<Mem::Main, double> > mkl_sm_csr_scale_test_double;
SMScaleTest2<Mem::Main, Algo::MKL, float, SparseMatrixELL<Mem::Main, float> > mkl_sm_ell_scale_test_float;
SMScaleTest2<Mem::Main, Algo::MKL, double, SparseMatrixELL<Mem::Main, double> > mkl_sm_ell_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
SMScaleTest2<Mem::CUDA, Algo::CUDA, float, SparseMatrixCOO<Mem::CUDA, float> > cuda_sm_coo_scale_test_float;
SMScaleTest2<Mem::CUDA, Algo::CUDA, double, SparseMatrixCOO<Mem::CUDA, double> > cuda_sm_coo_scale_test_double;
SMScaleTest2<Mem::CUDA, Algo::CUDA, float, SparseMatrixCSR<Mem::CUDA, float> > cuda_sm_csr_scale_test_float;
SMScaleTest2<Mem::CUDA, Algo::CUDA, double, SparseMatrixCSR<Mem::CUDA, double> > cuda_sm_csr_scale_test_double;
SMScaleTest2<Mem::CUDA, Algo::CUDA, float, SparseMatrixELL<Mem::CUDA, float> > cuda_sm_ell_scale_test_float;
SMScaleTest2<Mem::CUDA, Algo::CUDA, double, SparseMatrixELL<Mem::CUDA, double> > cuda_sm_ell_scale_test_double;
#endif

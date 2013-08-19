#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sum.hpp>
#include <kernel/lafem/scale.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVSumTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVSumTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_sum_test")
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
        ref(i, a_local(i) + b_local(i));
        ref2(i, a_local(i) + a_local(i));
      }

      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size);

      Sum<Algo_>::value(c, a, b);
      copy(result_local, c);
      TEST_CHECK_EQUAL(result_local, ref);

      Sum<Algo_>::value(a, a, b);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(a, a_local);
      Sum<Algo_>::value(b, a, b);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(b, b_local);
      Sum<Algo_>::value(a, a, a);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }
};
DVSumTest<Mem::Main, Algo::Generic, float> dv_sum_test_float;
DVSumTest<Mem::Main, Algo::Generic, double> dv_sum_test_double;
#ifdef FEAST_GMP
DVSumTest<Mem::Main, Algo::Generic, mpf_class> dv_sum_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_MKL
DVSumTest<Mem::Main, Algo::MKL, float> mkl_dv_sum_test_float;
DVSumTest<Mem::Main, Algo::MKL, double> mkl_dv_sum_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVSumTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_sum_test_float;
DVSumTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_sum_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class SMCSRSumTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  SMCSRSumTest()
    : TaggedTest<Arch_, DT_, Algo_>("smcsr_sum_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 3e2 ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
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
      SparseMatrixCSR<Arch_, DT_> a(a_local);
      SparseMatrixCSR<Arch_, DT_> b(a.clone());
      DenseVector<Arch_, DT_> b_val(b.used_elements(), b.val());
      Scale<Algo_>::value(b_val, b_val, DT_(2));

      SparseMatrixCSR<Arch_, DT_> r(a.clone());

      Sum<Algo_>::value(r, a, b);

      DenseVector<Arch_, DT_> a_val(a.used_elements(), a.val());
      DenseVector<Arch_, DT_> r_val(r.used_elements(), r.val());

      for (Index i(0) ; i < r_val.size() ; ++i)
      {
        TEST_CHECK_EQUAL(r_val(i), a_val(i) + b_val(i));
      }
    }
  }
};
SMCSRSumTest<Mem::Main, Algo::Generic, float> smcsr_sum_test_float;
SMCSRSumTest<Mem::Main, Algo::Generic, double> smcsr_sum_test_double;
#ifdef FEAST_GMP
SMCSRSumTest<Mem::Main, Algo::Generic, mpf_class> smcsr_sum_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
SMCSRSumTest<Mem::CUDA, Algo::CUDA, float> cuda_smcsr_sum_test_float;
SMCSRSumTest<Mem::CUDA, Algo::CUDA, double> cuda_smcsr_sum_test_double;
#endif

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/scale.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVScaleTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVScaleTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_scale_test")
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

      Scale<Algo_>::value(b, a, s);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      Scale<Algo_>::value(a, a, s);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);
    }
  }
};
DVScaleTest<Mem::Main, Algo::Generic, float> dv_scale_test_float;
DVScaleTest<Mem::Main, Algo::Generic, double> dv_scale_test_double;
#ifdef FEAST_BACKENDS_MKL
DVScaleTest<Mem::Main, Algo::MKL, float> mkl_dv_scale_test_float;
DVScaleTest<Mem::Main, Algo::MKL, double> mkl_dv_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVScaleTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_scale_test_float;
DVScaleTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_scale_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename SM_>
class SMScaleTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  SMScaleTest()
    : TaggedTest<Arch_, DT_, Algo_>("scale_test: " + SM_::type_name())
  {
  }

  virtual void run() const
  {
    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      SparseMatrixCOO<Mem::Main, DT_> ref_local(size, size + 2);
      for (unsigned long row(0) ; row < a_local.rows() ; ++row)
      {
        for (unsigned long col(0) ; col < a_local.columns() ; ++col)
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

      Scale<Algo_>::value(b, a, s);
      SparseMatrixCOO<Mem::Main, DT_> b_local(b);
      TEST_CHECK_EQUAL(b_local, ref_local);

      Scale<Algo_>::value(a, a, s);
      SparseMatrixCOO<Arch_, DT_> a_coo(a);
      a_local = a_coo;
      TEST_CHECK_EQUAL(a_local, ref_local);
    }
  }
};
SMScaleTest<Mem::Main, Algo::Generic, float, SparseMatrixCOO<Mem::Main, float> > sm_coo_scale_test_float;
SMScaleTest<Mem::Main, Algo::Generic, double, SparseMatrixCOO<Mem::Main, double> > sm_coo_scale_test_double;
SMScaleTest<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> > sm_csr_scale_test_float;
SMScaleTest<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> > sm_csr_scale_test_double;
SMScaleTest<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> > sm_ell_scale_test_float;
SMScaleTest<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> > sm_ell_scale_test_double;
#ifdef FEAST_BACKENDS_MKL
SMScaleTest<Mem::Main, Algo::MKL, float, SparseMatrixCOO<Mem::Main, float> > mkl_sm_coo_scale_test_float;
SMScaleTest<Mem::Main, Algo::MKL, double, SparseMatrixCOO<Mem::Main, double> > mkl_sm_coo_scale_test_double;
SMScaleTest<Mem::Main, Algo::MKL, float, SparseMatrixCSR<Mem::Main, float> > mkl_sm_csr_scale_test_float;
SMScaleTest<Mem::Main, Algo::MKL, double, SparseMatrixCSR<Mem::Main, double> > mkl_sm_csr_scale_test_double;
SMScaleTest<Mem::Main, Algo::MKL, float, SparseMatrixELL<Mem::Main, float> > mkl_sm_ell_scale_test_float;
SMScaleTest<Mem::Main, Algo::MKL, double, SparseMatrixELL<Mem::Main, double> > mkl_sm_ell_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
SMScaleTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixCOO<Mem::CUDA, float> > cuda_sm_coo_scale_test_float;
SMScaleTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixCOO<Mem::CUDA, double> > cuda_sm_coo_scale_test_double;
SMScaleTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixCSR<Mem::CUDA, float> > cuda_sm_csr_scale_test_float;
SMScaleTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixCSR<Mem::CUDA, double> > cuda_sm_csr_scale_test_double;
SMScaleTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixELL<Mem::CUDA, float> > cuda_sm_ell_scale_test_float;
SMScaleTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixELL<Mem::CUDA, double> > cuda_sm_ell_scale_test_double;
#endif

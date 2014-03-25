#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_band.hpp>
#include <kernel/util/binary_stream.hpp>

#include <kernel/util/random.hpp>
#include <sstream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the sparse matrix band class.
*
* \test test description missing
*
* \tparam Mem_
* description missing
*
* \tparam DT_
* description missing
*
* \author Christoph Lohmnn
*/
template<
  typename Mem_,
  typename DT_>
class SparseMatrixBandTest
  : public TaggedTest<Mem_, DT_>
{
public:
  SparseMatrixBandTest()
    : TaggedTest<Mem_, DT_>("SparseMatrixBandTest")
  {
  }

  typedef SparseMatrixBand<Mem_, DT_, FiniteElementType::fe_q1, Index> BM_;
  typedef Algo::Generic Algo_;

  virtual void run() const
  {
    Random rand;
    const Index npr(10);
    const Index npc(15);

    DenseVector<Mem_, DT_, Index> val(9 * npr * npc);
    for (Index i(0); i < val.size(); ++i)
    {
      val(i, rand(DT_(0), DT_(1)));
    }

    BM_ sys(npr, npc, val);

    auto x(sys.create_vector_r());
    auto y1(sys.create_vector_l());
    auto y2(sys.create_vector_l());

    for (Index i(0); i < x.size(); ++i)
    {
      x(i, rand(DT_(-1), DT_(1)));
      y2(i, DT_(0));
    }

    sys.template apply<Algo_>(y1, x);

    for(Index i(0); i < sys.rows(); ++i)
    {
      for(Index j(0); j < sys.columns(); ++j)
      {
        y2(i, y2(i) + sys(i, j) * x(j));
      }
    }

    // check, if the result is correct
    for (Index i(0) ; i < y1.size() ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(y1(i), y2(i), 1e-8);
    }

    sys.template apply<Algo_>(y2, x, y1, DT_(-1.0));

    // check, if the result is correct
    for (Index i(0) ; i < y1.size() ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(DT_(0.0), y2(i), 1e-8);
    }

    DenseVector<Mem_, DT_, Index> val2(9 * npr * npc);
    for (Index i(0); i < val2.size(); ++i)
    {
      val2(i, val(i));
    }
    BM_ sys2(npr, npc, val2);
    sys2.template scale<Algo_>(sys, DT_(2.0));
    sys2.template apply<Algo_>(y1, x);

    sys.template apply<Algo_>(y2, x, y1, DT_(-2.0));

    // check, if the result is correct
    for (Index i(0) ; i < y1.size() ; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(DT_(0.0), y2(i), 1e-8);
    }
}
};
SparseMatrixBandTest<Mem::Main, float> cpu_sparse_matrix_band_test_float;
SparseMatrixBandTest<Mem::Main, double> cpu_sparse_matrix_band_test_double;
// #ifdef FEAST_BACKENDS_CUDA
// SparseMatrixBandTest<Mem::CUDA, float> cuda_sparse_matrix_band_test_float;
// SparseMatrixBandTest<Mem::CUDA, double> cuda_sparse_matrix_band_test_double;
// #endif


// template<
//   typename Mem_,
//   typename Algo_,
//   typename DT_>
// class SparseMatrixBandApplyTest
//   : public TaggedTest<Mem_, DT_, Algo_>
// {
// public:
//   SparseMatrixBandApplyTest()
//     : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixBandApplyTest")
//   {
//   }

//   virtual void run() const
//   {
//   }
// };

// SparseMatrixBandApplyTest<Mem::Main, Algo::Generic, float> sm_band_apply_test_float;
// SparseMatrixBandApplyTest<Mem::Main, Algo::Generic, double> sm_band_apply_test_double;
// // #ifdef HONEI_BACKENDS_MKL
// // SparseMatrixBandApplyTest<Mem::Main, Algo::MKL, float> mkl_sm_band_apply_test_float;
// // SparseMatrixBandApplyTest<Mem::Main, Algo::MKL, double> mkl_sm_band_apply_test_double;
// // #endif
// // #ifdef FEAST_BACKENDS_CUDA
// // SparseMatrixBandApplyTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_band_apply_test_float;
// // SparseMatrixBandApplyTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_band_apply_test_double;
// // #endif

// template<
//   typename Mem_,
//   typename Algo_,
//   typename DT_>
// class SparseMatrixBandScaleTest
//   : public TaggedTest<Mem_, DT_, Algo_>
// {
// public:
//   SparseMatrixBandScaleTest()
//     : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixBandScaleTest")
//   {
//   }

//   virtual void run() const
//   {
//   }
// };
// SparseMatrixBandScaleTest<Mem::Main, Algo::Generic, float> sm_band_scale_test_float;
// SparseMatrixBandScaleTest<Mem::Main, Algo::Generic, double> sm_band_scale_test_double;
// // #ifdef FEAST_BACKENDS_MKL
// // SparseMatrixBandScaleTest<Mem::Main, Algo::MKL, float> mkl_sm_band_scale_test_float;
// // SparseMatrixBandScaleTest<Mem::Main, Algo::MKL, double> mkl_sm_band_scale_test_double;
// // #endif
// // #ifdef FEAST_BACKENDS_CUDA
// // SparseMatrixBandScaleTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_band_scale_test_float;
// // SparseMatrixBandScaleTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_band_scale_test_double;
// // #endif

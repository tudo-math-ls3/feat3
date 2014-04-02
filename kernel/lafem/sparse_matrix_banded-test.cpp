#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
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
* \author Christoph Lohmann
*/
template<
  typename Mem_,
  typename DT_>
class SparseMatrixBandedTest
  : public TaggedTest<Mem_, DT_>
{
public:
  SparseMatrixBandedTest()
    : TaggedTest<Mem_, DT_>("SparseMatrixBandedTest")
  {
  }

  typedef SparseMatrixBanded<Mem_, DT_, Index> BM_;
  typedef Algo::Generic Algo_;

  virtual void run() const
  {
    Random random;
    const Index size(9);

    DenseVector<Mem_, Index, Index> offsets(4);
    DenseVector<Mem_, DT_, Index> val(offsets.size() * size);

    offsets(0, 3);
    offsets(1, 4);
    offsets(2, 8);
    offsets(3, 13);

    for (Index i(0); i < val.size(); ++i)
    {
      val(i, random(DT_(0), DT_(10)));
    }

    BM_ sys(size, size - 1, val, offsets);

    auto x(sys.create_vector_r());
    auto y1(sys.create_vector_l());
    auto y2(sys.create_vector_l());

    for (Index i(0); i < x.size(); ++i)
    {
      x(i, random(DT_(-1), DT_(1)));
    }

    for (Index i(0); i < y2.size(); ++i)
    {
      y2(i, 0);
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

    // sys.template apply<Algo_>(y2, x, y1, DT_(-1.0));

    // // check, if the result is correct
    // for (Index i(0) ; i < y1.size() ; ++i)
    // {
    //   TEST_CHECK_EQUAL_WITHIN_EPS(DT_(0.0), y2(i), 1e-8);
    // }

    // DenseVector<Mem_, DT_, Index> val2(9 * npr * npc);
    // for (Index i(0); i < val2.size(); ++i)
    // {
    //   val2(i, val(i));
    // }
    // BM_ sys2(npr, npc, val2, offsets);
    // sys2.template scale<Algo_>(sys, DT_(2.0));
    // sys2.template apply<Algo_>(y1, x);

    // sys.template apply<Algo_>(y2, x, y1, DT_(-2.0));

    // // check, if the result is correct
    // for (Index i(0) ; i < y1.size() ; ++i)
    // {
    //   TEST_CHECK_EQUAL_WITHIN_EPS(DT_(0.0), y2(i), 1e-8);
    // }
  }
};
SparseMatrixBandedTest<Mem::Main, float> cpu_sparse_matrix_banded_test_float;
SparseMatrixBandedTest<Mem::Main, double> cpu_sparse_matrix_banded_test_double;
// #ifdef FEAST_BACKENDS_CUDA
// SparseMatrixBandedTest<Mem::CUDA, float> cuda_sparse_matrix_banded_test_float;
// SparseMatrixBandedTest<Mem::CUDA, double> cuda_sparse_matrix_banded_test_double;
// #endif


// template<
//   typename Mem_,
//   typename Algo_,
//   typename DT_>
// class SparseMatrixBandedApplyTest
//   : public TaggedTest<Mem_, DT_, Algo_>
// {
// public:
//   SparseMatrixBandedApplyTest()
//     : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixBandedApplyTest")
//   {
//   }

//   virtual void run() const
//   {
//   }
// };

// SparseMatrixBandedApplyTest<Mem::Main, Algo::Generic, float> sm_banded_apply_test_float;
// SparseMatrixBandedApplyTest<Mem::Main, Algo::Generic, double> sm_banded_apply_test_double;
// // #ifdef HONEI_BACKENDS_MKL
// // SparseMatrixBandedApplyTest<Mem::Main, Algo::MKL, float> mkl_sm_banded_apply_test_float;
// // SparseMatrixBandedApplyTest<Mem::Main, Algo::MKL, double> mkl_sm_banded_apply_test_double;
// // #endif
// // #ifdef FEAST_BACKENDS_CUDA
// // SparseMatrixBandedApplyTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_banded_apply_test_float;
// // SparseMatrixBandedApplyTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_banded_apply_test_double;
// // #endif

// template<
//   typename Mem_,
//   typename Algo_,
//   typename DT_>
// class SparseMatrixBandedScaleTest
//   : public TaggedTest<Mem_, DT_, Algo_>
// {
// public:
//   SparseMatrixBandedScaleTest()
//     : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixBandedScaleTest")
//   {
//   }

//   virtual void run() const
//   {
//   }
// };
// SparseMatrixBandedScaleTest<Mem::Main, Algo::Generic, float> sm_banded_scale_test_float;
// SparseMatrixBandedScaleTest<Mem::Main, Algo::Generic, double> sm_banded_scale_test_double;
// // #ifdef FEAST_BACKENDS_MKL
// // SparseMatrixBandedScaleTest<Mem::Main, Algo::MKL, float> mkl_sm_banded_scale_test_float;
// // SparseMatrixBandedScaleTest<Mem::Main, Algo::MKL, double> mkl_sm_banded_scale_test_double;
// // #endif
// // #ifdef FEAST_BACKENDS_CUDA
// // SparseMatrixBandedScaleTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_banded_scale_test_float;
// // SparseMatrixBandedScaleTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_banded_scale_test_double;
// // #endif

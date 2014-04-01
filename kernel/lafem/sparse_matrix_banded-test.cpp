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
* \author Christoph Lohmnn
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
    Random rand;
    const Index npr(5);
    const Index npc(4);

    DenseVector<Mem_, DT_, Index> val(9 * npr * npc);
    DenseVector<Mem_, Index> offsets(9);

    offsets(0, npr * npc - npr - 2);
    offsets(1, npr * npc - npr - 1);
    offsets(2, npr * npc - npr);
    offsets(3, npr * npc - 2);
    offsets(4, npr * npc - 1);
    offsets(5, npr * npc);
    offsets(6, npr * npc + npr - 2);
    offsets(7, npr * npc + npr - 1);
    offsets(8, npr * npc + npr);

    for (Index i(0); i < val.size(); ++i)
    {
      val(i, DT_(7)); // rand(DT_(0), DT_(1)));
    }

    BM_ sys(npr * npc, npr * npc - 7, val, offsets);

    std::cout << sys << std::endl;

    // auto x(sys.create_vector_r());
    // auto y1(sys.create_vector_l());
    // auto y2(sys.create_vector_l());

    // for (Index i(0); i < x.size(); ++i)
    // {
    //   x(i, rand(DT_(-1), DT_(1)));
    //   y2(i, DT_(0));
    // }

    // sys.template apply<Algo_>(y1, x);

    // for(Index i(0); i < sys.rows(); ++i)
    // {
    //   for(Index j(0); j < sys.columns(); ++j)
    //   {
    //     y2(i, y2(i) + sys(i, j) * x(j));
    //   }
    // }

    // // check, if the result is correct
    // for (Index i(0) ; i < y1.size() ; ++i)
    // {
    //   TEST_CHECK_EQUAL_WITHIN_EPS(y1(i), y2(i), 1e-8);
    // }

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

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/power_mirror.hpp>
#include <kernel/lafem/vector_mirror.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template<
  typename DataType_,
  typename IndexType_>
class MetaMirrorTest :
  public UnitTest
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;

  typedef DenseVector<DataType, IndexType> BufferVector;

  typedef DenseVector<DataType, IndexType> ScalarVector;
  typedef PowerVector<ScalarVector, 2> PowerVector2;
  typedef TupleVector<PowerVector2, ScalarVector> MetaVector;

  typedef VectorMirror<DataType, IndexType> ScalarMirror;
  typedef PowerMirror<ScalarMirror, 2> PowerMirror2;
  typedef TupleMirror<PowerMirror2, ScalarMirror> MetaMirror;

  explicit MetaMirrorTest(PreferredBackend backend) :
    UnitTest("MetaMirrorTest", Type::Traits<DataType>::name(), Type::Traits<IndexType>::name(), backend)
  {
  }

  static ScalarMirror gen_mirror_x(IndexType m)
  {
    ScalarMirror mir(m*m, m);
    Memory::TypedView<IndexType> ci(mir.indices_view_w());
    IndexType m2(m*(m-1));
    for(IndexType i(0); i < m; ++i)
      ci[i] = m2+i;
    return mir;
  }

  static ScalarMirror gen_mirror_y(IndexType m)
  {
    ScalarMirror mir(m*m, m);
    Memory::TypedView<IndexType> ci(mir.indices_view_w());
    for(IndexType i(0); i < m; ++i)
      ci[i] = i;
    return mir;
  }

  virtual void run() const override
  {
    const DataType tol = TestSystem::tol<DataType>();

    const Index m = 3;
    const Index n = m*m;

    // create the meta-mirrors
    MetaMirror mirror_x(PowerMirror2(gen_mirror_x(m)), gen_mirror_x(m));
    MetaMirror mirror_y(PowerMirror2(gen_mirror_y(m)), gen_mirror_y(m));

    // create meta-vectors
    MetaVector vec_x;
    vec_x.template at<0>().template at<0>() = ScalarVector(n, DataType(1));
    vec_x.template at<0>().template at<1>() = ScalarVector(n, DataType(2));
    vec_x.template at<1>() = ScalarVector(n, DataType(3));
    MetaVector vec_y;
    vec_y.template at<0>().template at<0>() = ScalarVector(n, DataType(-1));
    vec_y.template at<0>().template at<1>() = ScalarVector(n, DataType( 1));
    vec_y.template at<1>() = ScalarVector(n, DataType(-2));

    // create reference synced vectors
    MetaVector sync_x(vec_x.clone());
    MetaVector sync_y(vec_y.clone());
    {
      Memory::TypedView<DataType> vx00 = sync_x.template at<0>().template at<0>().elements_view_rw();
      Memory::TypedView<DataType> vy00 = sync_y.template at<0>().template at<0>().elements_view_rw();
      Memory::TypedView<DataType> vx01 = sync_x.template at<0>().template at<1>().elements_view_rw();
      Memory::TypedView<DataType> vy01 = sync_y.template at<0>().template at<1>().elements_view_rw();
      Memory::TypedView<DataType> vx1 = sync_x.template at<1>().elements_view_rw();
      Memory::TypedView<DataType> vy1 = sync_y.template at<1>().elements_view_rw();
      for(Index i(0); i < m; ++i)
      {
        Index k(n-m+i);
        vx00[k] = DataType_(0);
        vy00[i] = DataType_(0);
        vx01[k] = DataType_(3);
        vy01[i] = DataType_(3);
        vx1[k] = DataType_(1);
        vy1[i] = DataType_(1);
      }
    }

    // create two buffer-vectors
    BufferVector buf_x = mirror_x.create_buffer(vec_x);
    BufferVector buf_y = mirror_y.create_buffer(vec_y);

    // gather local vectors
    mirror_x.gather(buf_x, vec_x);
    mirror_y.gather(buf_y, vec_y);

    // scatter exchanged buffers
    mirror_x.scatter_axpy(vec_x, buf_y);
    mirror_y.scatter_axpy(vec_y, buf_x);

    // check against reference
    TEST_CHECK_LESS_THAN(vec_x.max_rel_diff(sync_x), tol);
    TEST_CHECK_LESS_THAN(vec_y.max_rel_diff(sync_y), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMirrorTest, double, std::uint32_t, PreferredBackend::cuda);
#endif

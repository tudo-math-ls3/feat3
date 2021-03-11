// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
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

template<typename MemType_, typename DataType_, typename IndexType_>
class MetaMirrorTest :
  public TestSystem::FullTaggedTest<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;

  typedef DenseVector<Mem::Main, DataType, IndexType> BufferVector;

  typedef DenseVector<MemType_, DataType, IndexType> ScalarVector;
  typedef PowerVector<ScalarVector, 2> PowerVector2;
  typedef TupleVector<PowerVector2, ScalarVector> MetaVector;

  typedef VectorMirror<MemType_, DataType, IndexType> ScalarMirror;
  typedef PowerMirror<ScalarMirror, 2> PowerMirror2;
  typedef TupleMirror<PowerMirror2, ScalarMirror> MetaMirror;

  MetaMirrorTest() :
    TestSystem::FullTaggedTest<MemType_, DataType_, IndexType_>("MetaMirrorTest")
  {
  }

  virtual ~MetaMirrorTest()
  {
  }

  static ScalarMirror gen_mirror_x(IndexType m)
  {
    ScalarMirror mir(m*m, m);
    IndexType* ci(mir.indices());
    IndexType m2(m*(m-1));
    for(IndexType i(0); i < m; ++i)
      ci[i] = m2+i;
    return mir;
  }

  static ScalarMirror gen_mirror_y(IndexType m)
  {
    ScalarMirror mir(m*m, m);
    IndexType* ci(mir.indices());
    for(IndexType i(0); i < m; ++i)
      ci[i] = i;
    return mir;
  }

  virtual void run() const override
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

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
    for(Index i(0); i < m; ++i)
    {
      Index k(n-m+i);
      sync_x.template at<0>().template at<0>()(k, DataType_(0));
      sync_y.template at<0>().template at<0>()(i, DataType_(0));
      sync_x.template at<0>().template at<1>()(k, DataType_(3));
      sync_y.template at<0>().template at<1>()(i, DataType_(3));
      sync_x.template at<1>()(k, DataType_(1));
      sync_y.template at<1>()(i, DataType_(1));
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

    // compute difference to reference
    vec_x.axpy(sync_x, vec_x, -DataType(1));
    vec_y.axpy(sync_y, vec_y, -DataType(1));

    // check difference norm
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_x.norm2(), DataType_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_y.norm2(), DataType_(0), tol);
  }
};

MetaMirrorTest<Mem::Main, float, Index> meta_mirror_test_generic_float_index;
MetaMirrorTest<Mem::Main, double, Index> meta_mirror_test_generic_double_index;

#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/power_mirror.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/transposition.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<typename AlgoType_, typename DataType_>
class MetaMirrorTest :
  public TestSystem::TaggedTest<typename AlgoType_::MemType, DataType_>
{
public:
  typedef AlgoType_ AlgoType;
  typedef typename AlgoType::MemType MemType;
  typedef DataType_ DataType;
  typedef DenseVector<MemType, DataType> ScalarVector;
  typedef PowerVector<ScalarVector, 2> PowerVector2;
  typedef TupleVector<PowerVector2, ScalarVector> MetaVector;

  typedef VectorMirror<MemType, DataType> ScalarMirror;
  typedef PowerMirror<ScalarMirror, 2> PowerMirror2;
  typedef TupleMirror<PowerMirror2, ScalarMirror> MetaMirror;

  typedef SparseMatrixCSR<MemType, DataType> ScalarMatrix;

  MetaMirrorTest() :
    TestSystem::TaggedTest<MemType, DataType>("MetaMirrorTest")
  {
  }

  static ScalarMatrix gen_mir_x(Index m)
  {
    DenseVector<Mem::Main, Index> row_ptr(m+1);
    DenseVector<Mem::Main, Index> col_idx(m);
    DenseVector<Mem::Main, DataType> data(m, DataType(1));

    Index* rp(row_ptr.elements());
    for(Index i(0); i <= m; ++i)
      rp[i] = i;

    Index* ci(col_idx.elements());
    Index m2(m*(m-1));
    for(Index i(0); i < m; ++i)
      ci[i] = m2+i;

    return SparseMatrixCSR<Mem::Main, DataType>(m, m*m, col_idx, data, row_ptr);
  }

  static ScalarMatrix gen_mir_y(Index m)
  {
    DenseVector<Mem::Main, Index> row_ptr(m+1);
    DenseVector<Mem::Main, Index> col_idx(m);
    DenseVector<Mem::Main, DataType> data(m, DataType(1));

    Index* rp(row_ptr.elements());
    for(Index i(0); i <= m; ++i)
      rp[i] = i;

    Index* ci(col_idx.elements());
    for(Index i(0); i < m; ++i)
      ci[i] = i;

    return SparseMatrixCSR<Mem::Main, DataType>(m, m*m, col_idx, data, row_ptr);
  }

  virtual void run() const
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.7)));

    const Index m = 3;
    const Index n = m*m;

    // create the two mirror matrices
    ScalarMatrix gather_x(gen_mir_x(m));
    ScalarMatrix gather_y(gen_mir_y(m));

    // transpose for scatter
    ScalarMatrix scatter_x(Transposition<AlgoType>::value(gather_x));
    ScalarMatrix scatter_y(Transposition<AlgoType>::value(gather_y));

    // create the meta-mirrors
    MetaMirror mirror_x(PowerMirror2(ScalarMirror(gather_x.clone(), scatter_x.clone())), ScalarMirror(gather_x.clone(), scatter_x.clone()));
    MetaMirror mirror_y(PowerMirror2(ScalarMirror(gather_y.clone(), scatter_y.clone())), ScalarMirror(gather_y.clone(), scatter_y.clone()));

    // create meta-vectors
    MetaVector vec_x;
    vec_x.template at<Index(0)>().template at<Index(0)>() = ScalarVector(n, DataType(1));
    vec_x.template at<Index(0)>().template at<Index(1)>() = ScalarVector(n, DataType(2));
    vec_x.template at<Index(1)>() = ScalarVector(n, DataType(3));
    MetaVector vec_y;
    vec_y.template at<Index(0)>().template at<Index(0)>() = ScalarVector(n, DataType(-1));
    vec_y.template at<Index(0)>().template at<Index(1)>() = ScalarVector(n, DataType( 1));
    vec_y.template at<Index(1)>() = ScalarVector(n, DataType(-2));

    // create reference synced vectors
    MetaVector sync_x(vec_x.clone());
    MetaVector sync_y(vec_y.clone());
    for(Index i(0); i < m; ++i)
    {
      Index k(n-m+i);
      sync_x.template at<Index(0)>().template at<Index(0)>()(k, DataType_(0));
      sync_y.template at<Index(0)>().template at<Index(0)>()(i, DataType_(0));
      sync_x.template at<Index(0)>().template at<Index(1)>()(k, DataType_(3));
      sync_y.template at<Index(0)>().template at<Index(1)>()(i, DataType_(3));
      sync_x.template at<Index(1)>()(k, DataType_(1));
      sync_y.template at<Index(1)>()(i, DataType_(1));
    }

    // create two buffer-vectors
    ScalarVector buf_x(mirror_x.size(), DataType(0));
    ScalarVector buf_y(mirror_y.size(), DataType(0));

    // gather local vectors
    mirror_x.gather_prim(buf_x, vec_x);
    mirror_y.gather_prim(buf_y, vec_y);

    // combine buffers: x,y <- x+y
    buf_x.template axpy<AlgoType>(buf_y, buf_x);
    buf_y.copy(buf_x);

    // scatter synced buffers
    mirror_x.scatter_prim(vec_x, buf_x);
    mirror_y.scatter_prim(vec_y, buf_y);

    // compute difference to reference
    vec_x.template axpy<AlgoType>(sync_x, vec_x, -DataType(1));
    vec_y.template axpy<AlgoType>(sync_y, vec_y, -DataType(1));

    // check difference norm
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_x.template norm2<AlgoType>(), DataType_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_y.template norm2<AlgoType>(), DataType_(0), tol);
  }
};

MetaMirrorTest<Algo::Generic, float> meta_mirror_test_generic_float;
MetaMirrorTest<Algo::Generic, double> meta_mirror_test_generic_double;

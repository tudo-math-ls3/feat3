#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template<typename DT_, typename IT_>
class MatrixMultTest
  : public FullTaggedTest<Mem::Main, DT_, IT_>
{
  typedef SparseMatrixCSR<Mem::Main, DT_, IT_> MatrixType;
  typedef VectorMirror<Mem::Main, DT_, IT_> VectorMirrorType;
  typedef MatrixMirror<Mem::Main, DT_, IT_> MatrixMirrorType;
  typedef MatrixMirrorBuffer<Mem::Main, DT_, IT_> MatMirBufType;

public:
  MatrixMultTest()
    : FullTaggedTest<Mem::Main, DT_, IT_>("MatrixMultTest")
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.6));

    // create 3 matrices a,b and d
    LAFEM::PointstarFactoryFD<DT_, IT_> psf(17, 2);
    MatrixType a = psf.matrix_csr();
    MatrixType b = a.clone();
    MatrixType d = a.clone();

    // fill a, b and d with random values
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    for(Index i(0); i < a.used_elements(); ++i)
    {
      a.val()[i] = DT_(rng(0.0, 1.0));
      b.val()[i] = DT_(rng(0.0, 1.0));
      d.val()[i] = DT_(rng(0.0, 1.0));
    }

    // create matrix structure for X = D*A*B
    Adjacency::Graph graph_da(Adjacency::rt_injectify, d, a);
    Adjacency::Graph graph_dab(Adjacency::rt_injectify, graph_da, b);
    graph_dab.sort_indices();
    MatrixType x(graph_dab);

    // create vector mirrors
    VectorMirrorType mir_b(b.transpose(), b.clone());
    VectorMirrorType mir_d(d.clone(), d.transpose());
    MatrixMirrorType mirror(mir_d, mir_b);

    // compute X := D*A*B using matrix mirror
    // This is a little tricky, as the matrix mirror works on MatrixMirrorBuffer
    // objects...
    MatMirBufType temp_c(Adjacency::Graph(Adjacency::rt_as_is, x), Index(1));
    mirror.gather(temp_c, a);
    const DT_* vt = temp_c.val();
    DT_* vx = x.val();
    for(Index i(0); i < x.used_elements(); ++i)
        vx[i] = vt[i];

    // subtract: X -= D*A*B
    x.add_double_mat_mult(d, a, b, -DT_(1));

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(x.norm_frobenius(), DT_(0), tol);
  }
};

MatrixMultTest<float, unsigned int> matrix_mult_test_float_uint;
MatrixMultTest<double, unsigned int> matrix_mult_test_double_uint;
MatrixMultTest<float, unsigned long> matrix_mult_test_float_ulong;
MatrixMultTest<double, unsigned long> matrix_mult_test_double_ulong;

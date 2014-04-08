#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<typename DataType_>
class PointstarFactoryTest :
  public FEAST::TestSystem::TaggedTest<Mem::Main, DataType_>
{
public:
  typedef DenseVector<Mem::Main, DataType_> VectorType;
  typedef SparseMatrixCSR<Mem::Main, DataType_> MatrixType;

public:
  PointstarFactoryTest() :
    FEAST::TestSystem::TaggedTest<Mem::Main, DataType_>("PointstarFactoryTest")
  {
  }

  virtual void run() const
  {
    const DataType_ tol(Math::pow(Math::eps<DataType_>(), DataType_(0.8)));

    // Test FD poinstars
    // dimension loop: d=1,2,3,4
    for(Index d(1); d < 5; ++d)
    {
      // point loop: m=2,3,4,5
      for(Index m(3); m < 6; ++m)
      {
        // generate FD matrix A
        PointstarFactoryFD<DataType_> factory(m, d);
        MatrixType a(factory.matrix_csr());

        // compute smallest and largest eigenvalues of A
        const DataType_ lambda_min(factory.lambda_min());
        //const DataType_ lambda_max(factory.lambda_max());

        // generate eigenvector
        const VectorType ev(factory.eigenvector_min());

        // compute w = A*ev - lambda_min*ev
        VectorType w(ev.size());
        a.template apply<Algo::Generic>(w, ev);
        w.template axpy<Algo::Generic>(ev, w, -lambda_min);

        // check norm of w
        TEST_CHECK_EQUAL_WITHIN_EPS(w.template norm2<Algo::Generic>(), DataType_(0), tol);
      }
    }

    // Test FE poinstars
    for(Index m(3); m < 9; ++m)
    {
      // generate 2D FE matrix A
      PointstarFactoryFE<DataType_> factory(m);
      MatrixType a(factory.matrix_csr());

      // compute smallest and largest eigenvalues of A
      const DataType_ lambda_min(factory.lambda_min());
      //const DataType_ lambda_max(factory.lambda_max());

      // generate eigenvector
      const VectorType ev(factory.eigenvector_min());

      // compute w = A*ev - lambda_min*ev
      VectorType w(ev.size());
      a.template apply<Algo::Generic>(w, ev);
      w.template axpy<Algo::Generic>(ev, w, -lambda_min);

      // check norm of w
      TEST_CHECK_EQUAL_WITHIN_EPS(w.template norm2<Algo::Generic>(), DataType_(0), tol);
      //std::cout << m << ": " << scientify(lambda_min,10) << " , " << scientify(w(0)/ev(0),10) << std::endl;
    }
  }
};

PointstarFactoryTest<double> pointstar_factory_test_double;









template<typename DataType_>
class PointstarStructureTest :
  public FEAST::TestSystem::TaggedTest<Mem::Main, DataType_>
{
public:
  typedef DenseVector<Mem::Main, DataType_> VectorType;
  typedef SparseMatrixBanded<Mem::Main, DataType_> MatrixType;

public:
  PointstarStructureTest() :
    FEAST::TestSystem::TaggedTest<Mem::Main, DataType_>("PointstarStructureTest")
  {
  }

  virtual void run() const
  {
    const DataType_ tol(Math::pow(Math::eps<DataType_>(), DataType_(0.8)));

    DenseVector<Mem::Main, Index> num_of_nodes(4);
    DenseVector<Mem::Main, DataType_> dimensions(3);

    num_of_nodes(0, 2);
    num_of_nodes(1, 9);
    num_of_nodes(2, 11);
    num_of_nodes(3, 24);

    dimensions(0, DataType_(3.0));
    dimensions(1, DataType_(0.47));
    dimensions(2, DataType_(4.0));

    // generate FD matrix A
    PointstarFactoryFD2<DataType_> factory(num_of_nodes, dimensions);
    MatrixType a(factory.matrix_banded());

    // compute smallest and largest eigenvalues of A
    const DataType_ lambda_min(factory.lambda_min());

    // generate eigenvector
    const VectorType ev(factory.eigenvector_min());

    // compute w = A*ev - lambda_min*ev
    VectorType w(ev.size());
    a.template apply<Algo::Generic>(w, ev);
    w.template axpy<Algo::Generic>(ev, w, -lambda_min);

    // check norm of w
    TEST_CHECK_EQUAL_WITHIN_EPS(w.template norm2<Algo::Generic>(), DataType_(0), tol);


    DenseVector<Mem::Main, Index> num_of_nodes2(3);
    num_of_nodes2(0, 5);
    num_of_nodes2(1, 4);
    num_of_nodes2(2, 3);

    MatrixType b(PointstarStructureFE<Algo::Generic>::value<DataType_>(3, num_of_nodes2));
  }
};

PointstarStructureTest<double> pointstar_structure_test_double;

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template<typename DataType_, typename IndexType_>
class PointstarFactoryTest :
  public FEAT::TestSystem::FullTaggedTest<Mem::Main, DataType_, IndexType_>
{
public:
  typedef DenseVector<Mem::Main, DataType_, IndexType_> VectorType;
  typedef SparseMatrixCSR<Mem::Main, DataType_, IndexType_> MatrixType;

public:
   PointstarFactoryTest() :
    FEAT::TestSystem::FullTaggedTest<Mem::Main, DataType_, IndexType_>("PointstarFactoryTest")
  {
  }

  virtual ~PointstarFactoryTest()
  {
  }

  virtual void run() const override
  {
    const DataType_ tol(Math::pow(Math::eps<DataType_>(), DataType_(0.7)));

    // Test FD poinstars
    // dimension loop: d=1,2,3,4
    for(Index d(1); d < 5; ++d)
    {
      // point loop: m=2,3,4,5
      for(Index m(3); m < 6; ++m)
      {
        // generate FD matrix A
        PointstarFactoryFD<DataType_, IndexType_> factory(m, d);
        MatrixType a(factory.matrix_csr());

        // compute smallest and largest eigenvalues of A
        const DataType_ lambda_min(factory.lambda_min());
        //const DataType_ lambda_max(factory.lambda_max());

        // generate eigenvector
        const VectorType ev(factory.eigenvector_min());

        // compute w = A*ev - lambda_min*ev
        VectorType w(ev.size());
        a.apply(w, ev);
        w.axpy(ev, w, -lambda_min);

        // check norm of w
        TEST_CHECK_EQUAL_WITHIN_EPS(w.norm2(), DataType_(0), tol);
      }
    }

    // Test FE poinstars
    for(Index m(3); m < 9; ++m)
    {
      // generate 2D FE matrix A
      PointstarFactoryFE<DataType_, IndexType_> factory(m);
      MatrixType a(factory.matrix_csr());

      // compute smallest and largest eigenvalues of A
      const DataType_ lambda_min(factory.lambda_min());
      //const DataType_ lambda_max(factory.lambda_max());

      // generate eigenvector
      const VectorType ev(factory.eigenvector_min());

      // compute w = A*ev - lambda_min*ev
      VectorType w(ev.size());
      a.apply(w, ev);
      w.axpy(ev, w, -lambda_min);

      // check norm of w
      TEST_CHECK_EQUAL_WITHIN_EPS(w.norm2(), DataType_(0), tol);
      //std::cout << m << ": " << stringify_fp_sci(lambda_min,10) << " , " << stringify_fp_sci(w(0)/ev(0),10) << std::endl;
    }
  }
};

PointstarFactoryTest<float, unsigned int> pointstar_factory_test_float_uint;
PointstarFactoryTest<float, unsigned long> pointstar_factory_test_float_ulong;
PointstarFactoryTest<double, unsigned int> pointstar_factory_test_double_uint;
PointstarFactoryTest<double, unsigned long> pointstar_factory_test_double_ulong;









template<typename DataType_, typename IndexType_>
class PointstarStructureTest :
  public FEAT::TestSystem::FullTaggedTest<Mem::Main, DataType_, IndexType_>
{
public:
  typedef DenseVector<Mem::Main, DataType_, IndexType_> VectorType;
  typedef SparseMatrixBanded<Mem::Main, DataType_, IndexType_> MatrixType;

public:
   PointstarStructureTest() :
    FEAT::TestSystem::FullTaggedTest<Mem::Main, DataType_, IndexType_>("PointstarStructureTest")
  {
  }

  virtual ~PointstarStructureTest()
  {
  }

  virtual void run() const override
  {
    const DataType_ tol(Math::pow(Math::eps<DataType_>(), DataType_(0.8)));


    std::vector<IndexType_> num_of_subintervalls;
    num_of_subintervalls.push_back(2);
    num_of_subintervalls.push_back(9);
    num_of_subintervalls.push_back(11);
    num_of_subintervalls.push_back(24);

    std::vector<DataType_> dimensions;
    dimensions.push_back(DataType_(3.0));
    dimensions.push_back(DataType_(0.47));
    dimensions.push_back(DataType_(4.0));

    // generate FD matrix A
    PointstarFactoryFD2<DataType_, IndexType_> factory(num_of_subintervalls, dimensions);
    MatrixType a(factory.matrix_banded());

    // compute smallest and largest eigenvalues of A
    const DataType_ lambda_min(factory.lambda_min());

    // generate eigenvector
    const VectorType ev(factory.eigenvector_min());

    // compute w = A*ev - lambda_min*ev
    VectorType w(ev.size());
    a.apply(w, ev);
    w.axpy(ev, w, -lambda_min);

    // check norm of w
    TEST_CHECK_EQUAL_WITHIN_EPS(w.norm2(), DataType_(0), tol);


    std::vector<IndexType_> num_of_subintervalls2;
    num_of_subintervalls2.push_back(5);
    num_of_subintervalls2.push_back(4);
    num_of_subintervalls2.push_back(3);

    MatrixType b(PointstarStructureFE::value<DataType_>(3, num_of_subintervalls2));
  }
};

PointstarStructureTest<float, unsigned int> pointstar_structure_test_float_uint;
PointstarStructureTest<float, unsigned long> pointstar_structure_test_float_ulong;
PointstarStructureTest<double, unsigned int> pointstar_structure_test_double_uint;
PointstarStructureTest<double, unsigned long> pointstar_structure_test_double_ulong;

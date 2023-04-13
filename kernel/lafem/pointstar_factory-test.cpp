// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template<
 typename DT_,
 typename IT_>
class PointstarFactoryTest
  : public UnitTest
{
public:
  typedef DenseVector<DT_, IT_> VectorType;
  typedef SparseMatrixCSR<DT_, IT_> MatrixType;

public:
   PointstarFactoryTest(PreferredBackend backend)
   : UnitTest("PointstarFactoryTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~PointstarFactoryTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.7)));

    // Test FD poinstars
    // dimension loop: d=1,2,3,4
    for(Index d(1); d < 5; ++d)
    {
      // point loop: m=2,3,4,5
      for(Index m(3); m < 6; ++m)
      {
        // generate FD matrix A
        PointstarFactoryFD<DT_, IT_> factory(m, d);
        MatrixType a(factory.matrix_csr());

        // compute smallest and largest eigenvalues of A
        const DT_ lambda_min(factory.lambda_min());
        //const DataType_ lambda_max(factory.lambda_max());

        // generate eigenvector
        const VectorType ev(factory.eigenvector_min());

        // compute w = A*ev - lambda_min*ev
        VectorType w(ev.size());
        a.apply(w, ev);
        w.axpy(ev, w, -lambda_min);

        // check norm of w
        TEST_CHECK_EQUAL_WITHIN_EPS(w.norm2(), DT_(0), tol);
      }
    }

    // Test FE poinstars
    for(Index m(3); m < 9; ++m)
    {
      // generate 2D FE matrix A
      PointstarFactoryFE<DT_, IT_> factory(m);
      MatrixType a(factory.matrix_csr());

      // compute smallest and largest eigenvalues of A
      const DT_ lambda_min(factory.lambda_min());
      //const DataType_ lambda_max(factory.lambda_max());

      // generate eigenvector
      const VectorType ev(factory.eigenvector_min());

      // compute w = A*ev - lambda_min*ev
      VectorType w(ev.size());
      a.apply(w, ev);
      w.axpy(ev, w, -lambda_min);

      // check norm of w
      TEST_CHECK_EQUAL_WITHIN_EPS(w.norm2(), DT_(0), tol);
      //std::cout << m << ": " << stringify_fp_sci(lambda_min,10) << " , " << stringify_fp_sci(w(0)/ev(0),10) << std::endl;
    }
  }
};

PointstarFactoryTest<float, unsigned int> pointstar_factory_test_float_uint(PreferredBackend::generic);
PointstarFactoryTest<float, unsigned long> pointstar_factory_test_float_ulong(PreferredBackend::generic);
PointstarFactoryTest<double, unsigned int> pointstar_factory_test_double_uint(PreferredBackend::generic);
PointstarFactoryTest<double, unsigned long> pointstar_factory_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
PointstarFactoryTest<float, unsigned long> mkl_pointstar_factory_test_float_ulong(PreferredBackend::mkl);
PointstarFactoryTest<double, unsigned long> mkl_pointstar_factory_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
PointstarFactoryTest<__float128, unsigned long> pointstar_factory_test_float128_ulong(PreferredBackend::generic);
PointstarFactoryTest<__float128, unsigned int> pointstar_factory_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
PointstarFactoryTest<Half, unsigned int> pointstar_factory_test_half_uint(PreferredBackend::generic);
PointstarFactoryTest<Half, unsigned long> pointstar_factory_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
PointstarFactoryTest<float, unsigned long> cuda_pointstar_factory_test_float_ulong(PreferredBackend::cuda);
PointstarFactoryTest<double, unsigned long> cuda_pointstar_factory_test_double_ulong(PreferredBackend::cuda);
PointstarFactoryTest<float, unsigned int> cuda_pointstar_factory_test_float_uint(PreferredBackend::cuda);
PointstarFactoryTest<double, unsigned int> cuda_pointstar_factory_test_double_uint(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class PointstarStructureTest
  :  public UnitTest
{
public:
  typedef DenseVector<DT_, IT_> VectorType;
  typedef SparseMatrixBanded<DT_, IT_> MatrixType;

public:
   PointstarStructureTest(PreferredBackend backend)
   : UnitTest("PointstarStructureTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~PointstarStructureTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol(Math::pow(Math::eps<DT_>(), DT_(0.8)));


    std::vector<IT_> num_of_subintervalls;
    num_of_subintervalls.push_back(2);
    num_of_subintervalls.push_back(9);
    num_of_subintervalls.push_back(11);
    num_of_subintervalls.push_back(24);

    std::vector<DT_> dimensions;
    dimensions.push_back(DT_(3.0));
    dimensions.push_back(DT_(0.47));
    dimensions.push_back(DT_(4.0));

    // generate FD matrix A
    PointstarFactoryFD2<DT_, IT_> factory(num_of_subintervalls, dimensions);
    MatrixType a(factory.matrix_banded());

    // compute smallest and largest eigenvalues of A
    const DT_ lambda_min(factory.lambda_min());

    // generate eigenvector
    const VectorType ev(factory.eigenvector_min());

    // compute w = A*ev - lambda_min*ev
    VectorType w(ev.size());
    a.apply(w, ev);
    w.axpy(ev, w, -lambda_min);

    // check norm of w
    TEST_CHECK_EQUAL_WITHIN_EPS(w.norm2(), DT_(0), tol);


    std::vector<IT_> num_of_subintervalls2;
    num_of_subintervalls2.push_back(5);
    num_of_subintervalls2.push_back(4);
    num_of_subintervalls2.push_back(3);

    MatrixType b(PointstarStructureFE::value<DT_>(3, num_of_subintervalls2));
  }
};

PointstarStructureTest<float, unsigned int> pointstar_structure_test_float_uint(PreferredBackend::generic);
PointstarStructureTest<float, unsigned long> pointstar_structure_test_float_ulong(PreferredBackend::generic);
PointstarStructureTest<double, unsigned int> pointstar_structure_test_double_uint(PreferredBackend::generic);
PointstarStructureTest<double, unsigned long> pointstar_structure_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
PointstarStructureTest<float, unsigned long> mkl_pointstar_structure_test_float_ulong(PreferredBackend::mkl);
PointstarStructureTest<double, unsigned long> mkl_pointstar_structure_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
PointstarStructureTest<__float128, unsigned long> pointstar_structure_test_float128_ulong(PreferredBackend::generic);
PointstarStructureTest<__float128, unsigned int> pointstar_structure_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
PointstarStructureTest<Half, unsigned int> pointstar_structure_test_half_uint(PreferredBackend::generic);
PointstarStructureTest<Half, unsigned long> pointstar_structure_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
PointstarStructureTest<float, unsigned long> cuda_pointstar_structure_test_float_ulong(PreferredBackend::cuda);
PointstarStructureTest<double, unsigned long> cuda_pointstar_structure_test_double_ulong(PreferredBackend::cuda);
PointstarStructureTest<float, unsigned int> cuda_pointstar_structure_test_float_uint(PreferredBackend::cuda);
PointstarStructureTest<double, unsigned int> cuda_pointstar_structure_test_double_uint(PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
PointstarStructureTest<Half, unsigned int> cuda_pointstar_structure_test_half_uint(PreferredBackend::cuda);
PointstarStructureTest<Half, unsigned long> cuda_pointstar_structure_test_half_ulong(PreferredBackend::cuda);
#endif
#endif

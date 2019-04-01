// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Meta-Vector component-product test class
 *
 * \test The 'component_product' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<typename MemType_, typename DataType_, typename IndexType_>
class MetaVectorCompProdTest
  : public MetaVectorTestBase<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<MemType_, DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

   MetaVectorCompProdTest() :
    BaseClass("MetaVectorCompProdTest")
  {
  }

  virtual ~MetaVectorCompProdTest()
  {
  }

  using BaseClass::fx00;
  using BaseClass::fx01;
  using BaseClass::fx1;
  using BaseClass::fy00;
  using BaseClass::fy01;
  using BaseClass::fy1;

  virtual void run() const override
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    MetaVector x(this->gen_vector_x(n00, n01, n1));
    MetaVector y(this->gen_vector_y(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // test: z <- x * y
    // purpose: general test
    z.component_product(x, y);
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), fx00(i) * fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), fx01(i) * fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), fx1(i) * fy1(i), tol);

    // test: z <- x; z <- z * y
    // purpose: z = x
    z.copy(x);
    z.component_product(z, y);
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), fx00(i) * fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), fx01(i) * fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), fx1(i) * fy1(i), tol);

    // test: z <- y; z <- x * z
    // purpose: z = y
    z.copy(y);
    z.component_product(x, z);
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), fx00(i) * fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), fx01(i) * fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), fx1(i) * fy1(i), tol);
  }
};

MetaVectorCompProdTest<Mem::Main, float, Index> meta_vector_comp_prod_test_generic_float;
MetaVectorCompProdTest<Mem::Main, double, Index> meta_vector_comp_prod_test_generic_double;
#ifdef FEAT_HAVE_CUDA
MetaVectorCompProdTest<Mem::CUDA, float, Index> meta_vector_comp_prod_test_cuda_float;
MetaVectorCompProdTest<Mem::CUDA, double, Index> meta_vector_comp_prod_test_cuda_double;
#endif

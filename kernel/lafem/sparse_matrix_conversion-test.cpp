#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/util/random.hpp>

#include <cstdio>
#include <sstream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for sparse matrix conversions.
*
* \test test description missing
*
* \tparam Mem_
* description missing
*
* \tparam DT_
* description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Mem_,
  typename DT_>
class SparseMatrixConversionTest
  : public TaggedTest<Mem_, DT_>
{

public:

  SparseMatrixConversionTest()
    : TaggedTest<Mem_, DT_>("sparse_matrix_conversion_test")
  {
  }

  virtual void run() const
  {

    SparseMatrixCOO<Mem::Main, DT_> a(121, 121);
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        if(row == col)
          a(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          a(row, col, DT_(-1));
      }
    }

    {
      SparseMatrixCOO<Mem::Main, DT_> coo_m1;
      coo_m1.convert(a);
      SparseMatrixELL<Mem::Main, DT_> ell_m1(coo_m1);
      SparseMatrixCOO<Mem::Main, DT_> coo_m2(ell_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m1(coo_m1);
      SparseMatrixCOO<Mem::Main, DT_> coo_m3(csr_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
    }

    {
      SparseMatrixELL<Mem::Main, DT_> ell_m1(a);
      SparseMatrixCOO<Mem::Main, DT_> coo_m1(ell_m1);
      SparseMatrixELL<Mem::Main, DT_> ell_m2(coo_m1);
      TEST_CHECK_EQUAL(ell_m2, ell_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m1(ell_m1);
      SparseMatrixELL<Mem::Main, DT_> ell_m3(csr_m1);
      TEST_CHECK_EQUAL(ell_m3, ell_m1);
    }

    {
      SparseMatrixCSR<Mem::Main, DT_> csr_m1(a);
      SparseMatrixCOO<Mem::Main, DT_> coo_m1(csr_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m2(coo_m1);
      TEST_CHECK_EQUAL(csr_m2, csr_m1);
      SparseMatrixELL<Mem::Main, DT_> ell_m1(csr_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m3(ell_m1);
      TEST_CHECK_EQUAL(csr_m3, csr_m1);
    }

  }
};
SparseMatrixConversionTest<Mem::Main, float> sparse_matrix_conversion_test_float;
SparseMatrixConversionTest<Mem::Main, double> sparse_matrix_conversion_test_double;

/**
* \brief Test class for sparse matrix conversions.
*
* \test test description missing
*
* \tparam Mem_
* description missing
*
* \tparam DT_
* description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Mem_,
  typename DT_>
class SparseMatrixCudaConversionTest
  : public TaggedTest<Mem_, DT_>
{

public:

  SparseMatrixCudaConversionTest()
    : TaggedTest<Mem_, DT_>("sparse_matrix_cuda_conversion_test")
  {
  }

  virtual void run() const
  {

    SparseMatrixCOO<Mem::Main, DT_> a(121, 121);
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        if(row == col)
          a(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          a(row, col, DT_(-1));
      }
    }

    {
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1;
      coo_m1.convert(a);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m2(ell_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m3(csr_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
    }

    {
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(a);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m2(coo_m1);
      TEST_CHECK_EQUAL(ell_m2, ell_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m3(csr_m1);
      TEST_CHECK_EQUAL(ell_m3, ell_m1);
    }

    {
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(a);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m2(coo_m1);
      TEST_CHECK_EQUAL(csr_m2, csr_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m3(ell_m1);
      TEST_CHECK_EQUAL(csr_m3, csr_m1);
    }

    {
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1;
      coo_m1.convert(a);
      SparseMatrixELL<Mem::Main, DT_> ell_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m2(ell_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m3(csr_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
    }

    {
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(a);
      SparseMatrixCOO<Mem::Main, DT_> coo_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m2(coo_m1);
      TEST_CHECK_EQUAL(ell_m2, ell_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m3(csr_m1);
      TEST_CHECK_EQUAL(ell_m3, ell_m1);
    }

    {
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(a);
      SparseMatrixCOO<Mem::Main, DT_> coo_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m2(coo_m1);
      TEST_CHECK_EQUAL(csr_m2, csr_m1);
      SparseMatrixELL<Mem::Main, DT_> ell_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m3(ell_m1);
      TEST_CHECK_EQUAL(csr_m3, csr_m1);
    }

  }
};
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCudaConversionTest<Mem::CUDA, float> sparse_matrix_cuda_conversion_test_float;
SparseMatrixCudaConversionTest<Mem::CUDA, double> sparse_matrix_cuda_conversion_test_double;
#endif

/**
 * \brief Test class for the convertion from SparseMatrixBanded to another matrix-format
 *
 * \test test description missing
 *
 * \tparam MT_
 * description missing
 *
 * \author Christoph Lohmann
 */
template<typename MT_>
class SparseMatrixBandedConversionTest
  : public TaggedTest<typename MT_::MemType, typename MT_::DataType>
{
public:
  SparseMatrixBandedConversionTest()
    : TaggedTest<typename MT_::MemType, typename MT_::DataType>("sparse_matrix_banded_conversion_test: " + MT_::name())
  {
  }

  typedef Algo::Generic Algo_;
  typedef typename MT_::MemType Mem_;
  typedef typename MT_::DataType DT_;

  virtual void run() const
  {
    Random random;
    const Index size(9);

    DenseVector<Mem_, Index, Index> offsets(4);
    DenseVector<Mem_, DT_, Index> val(offsets.size() * size);

    offsets(0, 3);
    offsets(1, 4);
    offsets(2, 9);
    offsets(3, 12);

    for (Index i(0); i < val.size(); ++i)
    {
      val(i, random(DT_(0), DT_(10)));
    }

    SparseMatrixBanded<Mem_, DT_, Index> sys_banded(size, size + 1, val, offsets);
    SparseMatrixCSR<Mem_, DT_, Index> sys_csr(sys_banded);

    DenseVector<Mem_, DT_, Index> x_banded(sys_banded.columns(), DT_(0.0));
    DenseVector<Mem_, DT_, Index> y_banded(sys_banded.rows());
    DenseVector<Mem_, DT_, Index> x_csr(sys_csr.columns(), DT_(0.0));
    DenseVector<Mem_, DT_, Index> y_csr(sys_csr.rows());

    for (Index i(0); i < sys_banded.columns(); ++i)
    {
      x_banded(i, DT_(1.0));
      x_csr(i, DT_(1.0));

      sys_banded.template apply<Algo_>(y_banded, x_banded);
      sys_csr.template apply<Algo_>(y_csr, x_csr);

      // check, if the result is correct
      for (Index j(0) ; j < sys_banded.rows() ; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(y_banded(j), y_csr(j), 1e-8);
      }

      x_banded(i, DT_(0.0));
      x_csr(i, DT_(0.0));
    }
  }
};

SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::Main, float> > sparse_matrix_banded_csr_conversion_test_float;
SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::Main, double> > sparse_matrix_banded_csr_conversion_test_double;

SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::Main, float> > sparse_matrix_banded_ell_conversion_test_float;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::Main, double> > sparse_matrix_banded_ell_conversion_test_double;

SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::Main, float> > sparse_matrix_banded_coo_conversion_test_float;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::Main, double> > sparse_matrix_banded_coo_conversion_test_double;

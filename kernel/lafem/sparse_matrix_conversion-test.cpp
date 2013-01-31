#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

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
    for (unsigned long row(0) ; row < a.rows() ; ++row)
    {
      for (unsigned long col(0) ; col < a.columns() ; ++col)
      {
        if(row == col)
          a(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          a(row, col, DT_(-1));
      }
    }

    {
      SparseMatrixCOO<Mem::Main, DT_> coo_m1(a);
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
    for (unsigned long row(0) ; row < a.rows() ; ++row)
    {
      for (unsigned long col(0) ; col < a.columns() ; ++col)
      {
        if(row == col)
          a(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          a(row, col, DT_(-1));
      }
    }

    {
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1(a);
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
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1(a);
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

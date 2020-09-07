// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;

/**
 * \brief Test class for the sparse matrix factory class
 *
 *\author Gesa Pottbrock
 *
 */
template<typename Mem_, typename DT_, typename IT_>
class SparseMatrixFactoryTest :
  public TestSystem::FullTaggedTest<Mem_, DT_, IT_>
{
public:
  SparseMatrixFactoryTest() :
    TestSystem::FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixFactoryTest")
  {
  }

  virtual void run() const override
  {
    /* 8*7 Matrix getting implemented in this test:
    13,5   - -    - -     - 200
    -      3 -    - -     - -
    -      - -    - 21.13 - -
    -      - -    - -     - -
    -      - -    - -     - -
    -      - -    - -     - -
    -      - 10.5 - -     - -
    1034.5 - -    - -     - 18.7

    */
    //Implement 8*7 sparsematrix
    SparseMatrixFactory<DT_, IT_> factory(IT_(8), IT_(7));
    //Adding Elements to our Factory
    factory.add(IT_(0), IT_(0), DT_(13.5));
    factory.add(IT_(7), IT_(6), DT_(18.7));
    factory.add(IT_(6), IT_(2), DT_(10.5));
    factory.add(IT_(2), IT_(4), DT_(21.13));
    factory.add(IT_(0), IT_(6), DT_(200));
    factory.add(IT_(1), IT_(1), DT_(3));
    factory.add(IT_(7), IT_(0), DT_(1034.5));
    //Converting  Factory to CSR Matrix
    SparseMatrixCSR<Mem_, DT_, IT_> matrix_csr(factory.make_csr());

    //Testing if CSR Matrixs has the correct dimension and NNZ
    TEST_CHECK_EQUAL(matrix_csr.rows(),8);
    TEST_CHECK_EQUAL(factory.rows(), 8);
    TEST_CHECK_EQUAL(matrix_csr.columns(), 7);
    TEST_CHECK_EQUAL(factory.columns(), 7);
    TEST_CHECK_EQUAL(factory.size(), IT_(7 * 8));
    TEST_CHECK_EQUAL(matrix_csr.used_elements(), 7);
    TEST_CHECK_EQUAL(factory.used_elements(), IT_(7));

    const IT_* row_ptr = matrix_csr.row_ptr();
    const IT_* col_ind = matrix_csr.col_ind();
    const DT_* val = matrix_csr.val();

    //Testing if the value array is implemented in the correkt order
    TEST_CHECK_EQUAL(val[0], DT_(13.5));
    TEST_CHECK_EQUAL(val[1], DT_(200));
    TEST_CHECK_EQUAL(val[2], DT_(3));
    TEST_CHECK_EQUAL(val[3], DT_(21.13));
    TEST_CHECK_EQUAL(val[5], DT_(1034.5));
    TEST_CHECK_EQUAL(val[6], DT_(18.7));

    //Testing if the column indize array is correctly implemented
    TEST_CHECK_EQUAL(col_ind[0], IT_(0));
    TEST_CHECK_EQUAL(col_ind[1], IT_(6));
    TEST_CHECK_EQUAL(col_ind[2], IT_(1));
    TEST_CHECK_EQUAL(col_ind[3], IT_(4));
    TEST_CHECK_EQUAL(col_ind[4], IT_(2));
    TEST_CHECK_EQUAL(col_ind[5], IT_(0));
    TEST_CHECK_EQUAL(col_ind[6], IT_(6));

    //Testing if the row pointer array is implemented correctly
    TEST_CHECK_EQUAL(row_ptr[0], IT_(0));
    TEST_CHECK_EQUAL(row_ptr[1], IT_(2));
    TEST_CHECK_EQUAL(row_ptr[2], IT_(3));
    TEST_CHECK_EQUAL(row_ptr[3], IT_(4));
    TEST_CHECK_EQUAL(row_ptr[4], IT_(4));
    TEST_CHECK_EQUAL(row_ptr[5], IT_(4));
    TEST_CHECK_EQUAL(row_ptr[6], IT_(4));
    TEST_CHECK_EQUAL(row_ptr[7], IT_(5));
    TEST_CHECK_EQUAL(row_ptr[8], IT_(7));

    //Testing
   /* TEST_CHECK_EQUAL(matrix_csr(0,0), DT_(13.5));
    TEST_CHECK_EQUAL(matrix_csr(7,6), DT_(18.7));
    TEST_CHECK_EQUAL(matrix_csr(6,2), DT_(10.5));
    TEST_CHECK_EQUAL(matrix_csr(2,4 ), DT_(21.13));
    TEST_CHECK_EQUAL(matrix_csr(0, 6), DT_(200));
    TEST_CHECK_EQUAL(matrix_csr(1,1 ), DT_(3));
    TEST_CHECK_EQUAL(matrix_csr(7, 0), DT_(1034.5));
    */

  }
}; // class SparseMatrixFactoryTest<...>

SparseMatrixFactoryTest<Mem::Main, double, Index> sparse_matrix_factory_test_main_double_index;

// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse matrix factory class
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 *\author Gesa Pottbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseMatrixFactoryTest
  : public UnitTest
{
public:
  explicit SparseMatrixFactoryTest(PreferredBackend backend)
    : UnitTest("SparseMatrixFactoryTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
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
    SparseMatrixCSR<DT_, IT_> matrix_csr(factory.make_csr());

    //Testing if CSR Matrixs has the correct dimension and NNZ
    TEST_CHECK_EQUAL(matrix_csr.num_rows(),8);
    TEST_CHECK_EQUAL(factory.num_rows(), 8);
    TEST_CHECK_EQUAL(matrix_csr.num_cols(), 7);
    TEST_CHECK_EQUAL(factory.num_cols(), 7);
    TEST_CHECK_EQUAL(matrix_csr.num_nzes(), 7);
    TEST_CHECK_EQUAL(factory.num_nzes(), IT_(7));

    const Memory::TypedView<IT_> row_ptr = matrix_csr.row_ptr_view_r();
    const Memory::TypedView<IT_> col_ind = matrix_csr.col_idx_view_r();
    const Memory::TypedView<DT_> val = matrix_csr.val_view_r();

    //Testing if the value array is implemented in the correct order
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

SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseMatrixFactoryTest, double, std::uint32_t, PreferredBackend::cuda);
#endif
